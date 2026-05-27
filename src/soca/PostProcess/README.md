# PostProcessIce: Sea-Ice Analysis Postprocessing for CICE Restarts

## Overview

`PostProcessIce` projects a SOCA sea-ice analysis (aggregated, single-category
`sea_ice_area_fraction` / `sea_ice_thickness` / `sea_ice_snow_thickness`) onto a
CICE per-category restart. It is the C++/atlas successor to the
[`Soca2Cice`](../VariableChange/Soca2Cice/README.md) Fortran variable change.

CICE represents sea ice with `ncat` thickness categories, each carrying its own
`aicen` (area fraction), `vicen` (volume), `vsnon` (snow volume), surface
temperature `Tsfcn`, per-layer enthalpy `qice00N`, per-layer salinity `sice00N`,
per-layer snow enthalpy `qsno00N`, and melt-pond fields (`apnd`, `hpnd`,
`ipnd`). SOCA's analysis only constrains the aggregates. PostProcessIce
bridges the two while preserving the per-cat structure and keeping the CICE
restart thermodynamically self-consistent enough that CICE will run from it
without crashing on the first timestep.

## Workflow

A single call:

```cpp
PostProcessIce ppIce(geom, config.getSubConfiguration("postprocess ice"));
soca::State aggregateOut = ppIce.postprocess(analysis);
```

does three things:

1. **Read** the per-category CICE background restart at
   `postprocess ice: cice restart: input:` into a `soca::State`. The
   per-cat/per-(cat,layer) variable list (~115 names) is auto-injected from
   `ncat / ice_lev / sno_lev`; the caller does not enumerate them.
2. **Process** the per-cell pass below, mutating `pproc` (initialised from
   `restart`) toward `analysis`.
3. **Write** the postprocessed restart at `postprocess ice: cice restart:
   output:` in **update mode**: the input file is the template, the writer
   byte-copies it to the output path and overwrites only the variables SOCA
   models. Roughly 40 CICE restart variables PostProcessIce does not touch
   (dynamics state, tracers it doesn't analyse, etc.) pass through unchanged.

The returned `aggregateOut` State carries the `sea_ice_area_fraction`,
`sea_ice_thickness`, and `sea_ice_snow_thickness` fields *actually written* to
the restart (post min-vice cleanup, post freeboard, etc.). It is on
`analysis`'s geometry. Standalone callers normally discard it; the
`gdasapp` increment handler can use it to recover the effective ice
analysis if needed.

PostProcessIce does **not** read or write ocean variables: ocean T/S on the
analysis is unused. Any surface-ocean adjustment on ice removal (e.g.
warming SST toward freezing when ice2noice) is the DA increment generator's
job, not the restart writer's.

## Per-cell pass

For each owned cell:

1. **Clamp** analysis bounds (`a_aice` to `[0, 1]`; `a_hice`, `a_hsno` to
   `[0, inf)`). Drop near-zero analysis aice cells below `min aice output`
   (prevents the rebin from leaving thick-ice-on-edge artifacts).

2. **Classify** by `(ai_bg, ai_an)`:
   - **LAND** (mask=0): always zero per-cat output.
   - **ICE2NOICE** (`ai_an ~ 0`, `ai_bg > 0`): zero per-cat output.
   - **NOICE2ICE** (`ai_bg ~ 0`, `ai_an > 0`): seed all area in cat 0; set
     volume target from analysis hice if positive, else from the KDTree
     donor mean hice, else from `min new ice thickness * ai_an`.
   - **ICE2ICE**: alpha-rescale per-cat aicen and vicen by `ai_an / ai_bg`;
     set volume target to `analysis.hice * ai_an`.

3. **ITD rebin** (`itd.rebin: true`, default): solve for `(aicen, vicen)`
   such that `aicen[k] == 0` or `hicat[k] <= vicen[k]/aicen[k] <= hicat[k+1]`
   for every category, while hitting both `Sum(aicen) = a_aice` and
   `Sum(vicen) = vtot_target`. See [`IcePhysics.h`](IcePhysics.h)
   `adjustThicknessCategories` for the algorithm.

4. **Snow distribution**: take the cell-mean snow thickness target from
   analysis hsno (or fall back to background cell-mean), clamp to
   `[min snow thickness, max snow thickness]`, and distribute per
   category proportional to `aicen` (so all surviving cats have the same
   per-area snow depth, matching the standard SOCA convention).

5. **Freeboard enforcement** (`freeboard.enforce: true`, off by default):
   require `rho_ice*hi + rho_snow*hs <= rho_ocean*(hi + hs)` per category.
   First redistribute snow off flooded cats onto cats with headroom; if any
   cat remains negative-freeboard, grow its ice volume to lift the snow-ice
   interface back to sea level.

6. **Min-vice cleanup**: sweep cats with `vicen < min cat ice volume` and
   redistribute their `(aicen, vicen, vsnon)` mass into surviving cats
   proportionally to `aicen`. Without this the rebin's marginal-aice
   solutions get clipped at one cat's upper bin edge, producing
   thick-ice-on-edge artifacts.

7. **Aggregate diagnostics**: compute the cell totals
   `sea_ice_area_fraction`, `sea_ice_thickness`, `sea_ice_snow_thickness`
   that go into the returned State.

## Thermo / pond pass

After the per-cell pass, the per-layer thermo and melt-pond fields are
updated in place:

- On cats that lost ice (LAND, ICE2NOICE, rebin-emptied, min-vice cleanup):
  zero `Tsfcn`, all `qice`, all `sice`, all `qsno`, all pond fields.
  Avoids leaving stale background values on empty slots.

- On cats with ice and `thermo.update snow thermo: true` (default):
  clamp `Tsfcn` to `[min Tsfc, max Tsfc]`; rebuild
  `qsno = snowEnthalpy(Tsfcn)` (so qsno is consistent with the prevailing
  surface temperature, not whatever stale background value happened to be
  there); cap the surface ice layer enthalpy by
  `iceEnthalpyBL99(Tsfcn, sice)`.

- On cats where the snow distribution actually modified `vsnon` and
  `thermo.reset ponds: true` (default): zero `apnd` and `hpnd` (the new
  snow load is inconsistent with the bg pond state). `ipnd` (refrozen-lid
  thickness) is preserved.

## New-ice seeding (`thermo.seed new ice: true`, default)

For every (jnode, cat) slot that transitioned from `bg_aicen=0` to
`new_aicen>0`, SOCA needs per-layer thermo on a slot CICE never wrote.
PostProcessIce:

1. Picks a donor cell from the global lat/lon KDTree of cells that carry
   any ice within `thermo.seed search neighbors` (default 64) nearest
   neighbours. The donor sits in the same `donorCache` built once at the
   top of `runPostprocess` via a sparse halo exchange, so this step is
   purely local.
2. Seeds `Tsfcn` from the donor's area-weighted mean Tsfc, clamped to
   `[min Tsfc, max Tsfc]`, and `qsno = snowEnthalpy(Tsfcn)`.
3. Synthesizes sub-surface `sice` from the CICE5/CICE6 BZ99 layer profile
   and `qice` from `iceEnthalpyBL99(Tfrz - Tice_seed_offset, sice)`. The
   sub-surface `qice` is *not* seeded from the surface temperature: that
   would give qice values several MJ/m^3 too warm.
4. Falls back to Tfrz when no donor with ice is found within K (logged as a
   warning).

## Configuration

The full schema (showing defaults):

```yaml
postprocess ice:

  # Required: CICE category and layer counts.
  ncat: 5
  ice_lev: 7
  sno_lev: 1

  # Required: I/O paths.
  cice restart:
    input:  /path/to/cice_background.res.nc   # also the update-mode template
    output: /path/to/cice_postprocessed.res.nc

  # Cells with analysis aice below this threshold are treated as no-ice in
  # the output (all per-cat fields zeroed). Set to 0 to disable.
  min aice output: 1.0e-4

  # Per-cat ice volume floor (m^3/m^2-cell). Mass-conserving sweep after
  # the rebin. Distinct from `min new ice thickness`.
  min cat ice volume: 1.0e-5

  # Cell-mean new-ice thickness used for the rebin's volume target when an
  # noice->ice cell has no analysis hice and no donor with ice is found.
  min new ice thickness: 0.1

  # Which ice analysis fields are actually produced by the DA. Required:
  # sea_ice_area_fraction. For ice volume, list one of
  # `sea_ice_thickness` (per-ice-area, m) or `sea_ice_volume`
  # (per-cell-area, m); for snow, similarly. Missing pairs fall back to bg.
  analysis variables:
  - sea_ice_area_fraction
  - sea_ice_thickness
  - sea_ice_snow_thickness

  itd:
    rebin: true                                # default
    category bounds: [0.0, 0.6445072, 1.391433, 2.470179, 4.567288, 9.333887]
    min thickness gap: 0.01                    # numerical robustness, m

  snow:
    min snow thickness: 0.0                    # per-ice-area, m
    max snow thickness: 5.0                    # per-ice-area, m

  freeboard:
    enforce: false                             # off by default

  thermo:
    update snow thermo: true
    reset ponds: true
    seed new ice: true
    seed search neighbors: 64
    max Tsfc: -1.0
    min Tsfc: -100.0
```

`ncat`, `ice_lev`, `sno_lev`, and the `cice restart: input/output` block
are required. Everything else has a default that matches the production
GDAS configuration; you only need to set what differs.

The freeboard pass uses fixed densities from `icephysics::Constants`
(`rho_ice = 917`, `rho_snow = 330`, `rho_ocean = 1025`, all in kg/m^3) -
the CICE5/CICE6 defaults. Making them configurable would let users write
a restart inconsistent with CICE's own physics.

## Standalone driver

The standalone application `soca_postproc.x` (see [`Postproc.h`](../../mains/Postproc.h))
reads `background` + `increment`, forms `analysis = bg + incr`, and calls
`postprocess`. The standalone YAML adds top-level `geometry`, `background`,
`increment` blocks alongside `postprocess ice`. Alternatively the user can
provide an explicit `analysis` block to skip the bg+incr step.

`PostProcessIce` is also invoked per ensemble member from
[`AnalysisPostproc.h`](../../mains/AnalysisPostproc.h) (`soca_ensanpproc.x`)
and from `gdasapp`'s increment handler.

## Implementation notes

- Per-cat reads / writes go through the `ppiCatViews` / `ppiCatLevViews`
  helpers in `PostProcessIce.cc`, which gather atlas array views for the
  expanded category fields. Field-name conventions
  (`sea_ice_category<N>_*`, `sea_ice_category<N>_layer<L>_*`) come from
  `fields_metadata.yml` and are auto-expanded by the `<CATEGORY>` and
  `<LEVEL>` placeholders.
- The global KDTree is built once at construction by all-gathering owned
  `(lon, lat, gidx)` triples across all ranks. The sparse halo exchange
  (`gatherDonorHalo`) then routes per-cat donor records through two
  `comm.allToAll` rounds, so memory scales with the halo size, not the
  global grid.
- The pure column-physics helpers live in [`IcePhysics.h`](IcePhysics.h):
  `adjustThicknessCategories` (ITD rebin), `enforceFreeboard`,
  `iceEnthalpyBL99`, `siceLayerCice4` (BZ99 salinity profile),
  `snowEnthalpy`. They have no atlas / MPI / yaml dependencies and are
  exercised by `TestIcePhysics.cc`.
