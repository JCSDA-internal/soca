# SOCA ensemble postprocessing (`soca_anpproc.x`)

The behavior below follows the execution order in `soca/src/mains/AnalysisPostproc.h`.

## Sample base YAML for core inputs

```yaml
geometry: ...
nens: <N>
nens per MPI task: <n_local>
increment variables: [<vars used in increments>]

backgrounds:         # ensemble backgrounds in
  members from template: ...

output increments:   # processed increments out
  datadir: <out_dir>/
```

## Regime 1: background-perturbation regime

This is active when `analysis increments` is **absent**.

What the code does:
- Reads `backgrounds`, computes ensemble mean and perturbations.
- If `ensemble inflation` is present:
  - `value`: multiplies perturbations by scalar.
  - `field`: applies Schur product with a weight field.
  - Converts inflated perturbations into increments by subtracting original perturbations.
- If `ensemble inflation` is absent: increments are set to zero.

Sample YAML for ensemble inflation:

```yaml
# no "analysis increments" block
ensemble inflation:
  value: 1.05
```

## Regime 2: analysis-increment regime

This is active when `analysis increments` is present.

What the code does:
- Reads member-wise `analysis increments`.
- If `ensemble inflation` exists, applies inflation to those analysis increments.
- Computes analysis ensemble mean as: `background mean + mean(analysis increments)`.
- Sets working increments to the analysis increments (these flow to later regimes).

Sample YAML for analysis increments:

```yaml
analysis increments:
  members from template:
    template:
      basename: <analysis_incr_dir>/
      ocn_filename: <incr>.%mem%.nc
```

### Inflation in Regime 2

When `analysis increments` and `ensemble inflation` are present, inflation is
applied directly to those increments before recentering/postprocessing,
using OOPS inflation methods.

Supported `method` values in this path are:
- `RTPS` (relaxation to prior spread): rescales analysis perturbation amplitude so
  the analysis spread is relaxed toward background spread.
- `RTPP` (relaxation to prior perturbations): blends analysis perturbations with
  background perturbations using `factor`.
- `Multiplicative`: scales analysis increments by configured factor(s)
  (globally with `factor`, or by level using `levels` + `factors`).

Notes:
- Regime 1 style inflation settings and Regime 2 style inflation settings are
  different and not interchangeable.
- If `ensemble inflation` is absent, analysis increments are used unchanged.

Sample YAML configurations for different kinds of inflation:

```yaml
ensemble inflation:
  method: RTPS
  factor: 0.8
```

```yaml
ensemble inflation:
  method: RTPP
  factor: 0.5
```

```yaml
ensemble inflation:
  method: Multiplicative
  factor: 1.05
```

## Regime 3: recentering regime

This is active when `recentering state` is present.

What the code does:
- Reads deterministic `recentering state`.
- Builds recentering increment as:
  - `recentering state - current ensemble mean` (current ensemble mean is
    either background mean for Regime 1, or analysis mean for Regime 2).
- Adds this recentering increment to each member increment.

A sample YAML for the deterministic recentering state:

```yaml
recentering state:
  basename: <center_dir>/
  ocn_filename: MOM.res.nc
  ice_filename: cice.res.nc
```

## Regime 4: increment postprocessing regime

This is active when `increment postprocessing` is present.

What the code does (in this order):
1. `append vertical geometry`: appends a layer field (for MOM6 IAU use).
2. `set increment variables to zero`: zeroes selected variables (must be subset of `increment variables`).
3. `change precision`: rounds selected variables to configured precision.
4. `bounds check`: runs `soca::incrqc::qcIncrement` member-by-member.

For details of `bounds check` options and behavior, see
[`soca::incrqc` documentation](../src/soca/Utils/incrqc/README.md).

A sample yaml:

```yaml
increment postprocessing:
  append vertical geometry:
    layers variable: sea_water_cell_thickness
    vertical geometry:
      basename: <geom_dir>/
      ocn_filename: MOM.res.nc
  set increment variables to zero:
  - eastward_sea_water_velocity
  - northward_sea_water_velocity
  change precision:
    variables: [sea_water_salinity]
    precision: 1.0e-5
  # Keep the analysis within physical bounds
  bounds check:
    state bounds:
      sea_water_potential_temperature: [-2.5, 36.0]
      sea_water_salinity: [0.0, 44.0]
    absolute steric increment max: 0.5
    steric variable change:
      linear variable changes:
      - linear variable change name: BalanceSOCA
    coastal increment filter:
      min distance: 0.0       # [m] zero increments at the coast
      max distance: 200000.0  # [m] full increments beyond 200 km from coast
```

## Regime 5: analysis sea-ice postprocessing regime

This is active when `analysis postprocessing` is present.

What the code does:
- For each ensemble member:
  - Adds final increment to the background state.
  - Applies `Soca2Cice` variable change and writes CICE restart outputs, the
    configuration comes from the `analysis postprocessing.sea ice variable change`.

Sample YAML:

```yaml
analysis postprocessing:
  sea ice variable change:
    pattern: "%mem%"
    variable change name: Soca2Cice
    cice background state:
      restart: <cice_bg>
      ncat: 5
      ice_lev: 7
      sno_lev: 1
    cice output:
      restart: <out_dir>/mem.%mem%/iced.nc
```

## Three common ensemble use cases

1. **Recenter + postprocess ensemble backgrounds**
   - Use Regimes 1 + 3 + 4 + 5.
   - Omit `analysis increments`, include `recentering state`.

2. **Recenter + postprocess ensemble analyses**
   - Use Regimes 2 + 3 + 4 + 5.
   - Include `analysis increments` and `recentering state`.

3. **Postprocess ensemble analyses only**
   - Use Regimes 2 + 3 + 5.
   - Omit `recentering state`.

## Run

```bash
mpiexec -n <MPI_TASKS> build/bin/soca_anpproc.x <config>.yml
```
