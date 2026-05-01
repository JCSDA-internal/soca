# soca2cice comparison study

Side-by-side runs of the legacy Fortran `Soca2Cice` variable change vs the new
C++ `PostProcessIce` pipeline on real GDAS cycles.

This is a **study tool, not a CI test**.

## Layout

```
postproc_compare/
├── templates/
│   ├── convertstate.yml.tmpl   # legacy path (soca_convertstate.x)
│   └── postproc.yml.tmpl       # new path    (soca_postproc.x)
├── make_yamls.py               # render both yamls + symlink inputs for one cycle
└── run_all.sh                  # loop over cycles, run both binaries
```

Generated outputs live **outside the repo** under whatever `<out_root>` you pass
to `run_all.sh`.

## Inputs expected per cycle (GDAS layout)

For `cycle_tag = gdas.YYYYMMDD`:

| role             | path                                                                     |
|------------------|--------------------------------------------------------------------------|
| bkg CICE restart | `gdas.YYYYMMDD/00/model/ice/restart/YYYYMMDD.030000.cice_model.res.nc`   |
| bkg MOM restart  | `gdas.YYYYMMDD/00/model/ocean/restart/YYYYMMDD.030000.MOM.res.nc`        |
| ice analysis     | `gdas.YYYYMMDD/06/analysis/ice/gdas.t06z.jedi_analysis.a006.nc`          |
| ocn analysis     | `gdas.YYYYMMDD/06/analysis/ocean/gdas.t06z.jedi_analysis.a006.nc`        |

Background date in the yaml = `YYYY-MM-DDT03:00:00Z`, analysis date =
`YYYY-MM-DDT06:00:00Z`.

If your tree doesn't match (e.g. MOM restart path differs), edit the source
paths near the top of `make_yamls.py`.

## Geometry (provide once, shared across cycles)

`run_all.sh` expects three files in `<out_root>/_shared/`:

- `soca_gridspec.nc`
- `input.nml`
- `fields_metadata.yml`

Symlinks are fine. They're symlinked into each per-cycle workdir.

## Run

```bash
cd soca/tools/postproc_compare
./run_all.sh <gdas_root> <out_root> <mpi_ranks> [cycle_tag ...]
```

With no `cycle_tag` args, it runs the default 7-cycle list:

```
gdas.20240715  gdas.20250115  gdas.20250330
gdas.20250715  gdas.20250920  gdas.20251221
gdas.20260322
```

Example:

```bash
./run_all.sh /scratch/gdas /scratch/soca2cice_compare 8
```

## Per-cycle workdir

Each cycle gets `<out_root>/YYYYMMDD/` containing:

```
soca_gridspec.nc       -> ../_shared/soca_gridspec.nc
input.nml              -> ../_shared/input.nml
fields_metadata.yml    -> ../_shared/fields_metadata.yml
background.cice.res.nc -> <gdas tree>            # bkg CICE restart
background.MOM.res.nc  -> <gdas tree>            # bkg MOM restart
analysis.ice.nc        -> <gdas tree>            # SOCA-grid analysis
analysis.ocn.nc        -> <gdas tree>            # SOCA-grid analysis
convertstate.yml                                 # generated
postproc.yml                                     # generated
iced.legacy.nc                                   # legacy Fortran output
iced.new.nc                                      # new C++ output
convertstate.log
postproc.log
```

## Tunables

In `make_yamls.py`:

- `NCAT`, `ICE_LEV`, `SNO_LEV` — CICE category / vertical layer counts (default 5/7/1)
- `BKG_RESTART_HHMMSS`, `BKG_DATE`, `AN_VALID` — cycle timestamps
- the `bkg_cice`, `bkg_mom`, `ice_an`, `ocn_an` path expressions — the GDAS layout

## Run a single cycle by hand

```bash
python3 make_yamls.py /scratch/gdas gdas.20240715 /scratch/soca2cice_compare/20240715
cd /scratch/soca2cice_compare/20240715
mpiexec -n 8 $BUNDLE/build/bin/soca_convertstate.x convertstate.yml
mpiexec -n 8 $BUNDLE/build/bin/soca_postproc.x     postproc.yml
```
