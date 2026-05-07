"""Paths and constants for the Fortran-vs-C++ PostProcessIce comparison.

These scripts are driven by `compare_cycles.py`, which sets all input/output
paths via env vars. Direct invocation requires the env vars below:

- POSTPROC_COMPARE_FORTRAN     : legacy Fortran output (.nc)
- POSTPROC_COMPARE_CPP         : new C++ output (.nc)
- POSTPROC_COMPARE_GRIDSPEC    : gridspec (.nc)
- POSTPROC_COMPARE_PLOTS_DIR   : output base for plots/
- POSTPROC_COMPARE_REPORTS_DIR : output base for diff_stats/ (defaults to PLOTS_DIR)
"""
import os
from pathlib import Path


def _required(name):
    val = os.environ.get(name)
    if not val:
        raise RuntimeError(
            f"{name} is not set. Run via compare_cycles.py, or set "
            f"POSTPROC_COMPARE_{{FORTRAN,CPP,GRIDSPEC,PLOTS_DIR}} explicitly."
        )
    return Path(val)


FORTRAN_OUT = _required("POSTPROC_COMPARE_FORTRAN")
CPP_OUT     = _required("POSTPROC_COMPARE_CPP")
GRIDSPEC    = _required("POSTPROC_COMPARE_GRIDSPEC")
PLOTS_DIR   = _required("POSTPROC_COMPARE_PLOTS_DIR")
REPORTS_DIR = Path(os.environ.get("POSTPROC_COMPARE_REPORTS_DIR", PLOTS_DIR))

RAW_DIR       = PLOTS_DIR / "raw"
AGG_DIR       = PLOTS_DIR / "aggregate"
DIFFSTATS_DIR = REPORTS_DIR / "diff_stats"

for d in (RAW_DIR, AGG_DIR, DIFFSTATS_DIR):
    d.mkdir(parents=True, exist_ok=True)

PER_CAT_VARS = ["aicen", "vicen", "vsnon", "Tsfcn", "sice001", "qice001", "qsno001"]
DEFAULT_CAT  = 0
REGIONS      = ["Arctic", "Antarctic"]


def get_cat_dim(ds):
    for d in ("cat", "ncat", "zaxis_1"):
        if d in ds.dims:
            return d
    return None


CANON_DIMS = {
    "ncat": "cat", "zaxis_1": "cat",
    "nj": "y", "yaxis_1": "y",
    "ni": "x", "xaxis_1": "x",
}


def open_squeezed(path, variables=None):
    """Open lazily; if `variables` given, drop all other data variables to bound memory.

    Renames dims to a canonical (cat, y, x) so arrays from the two files are
    broadcast-compatible — without this, xarray will broadcast across all 6
    dims and explode memory.
    """
    import netCDF4
    import xarray as xr
    drop = None
    if variables is not None:
        keep = set(variables)
        with netCDF4.Dataset(path) as nc:
            drop = [v for v in nc.variables.keys() if v not in keep and v not in nc.dimensions]
    ds = xr.open_dataset(path, drop_variables=drop)
    if "Time" in ds.dims:
        ds = ds.isel(Time=0)
    rename = {old: new for old, new in CANON_DIMS.items() if old in ds.dims}
    if rename:
        ds = ds.rename(rename)
    return ds
