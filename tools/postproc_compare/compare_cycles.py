#!/usr/bin/env python3
"""Per-cycle driver for the Fortran-vs-C++ comparison plots/reports.

Iterates over cycle workdirs produced by `run_all.sh` (each contains
`iced.legacy.nc` and `iced.new.nc`), and runs the three scripts in
`comparison_fortran_vs_cpp/scripts/` for each one. Outputs land under
`soca/tools/postproc_compare/{plots,reports}/<cycle>/`.

Usage:
    python3 compare_cycles.py <out_root> [--gridspec PATH] [--cycles tag1 tag2 ...]

`<out_root>` is the same directory passed to `run_all.sh`. Without `--cycles`,
the driver picks up every subdirectory of `<out_root>` containing both
`iced.legacy.nc` and `iced.new.nc`. The default gridspec is
`<out_root>/_shared/soca_gridspec.nc`.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = THIS_DIR / "scripts"
PLOTS_BASE = THIS_DIR / "plots"
REPORTS_BASE = THIS_DIR / "reports"
DEFAULT_PYTHON = "python3"

CYCLE_SCRIPTS = ["plot_panels.py", "plot_aggregates.py", "diff_stats.py"]


def discover_cycles(out_root):
    cycles = []
    for entry in sorted(out_root.iterdir()):
        if not entry.is_dir() or entry.name.startswith("_"):
            continue
        if (entry / "iced.legacy.nc").exists() and (entry / "iced.new.nc").exists():
            cycles.append(entry)
    return cycles


def run_cycle(cycle_dir, gridspec, python_bin):
    cycle_tag = cycle_dir.name
    plots_dir = PLOTS_BASE / cycle_tag
    reports_dir = REPORTS_BASE / cycle_tag
    plots_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["POSTPROC_COMPARE_FORTRAN"]     = str(cycle_dir / "iced.legacy.nc")
    env["POSTPROC_COMPARE_CPP"]         = str(cycle_dir / "iced.new.nc")
    env["POSTPROC_COMPARE_GRIDSPEC"]    = str(gridspec)
    env["POSTPROC_COMPARE_PLOTS_DIR"]   = str(plots_dir)
    env["POSTPROC_COMPARE_REPORTS_DIR"] = str(reports_dir)

    print(f"\n=== {cycle_tag} ===")
    print(f"  legacy:   {env['POSTPROC_COMPARE_FORTRAN']}")
    print(f"  new:      {env['POSTPROC_COMPARE_CPP']}")
    print(f"  plots:    {plots_dir}")
    print(f"  reports:  {reports_dir}")

    failed = []
    for script in CYCLE_SCRIPTS:
        cmd = [str(python_bin), str(SCRIPTS_DIR / script)]
        rc = subprocess.run(cmd, env=env, cwd=str(SCRIPTS_DIR)).returncode
        if rc != 0:
            failed.append(script)
            print(f"  ! {script} failed (exit {rc})")
    return failed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_root", type=Path,
                    help="root passed to run_all.sh; expects per-cycle subdirs.")
    ap.add_argument("--gridspec", type=Path, default=None,
                    help="gridspec netCDF (default: <out_root>/_shared/soca_gridspec.nc)")
    ap.add_argument("--cycles", nargs="*", default=None,
                    help="restrict to these cycle tags (subdir names)")
    ap.add_argument("--python", default=DEFAULT_PYTHON,
                    help=f"python interpreter (default: {DEFAULT_PYTHON})")
    args = ap.parse_args()

    out_root = args.out_root.resolve()
    if not out_root.is_dir():
        ap.error(f"{out_root} is not a directory")

    gridspec = args.gridspec or (out_root / "_shared" / "soca_gridspec.nc")
    if not gridspec.exists():
        ap.error(f"gridspec not found: {gridspec}")

    if args.cycles:
        cycles = []
        for tag in args.cycles:
            d = out_root / tag
            if not (d / "iced.legacy.nc").exists() or not (d / "iced.new.nc").exists():
                print(f"skip {tag}: missing iced.legacy.nc or iced.new.nc", file=sys.stderr)
                continue
            cycles.append(d)
    else:
        cycles = discover_cycles(out_root)

    if not cycles:
        ap.error(f"no cycle dirs with iced.legacy.nc + iced.new.nc under {out_root}")

    print(f"Driver: {len(cycles)} cycle(s)")
    print(f"Gridspec: {gridspec}")
    print(f"Plots base:   {PLOTS_BASE}")
    print(f"Reports base: {REPORTS_BASE}")

    any_failed = False
    for cycle_dir in cycles:
        fails = run_cycle(cycle_dir, gridspec, args.python)
        if fails:
            any_failed = True
            print(f"  cycle {cycle_dir.name}: {len(fails)} script(s) failed: {fails}")

    sys.exit(1 if any_failed else 0)


if __name__ == "__main__":
    main()
