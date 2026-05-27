#!/usr/bin/env python3
"""Generate convertstate.yml (legacy) and postproc.yml (new) for one cycle.

Usage: make_yamls.py <gdas_root> <cycle_tag> <output_dir>

  gdas_root   parent dir containing gdas.YYYYMMDD/...
  cycle_tag   e.g. gdas.20240715 (background at 00Z, analysis at 06Z)
  output_dir  where to write the two yamls (also the workdir for runs)
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TMPL = HERE / "templates"

# CICE category / vertical layer counts (GDAS production = 5/7/1).
NCAT = 5
ICE_LEV = 7
SNO_LEV = 1

# Time offsets relative to YYYYMMDD: bkg cycle 00Z with restart written at 03Z,
# analysis cycle 06Z with valid time 06Z.
BKG_CYCLE = "00"
BKG_RESTART_HHMMSS = "030000"
BKG_DATE = "T03:00:00Z"
AN_CYCLE = "06"
AN_VALID = "T06:00:00Z"


def render(tmpl_path: Path, subs: dict) -> str:
    txt = tmpl_path.read_text()
    for k, v in subs.items():
        txt = txt.replace("{{" + k + "}}", str(v))
    return txt


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    gdas_root = Path(sys.argv[1]).resolve()
    cycle_tag = sys.argv[2]                   # gdas.20240715
    out_dir = Path(sys.argv[3]).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ymd = cycle_tag.split(".")[1]             # 20240715
    bkg_dt = f"{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}{BKG_DATE}"
    an_dt = f"{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}{AN_VALID}"

    # Source paths inside gdas tree
    bkg_cice = (gdas_root / cycle_tag / BKG_CYCLE / "model" / "ice" / "restart"
                / f"{ymd}.{BKG_RESTART_HHMMSS}.cice_model.res.nc")
    bkg_mom = (gdas_root / cycle_tag / BKG_CYCLE / "model" / "ocean" / "restart"
               / f"{ymd}.{BKG_RESTART_HHMMSS}.MOM.res.nc")
    ice_an = (gdas_root / cycle_tag / AN_CYCLE / "analysis" / "ice"
              / f"gdas.t{AN_CYCLE}z.jedi_analysis.a006.nc")
    ocn_an = (gdas_root / cycle_tag / AN_CYCLE / "analysis" / "ocean"
              / f"gdas.t{AN_CYCLE}z.jedi_analysis.a006.nc")

    # Symlink (or note) inputs into out_dir under simple basenames so the
    # yamls can be relative.  Geometry files (soca_gridspec.nc, input.nml,
    # fields_metadata.yml) are user-supplied and expected in out_dir already.
    bkg_cice_lnk = "background.cice.res.nc"
    bkg_mom_lnk = "background.MOM.res.nc"
    ice_an_lnk = "analysis.ice.nc"
    ocn_an_lnk = "analysis.ocn.nc"
    for src, name in [(bkg_cice, bkg_cice_lnk),
                      (bkg_mom, bkg_mom_lnk),
                      (ice_an, ice_an_lnk),
                      (ocn_an, ocn_an_lnk)]:
        dst = out_dir / name
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        if src.exists():
            dst.symlink_to(src)
        else:
            print(f"WARNING: missing input {src}", file=sys.stderr)
            dst.symlink_to(src)  # leave dangling so the failure is obvious

    legacy_subs = {
        "BKG_CICE_RESTART": bkg_cice_lnk,
        "NCAT": NCAT, "ICE_LEV": ICE_LEV, "SNO_LEV": SNO_LEV,
        "OUT_LEGACY": "iced.legacy.nc",
        "OCN_ANALYSIS": ocn_an_lnk,
        "ICE_ANALYSIS": ice_an_lnk,
        "ANALYSIS_DATE": an_dt,
    }
    new_subs = {
        "BKG_CICE_RESTART": bkg_cice_lnk,
        "OCN_ANALYSIS": ocn_an_lnk,
        "ICE_ANALYSIS": ice_an_lnk,
        "ANALYSIS_DATE": an_dt,
        "NCAT": NCAT, "ICE_LEV": ICE_LEV, "SNO_LEV": SNO_LEV,
        "OUT_NEW": "iced.new.nc",
    }

    (out_dir / "convertstate.yml").write_text(
        render(TMPL / "convertstate.yml.tmpl", legacy_subs))
    (out_dir / "postproc.yml").write_text(
        render(TMPL / "postproc.yml.tmpl", new_subs))

    print(f"wrote {out_dir}/convertstate.yml")
    print(f"wrote {out_dir}/postproc.yml")


if __name__ == "__main__":
    main()
