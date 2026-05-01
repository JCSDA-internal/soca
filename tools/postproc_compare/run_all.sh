#!/usr/bin/env bash
# Loop over GDAS cycles, generate yamls, and run both legacy + new soca2cice.
#
# Usage:
#   run_all.sh <gdas_root> <out_root> <mpi_ranks> [cycle_tag ...]
#
# Example:
#   run_all.sh /scratch/gdas /scratch/soca2cice_compare 8
#
# Expects soca_gridspec.nc, input.nml, fields_metadata.yml to live in
# <out_root>/_shared/ — they are symlinked into each per-cycle dir.

set -euo pipefail

if [[ $# -lt 3 ]]; then
  sed -n '2,15p' "$0"
  exit 1
fi

GDAS_ROOT=$1; shift
OUT_ROOT=$1; shift
MPI=$1; shift

DEFAULT_CYCLES=(
  gdas.20240715 gdas.20250115 gdas.20250330
  gdas.20250715 gdas.20250920 gdas.20251221
  gdas.20260322
)
if [[ $# -gt 0 ]]; then
  CYCLES=("$@")
else
  CYCLES=("${DEFAULT_CYCLES[@]}")
fi

HERE="$(cd "$(dirname "$0")" && pwd)"
BUNDLE="$(cd "$HERE/../../.." && pwd)"
BIN="$BUNDLE/build/bin"

SHARED="$OUT_ROOT/_shared"
if [[ ! -f "$SHARED/soca_gridspec.nc" || ! -f "$SHARED/input.nml" || ! -f "$SHARED/fields_metadata.yml" ]]; then
  echo "ERROR: place soca_gridspec.nc, input.nml, fields_metadata.yml in $SHARED/" >&2
  exit 2
fi

mkdir -p "$OUT_ROOT"

for tag in "${CYCLES[@]}"; do
  ymd=${tag#gdas.}
  case_dir="$OUT_ROOT/${ymd}"
  echo "=== $tag → $case_dir ==="
  mkdir -p "$case_dir"
  for f in soca_gridspec.nc input.nml fields_metadata.yml; do
    ln -sfn "$SHARED/$f" "$case_dir/$f"
  done

  python3 "$HERE/make_yamls.py" "$GDAS_ROOT" "$tag" "$case_dir"

  pushd "$case_dir" >/dev/null

  echo "--- legacy (soca_convertstate.x) ---"
  mpiexec -n "$MPI" "$BIN/soca_convertstate.x" convertstate.yml \
      > convertstate.log 2>&1 || echo "  legacy FAILED (see convertstate.log)"

  echo "--- new (soca_postproc.x) ---"
  mpiexec -n "$MPI" "$BIN/soca_postproc.x" postproc.yml \
      > postproc.log 2>&1 || echo "  new FAILED (see postproc.log)"

  popd >/dev/null
done

echo "done. outputs under $OUT_ROOT/<ymd>/{iced.legacy.nc, iced.new.nc}"
