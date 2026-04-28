#!/usr/bin/env bash
#
# Regression test: compare CICE restart files produced by the legacy Fortran
# Soca2Cice variable change (run via ConvertState) and the new C++
# PostProcessIce / CiceRestartIO path (run via Postproc) on identical inputs.
#
# Both paths are expected to agree to high precision on aicen, vicen, vsnon.
# Small differences in the marginal ice zone are tolerated because the two
# implementations pick neighbour donors slightly differently.
#
# Args:
#   $1  path to the Fortran-path output CICE restart
#   $2  path to the C++-path output CICE restart
#   $3  per-cell tolerance (passed to nccmp -t)
#   $4  max number of cells allowed to differ at the given tolerance
#
# Exits 0 on success, non-zero on regression.

set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 <fortran.nc> <cpp.nc> <tol> <max_diff_cells>" >&2
  exit 2
fi

FORTRAN_NC=$1
CPP_NC=$2
TOL=$3
MAX_DIFF=$4

if [[ ! -f "$FORTRAN_NC" ]]; then
  echo "ERROR: Fortran output not found: $FORTRAN_NC" >&2
  exit 1
fi
if [[ ! -f "$CPP_NC" ]]; then
  echo "ERROR: C++ output not found: $CPP_NC" >&2
  exit 1
fi

# Locate nccmp. Prefer one already on PATH; fall back to known spack-stack path.
NCCMP=$(command -v nccmp || true)
if [[ -z "$NCCMP" ]]; then
  for c in /home/annash/spack-stack/envs/unified-env.gcc/install/gcc/*/nccmp-*/bin/nccmp; do
    [[ -x "$c" ]] && NCCMP=$c && break
  done
fi
if [[ -z "$NCCMP" ]]; then
  echo "ERROR: nccmp not found on PATH" >&2
  exit 1
fi

echo "Using nccmp: $NCCMP"
echo "Comparing:"
echo "  Fortran: $FORTRAN_NC"
echo "  C++:     $CPP_NC"
echo "  Tolerance: $TOL"
echo "  Max differing cells allowed: $MAX_DIFF"

# Run nccmp with -d (data) -N (continue past first diff) -t (tolerance), and
# limit to the variables both paths actually update. Capture the per-line
# diff report and count.
TMPOUT=$(mktemp)
trap 'rm -f "$TMPOUT"' EXIT
"$NCCMP" -d -N -t "$TOL" -v aicen,vicen,vsnon "$FORTRAN_NC" "$CPP_NC" > "$TMPOUT" 2>&1 \
  || true   # nccmp returns non-zero when files differ; we evaluate manually

DIFF_LINES=$(grep -c '^DIFFER' "$TMPOUT" || true)
echo "Cells differing beyond tolerance $TOL: $DIFF_LINES"

if [[ "$DIFF_LINES" -gt "$MAX_DIFF" ]]; then
  echo "REGRESSION: too many differing cells ($DIFF_LINES > $MAX_DIFF)" >&2
  echo "First lines of diff report:" >&2
  head -20 "$TMPOUT" >&2
  exit 1
fi

echo "PASS"
