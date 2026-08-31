#!/usr/bin/env bash
# Run E03's own summariser over the Tier-A table PLUS this batch's n = 17..20 rows (task 4).
#
#   bash server/jobs/e03_summary_merged.sh [OUTDIR]
#
# `summarise.py --md` writes `E03_summary.md` into $AIM3_RESULTS.  The committed
# `results/E03_summary.md` is a Tier-A artefact and is deliberately NOT overwritten by this batch
# (the batch may only add files), so the merged table is written to
# `results/server/<date>_E03/E03_summary_with_n17_20.md` instead, and the exit status of
# summarise.py (non-zero iff the calibrated rule has a mismatch) is reported.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
OUT=${1:-results/server/$(date +%Y%m%d)_E03}
mkdir -p "$OUT"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp results/E03_sbox_degree_fall.csv "$TMP/E03_sbox_degree_fall.csv"
if [ -f results/E03_server_n17_20.csv ]; then
  tail -n +2 results/E03_server_n17_20.csv >> "$TMP/E03_sbox_degree_fall.csv"
else
  echo "WARNING: results/E03_server_n17_20.csv is missing -- summarising Tier A only" >&2
fi
rows=$(( $(wc -l < "$TMP/E03_sbox_degree_fall.csv") - 1 ))

AIM3_RESULTS="$TMP" python3 experiments/E03_syzygy_and_degree_fall/summarise.py --md > /dev/null
rc=$?
cp "$TMP/E03_summary.md" "$OUT/E03_summary_with_n17_20.md"
echo "merged table: $rows runs -> $OUT/E03_summary_with_n17_20.md"
if [ $rc -eq 0 ]; then
  echo "summarise.py exit 0: the calibrated fall-count rule reproduces every completed run"
else
  echo "summarise.py exit $rc: the rule has at least one mismatch (listed in the file)"
fi
exit $rc
