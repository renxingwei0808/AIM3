#!/usr/bin/env bash
# Run E03's own summariser over the Tier-A table PLUS this batch's n = 17..20 rows (task 4).
#
#   bash server/jobs/e03_summary_merged.sh [OUTDIR]
#
# `summarise.py --md` writes `E03_summary.md` into $AIM3_RESULTS.  The committed
# `results/E03_summary.md` is a Tier-A artefact and is deliberately NOT overwritten by this batch
# (the batch may only add files), so the merged table goes to
# `results/server/<date>_E03/E03_summary_with_n17_20.md`, and the exit status of summarise.py
# (non-zero iff the calibrated rule has a mismatch) is reported.
#
# Rows produced with `--no-points` have empty point-set columns, which `summarise.py` reads as
# integers; they are filtered out here rather than by rewriting summarise.py, which this batch may
# not touch.  Those rows still appear in `results/REPORT_server_batch1.md` section A.3.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
OUT=${1:-results/server/$(date +%Y%m%d)_E03}
mkdir -p "$OUT"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
MERGED="$TMP/E03_sbox_degree_fall.csv"
cp results/E03_sbox_degree_fall.csv "$MERGED"

skipped=0
if [ -f results/E03_server_n17_20.csv ]; then
  skipped=$(python3 server/jobs/_e03_merge_rows.py results/E03_server_n17_20.csv "$MERGED")
else
  echo "WARNING: results/E03_server_n17_20.csv is missing -- summarising Tier A only" >&2
fi
rows=$(( $(wc -l < "$MERGED") - 1 ))

AIM3_RESULTS="$TMP" python3 experiments/E03_syzygy_and_degree_fall/summarise.py --md > /dev/null
rc=$?
cp "$TMP/E03_summary.md" "$OUT/E03_summary_with_n17_20.md"
echo "merged table: $rows runs -> $OUT/E03_summary_with_n17_20.md"
if [ "${skipped:-0}" -gt 0 ]; then
  echo "($skipped server rows excluded: run with --no-points, so summarise.py cannot read their"
  echo " point-set columns; their fall counts are checked in REPORT_server_batch1.md section A.3)"
fi
if [ $rc -eq 0 ]; then
  echo "summarise.py exit 0: the calibrated fall-count rule reproduces every merged run"
else
  echo "summarise.py exit $rc: the rule has at least one mismatch (listed in the file)"
fi
exit $rc
