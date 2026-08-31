#!/usr/bin/env bash
# E03 sweep: every admissible e (2<=e<=n-2, gcd(e,n)=1) for the given n, with a per-run time budget.
# A run that exceeds the budget is re-run with --no-mac4 --censored, so its D=4 columns stay empty
# (a timed-out run is a lower bound, never a value -- CLAUDE.md ground rule 4).
# Usage: bash sweep.sh [n ...]        (default 11 12 13 14 15 16)
set -uo pipefail
cd "$(dirname "$0")/../.."
NS=${@:-"11 12 13 14 15 16"}
for n in $NS; do
  # budget: 20 min up to n=13, 40 min for n>=14
  if [ "$n" -ge 14 ]; then BUDGET=2400; else BUDGET=1200; fi
  for e in $(seq 2 $((n-2))); do
    python3 -c "from math import gcd; import sys; sys.exit(0 if gcd($e,$n)==1 else 1)" || continue
    echo "=== n=$n e=$e (budget ${BUDGET}s) ==="
    if ! timeout "$BUDGET" python3 experiments/E03_syzygy_and_degree_fall/sbox_degree_fall.py "$n" "$e"; then
      echo "!!! n=$n e=$e exceeded ${BUDGET}s at D=4 -> re-running with --no-mac4 --censored"
      timeout "$BUDGET" python3 experiments/E03_syzygy_and_degree_fall/sbox_degree_fall.py "$n" "$e" --no-mac4 --censored \
        || echo "!!! n=$n e=$e censored at D=3 as well (no row written)"
    fi
  done
done
echo "sweep done"
