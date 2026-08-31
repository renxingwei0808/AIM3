#!/usr/bin/env bash
# E04 Tier-A batch: n = 7, 9, 11 for l=2 (three S-boxes) and l=1 (two S-boxes), 2 seeds each.
# --per-block adds the direct-sum control (each S-box's 3n quadratics alone in the same variables).
# Every full-system D=4 run at these sizes turns out to be SATURATED (see README); the runs are kept
# because that saturation, and the fact that the random control saturates identically, is the result.
set -uo pipefail
cd "$(dirname "$0")/../.."
E=experiments/E04_full_system_hilbert/full_system_falls.py
for l in 2 1; do
  for n in 7 9; do
    python3 $E $n -l $l --seeds 1,2 --systems S_quad,P1only,random --budget 1800 --per-block
  done
done
# n=11: the structured systems first (they are sparse and cheap), the dense random control last
for l in 2 1; do
  python3 $E 11 -l $l --seeds 1,2 --systems S_quad,P1only --budget 3600 --per-block
done
for l in 2 1; do
  python3 $E 11 -l $l --seeds 1 --systems random --budget 1800
done
echo "E04 tier-A batch done"
