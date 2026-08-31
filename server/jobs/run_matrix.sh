#!/usr/bin/env bash
# Farm one experiment over a parameter grid with GNU parallel. Example:
#   ./run_matrix.sh E03 "python3 experiments/E03_syzygy_and_degree_fall/sbox_degree_fall.py {1} {2}" "11 12 13 14" "3 5 7"
set -euo pipefail
EXP=$1; CMD=$2; NS=$3; ES=$4
OUT=results/server/$(date +%Y%m%d)_$EXP; mkdir -p "$OUT"
parallel --jobs "${JOBS:-40}" --joblog "$OUT/joblog.txt" --results "$OUT/{1}_{2}" \
  "ulimit -v ${MEM_KB:-33554432}; $CMD" ::: $NS ::: $ES
