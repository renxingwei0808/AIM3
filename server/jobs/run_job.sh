#!/usr/bin/env bash
# One accounted server job (batch rule 0.2): tmux + nohup + ulimit + /usr/bin/time -v.
#
#   server/jobs/run_job.sh EXP NAME MEM_GB -- command args...
#
# * stdout/stderr  -> results/server/<YYYYMMDD>_EXP/NAME.log
# * /usr/bin/time -v (wall clock, peak RSS, exit status) -> .../NAME.time
# * a one-line provenance stamp (host, commit, ulimit, command) -> .../NAME.meta
# * MEM_GB sets `ulimit -v`; a job killed by it exits non-zero and the caller must record
#   censored=1 (its numbers are lower bounds, never values -- CLAUDE.md ground rule 4).
#
# Runs in the foreground; the caller is expected to background it (`tmux new-session -d`).
set -uo pipefail

EXP=$1; NAME=$2; MEM_GB=$3; shift 3
[ "${1:-}" = "--" ] && shift

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
# AIM3_DATE must match batch1.sh, or the launcher writes its JOB_DONE stamp into a
# different directory from the one holding the log
OUT="$ROOT/results/server/${AIM3_DATE:-$(date +%Y%m%d)}_$EXP"
mkdir -p "$OUT"

LOG="$OUT/$NAME.log"
TIMEF="$OUT/$NAME.time"
META="$OUT/$NAME.meta"

{
  echo "host    : $(hostname)"
  echo "commit  : $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "branch  : $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
  echo "started : $(date -Is)"
  echo "ulimit -v: ${MEM_GB} GB"
  echo "results : ${AIM3_RESULTS:-<repo results/>}"
  echo "scratch : ${AIM3_SCRATCH:-<tmp>}"
  echo "command : $*"
} > "$META"

ulimit -v $((MEM_GB * 1024 * 1024)) 2>/dev/null || echo "WARNING: could not set ulimit -v" >> "$META"

/usr/bin/time -v -o "$TIMEF" "$@" > "$LOG" 2>&1
rc=$?

{
  echo "finished: $(date -Is)"
  echo "exit    : $rc"
  grep -E "Elapsed \(wall clock\)|Maximum resident set size" "$TIMEF" 2>/dev/null | sed 's/^\s*/time    : /'
} >> "$META"
exit $rc
