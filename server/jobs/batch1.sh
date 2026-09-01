#!/usr/bin/env bash
# Server batch 1 launcher.  Each stage starts one detached tmux session per independent job, so a
# job's CSV writes never race another's (every job gets its own $AIM3_RESULTS directory).
#
#   server/jobs/batch1.sh <stage> [args]
#
# Stages
#   e03 N...            E03 sweep, one session per n            (task 4)
#   e03d5 N...          E03 D=5 probe at n                      (task 4, second half)
#   gate1 ELL N...      E04 --per-block Gate-1 runs             (task 3)
#   e12 N...            E12 first-layer D=4                     (task 5)
#   e06 N...            E06 kernel model D=4                    (task 5)
#   e06 N...            E06 kernel model D=4                    (task 5)
#   f4export            export every system for the F4 runs      (task 6)
#   f4 TAGSUBSTR        F4 step degrees, msolve + Magma          (task 6 + E10)
#   e05 NS              E05 resultants, Magma                    (task 7)
#   cmd NAME MEM ...    run an arbitrary command as an accounted job (EXP=<exp> to file it)
#   status              list sessions and the tail of every log
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE=${AIM3_BASE:-/home/luck/rxw/AIM3_remote}
export AIM3_SAGE=${AIM3_SAGE:-$BASE/opt/envs/sage10/bin/sage}
DATE=$(date +%Y%m%d)

launch() {          # launch <session> <exp> <name> <mem_gb> <command...>
  local sess=$1 exp=$2 name=$3 mem=$4; shift 4
  local out="$BASE/out/$name"
  local runner="$BASE/out/$name.sh"
  mkdir -p "$out" "$BASE/scratch/$name" "$ROOT/results/server/${DATE}_${exp}"
  # a generated runner script keeps all the quoting in one place (a tmux command string with
  # nested quotes silently dies before run_job.sh is ever reached)
  {
    echo "#!/usr/bin/env bash"
    echo "export AIM3_RESULTS=$(printf %q "$out")"
    echo "export AIM3_SCRATCH=$(printf %q "$BASE/scratch/$name")"
    echo "export AIM3_SAGE=$(printf %q "$AIM3_SAGE")"
    echo "export PATH=$(printf %q "$BASE/opt/bin"):\$PATH"
    echo "cd $(printf %q "$ROOT")"
    printf 'bash server/jobs/run_job.sh %q %q %q --' "$exp" "$name" "$mem"
    printf ' %q' "$@"
    printf '\n'
    echo "echo JOB_DONE_rc=\$? >> $(printf %q "$ROOT/results/server/${DATE}_${exp}/$name.meta")"
  } > "$runner"
  chmod +x "$runner"
  tmux kill-session -t "$sess" 2>/dev/null
  tmux new-session -d -s "$sess" "bash $runner"
  echo "launched $sess -> results/server/${DATE}_${exp}/$name.log  (runner $runner)"
}

stage=$1; shift
case "$stage" in
  e03)
    for n in "$@"; do
      launch "e03n$n" E03 "E03_sweep_n$n" 48 \
        env AIM3_BACKEND=sage E03_BUDGET="${E03_BUDGET:-14400}" E03_EXTRA="${E03_EXTRA:-}" \
        bash experiments/E03_syzygy_and_degree_fall/sweep.sh "$n"
    done ;;
  e03d5)
    for n in "$@"; do
      launch "e03d5n$n" E03 "E03_D5_n$n" 120 \
        python3 server/jobs/e03_d5.py "$n" --backend sage
    done ;;
  gate1)
    ell=$1; shift
    for n in "$@"; do
      launch "g1l${ell}n$n" E04 "E04_gate1_l${ell}_n$n" 120 \
        python3 experiments/E04_full_system_hilbert/full_system_falls.py "$n" -l "$ell" \
          --seeds "${SEEDS:-1,2,3}" --systems S_quad,P1only,random --per-block \
          --backend sage --budget "${BUDGET:-86400}"
    done ;;
  e12)
    for n in "$@"; do
      launch "e12n$n" E12 "E12_first_layer_n$n" 120 \
        python3 experiments/E12_cross_sbox_relations/first_layer_ideal.py "$n" --D4 \
          --backend sage --budget "${BUDGET:-86400}" ${E12_EXTRA:-}
    done ;;
  e06)
    for n in "$@"; do
      launch "e06n$n" E06 "E06_kernel_n$n" 120 \
        python3 experiments/E06_kernel_bilinear_model/kernel_bilinear_toy.py "$n" "${DMAX:-4}" \
          --backend sage --budget "${BUDGET:-86400}"
    done ;;
  f4export)
    D="$BASE/out/E04_systems"
    mkdir -p "$D"
    export AIM3_RESULTS="$BASE/out/f4export" AIM3_SCRATCH="$BASE/scratch/f4export"
    mkdir -p "$AIM3_RESULTS" "$AIM3_SCRATCH"
    cd "$ROOT"
    python3 experiments/E04_full_system_hilbert/export_systems.py --out "$D" \
        --ns "${NS:-7,9,11,13,15}" --ells "${ELLS:-1,2}" --seeds "${SEEDS:-1,2}"
    for n in ${E06NS:-7 9 11 13}; do
      python3 experiments/E06_kernel_bilinear_model/kernel_bilinear_toy.py "$n" 2 --export "$D" \
        | grep -E "export|params"
    done
    echo "systems in $D:"; ls "$D"/*.json | wc -l ;;
  f4)
    tag=${1:-}
    launch "f4${tag:-all}" E04 "E04_f4_degrees${tag:+_$tag}" 120 \
      python3 experiments/E04_full_system_hilbert/run_f4_degrees.py \
        --systems "$BASE/out/E04_systems" ${tag:+--tags "$tag"} \
        --tools "${TOOLS:-msolve,magma}" \
        --orders "${ORDERS:-grevlex,blockPT,blockT,blockU,blockW}" \
        --threads "${THREADS:-32}" --timeout "${TIMEOUT:-3600}" \
        --csv E04_f4_degrees.csv --logdir "$BASE/out/f4logs" ;;
  e05)
    launch e05 E05 "E05_resultants" 150 \
      python3 experiments/E05_two_valued_inverse/run_resultants.py \
        --ns "${NS:-8,9,10,11,12,13,14,15,16,17,18,19,20}" --seed 1 \
        --timeout "${TIMEOUT:-7200}" --mem-gb "${MEMGB:-120}" \
        --workdir "$BASE/scratch/E05" --csv E05_server_resultants.csv --factor ;;
  cmd)
    name=$1; mem=$2; shift 2
    launch "$name" "${EXP:-MISC}" "$name" "$mem" "$@" ;;
  status)
    tmux ls 2>/dev/null || echo "(no tmux sessions)"
    for f in "$ROOT"/results/server/*/*.log; do
      [ -e "$f" ] || continue
      echo "--- $f"; tail -3 "$f"
    done ;;
  *) echo "unknown stage $stage"; exit 2 ;;
esac
