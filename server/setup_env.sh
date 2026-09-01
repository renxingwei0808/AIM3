#!/usr/bin/env bash
# Idempotent Tier-B/C toolchain bootstrap for the AIM3 compute server.
#
#   bash server/setup_env.sh [stage ...]      stages: micromamba sage msolve parallel check
#   (no argument = every stage)
#
# Everything is installed under $AIM3_OPT (default /home/luck/rxw/AIM3_remote/opt) because the
# account has no passwordless sudo; nothing outside that prefix and ~/.local is touched.
# Each stage logs to $AIM3_OPT/logs/<stage>.log and is skipped when its binary already exists.
set -uo pipefail

OPT="${AIM3_OPT:-/home/luck/rxw/AIM3_remote/opt}"
BIN="$OPT/bin"
LOGS="$OPT/logs"
SRC="$OPT/src"
JOBS="${BUILD_JOBS:-32}"
mkdir -p "$BIN" "$LOGS" "$SRC"

say() { printf '\n=== [%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }

# ---------------------------------------------------------------- micromamba
stage_micromamba() {
  if [ -x "$BIN/micromamba" ]; then say "micromamba already present"; return 0; fi
  say "installing micromamba into $BIN"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "$OPT" bin/micromamba >>"$LOGS/micromamba.log" 2>&1
  [ -x "$BIN/micromamba" ]
}

# ---------------------------------------------------------------- Sage >= 10
# The system Sage is 9.5; server/README.md asks for >= 10 (Ore polynomials, msolve bindings).
stage_sage() {
  if [ -x "$OPT/envs/sage10/bin/sage" ]; then say "sage10 env already present"; return 0; fi
  say "creating conda-forge env sage10 (this downloads ~2 GB)"
  export MAMBA_ROOT_PREFIX="$OPT"
  "$BIN/micromamba" create -y -p "$OPT/envs/sage10" -c conda-forge \
      "sage>=10" python=3.11 >>"$LOGS/sage.log" 2>&1
}

# ---------------------------------------------------------------- msolve
# Source build (server/README.md). msolve needs GMP + FLINT; the system FLINT is 2.8.4, which is
# older than msolve's minimum, so we build FLINT 2.9.0 into the same prefix first when needed.
stage_flint() {
  if [ -f "$OPT/lib/libflint.so" ]; then say "local FLINT already built"; return 0; fi
  say "building FLINT 2.9.0 (system 2.8.4 is too old for msolve)"
  cd "$SRC" || return 1
  [ -d flint2 ] || git clone --depth 1 -b v2.9.0 https://github.com/flintlib/flint2.git \
      >>"$LOGS/flint.log" 2>&1
  cd flint2 || return 1
  ./configure --prefix="$OPT" --disable-static >>"$LOGS/flint.log" 2>&1 &&
  make -j"$JOBS" >>"$LOGS/flint.log" 2>&1 &&
  make install >>"$LOGS/flint.log" 2>&1
}

stage_msolve() {
  if [ -x "$BIN/msolve" ]; then say "msolve already present"; return 0; fi
  say "building msolve from source"
  cd "$SRC" || return 1
  [ -d msolve ] || git clone --depth 1 https://github.com/algebraic-solving/msolve.git \
      >>"$LOGS/msolve.log" 2>&1
  cd msolve || return 1
  ./autogen.sh >>"$LOGS/msolve.log" 2>&1
  # msolve HEAD needs FLINT >= 3 (FLINT_FORCE_INLINE); the system FLINT is 2.8.4 and the
  # locally built 2.9.0 is still too old, so build against the conda-forge FLINT 3 that came
  # with the sage10 env.  Everything stays inside $OPT.
  local FL="$OPT/envs/sage10"
  if [ ! -f "$FL/include/flint/flint.h" ]; then
    say "no FLINT 3 headers under $FL -- run the sage stage first"; return 1
  fi
  make distclean >>"$LOGS/msolve.log" 2>&1
  CPPFLAGS="-I$FL/include" LDFLAGS="-L$FL/lib -Wl,-rpath,$FL/lib" \
    ./configure --prefix="$OPT" >>"$LOGS/msolve.log" 2>&1 &&
  make -j"$JOBS" >>"$LOGS/msolve.log" 2>&1 &&
  make install >>"$LOGS/msolve.log" 2>&1
  [ -x "$BIN/msolve" ]
}

# ---------------------------------------------------------------- GNU parallel
stage_parallel() {
  if [ -x "$BIN/parallel" ]; then say "GNU parallel already present"; return 0; fi
  say "building GNU parallel from source"
  cd "$SRC" || return 1
  curl -Ls -o parallel-latest.tar.bz2 https://ftpmirror.gnu.org/parallel/parallel-latest.tar.bz2 \
      >>"$LOGS/parallel.log" 2>&1
  tar xjf parallel-latest.tar.bz2 >>"$LOGS/parallel.log" 2>&1
  cd "$(ls -d parallel-2* | tail -1)" || return 1
  ./configure --prefix="$OPT" >>"$LOGS/parallel.log" 2>&1 &&
  make -j"$JOBS" >>"$LOGS/parallel.log" 2>&1 &&
  make install >>"$LOGS/parallel.log" 2>&1
  mkdir -p "$HOME/.parallel" && touch "$HOME/.parallel/will-cite"
}

# ---------------------------------------------------------------- report
stage_check() {
  say "versions"
  export PATH="$BIN:$PATH"
  printf 'host          : %s\n' "$(hostname)"
  printf 'os            : %s\n' "$(. /etc/os-release; echo "$PRETTY_NAME")"
  printf 'kernel        : %s\n' "$(uname -r)"
  printf 'cpu           : %s x %s threads\n' "$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)" "$(nproc)"
  printf 'memory        : %s\n' "$(free -g | awk '/^Mem:/{print $2" GB"}')"
  printf 'magma         : %s\n' "$(echo 'a,b,c := GetVersion(); printf "V%o.%o-%o", a,b,c; quit;' > /tmp/_v.m && magma -b /tmp/_v.m </dev/null 2>&1 | tail -1)"
  printf 'magma path    : %s\n' "$(command -v magma)"
  printf 'python3       : %s\n' "$(python3 --version 2>&1)"
  printf 'sage (system) : %s\n' "$(sage --version 2>&1 | head -1)"
  printf 'sage (conda)  : %s\n' "$( [ -x "$OPT/envs/sage10/bin/sage" ] && "$OPT/envs/sage10/bin/sage" --version 2>&1 | head -1 || echo 'NOT INSTALLED')"
  printf 'msolve        : %s\n' "$(have msolve && msolve -h 2>&1 | head -1 || echo 'NOT INSTALLED')"
  printf 'parallel      : %s\n' "$(have parallel && parallel --version 2>&1 | head -1 || echo 'NOT INSTALLED')"
  printf 'gcc           : %s\n' "$(gcc --version | head -1)"
  printf 'libm4ri       : %s\n' "$(ls /usr/lib/x86_64-linux-gnu/libm4ri*.so 2>/dev/null | head -1)"
  printf 'libflint sys  : %s\n' "$(ls /usr/lib/x86_64-linux-gnu/libflint-*.so 2>/dev/null | head -1)"
  printf 'libflint local: %s\n' "$(ls "$OPT"/lib/libflint.so* 2>/dev/null | head -1 || echo '-')"
  printf 'tmux          : %s\n' "$(tmux -V)"
  printf 'git           : %s\n' "$(git --version)"
}

stages=("$@")
[ ${#stages[@]} -eq 0 ] && stages=(micromamba sage msolve parallel check)
rc=0
for s in "${stages[@]}"; do
  if "stage_$s"; then say "stage $s OK"; else say "stage $s FAILED (see $LOGS/$s.log)"; rc=1; fi
done
say "setup finished with rc=$rc"
exit $rc
