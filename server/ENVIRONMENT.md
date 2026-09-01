# Server environment — Tier B/C toolchain, verified 2026-08-31

Bootstrapped by `bash server/setup_env.sh` (idempotent; stages `micromamba sage msolve parallel
check`).  The account has **no passwordless sudo**, so nothing is installed system-wide: every tool
built for this project lives under `$AIM3_OPT = /home/luck/rxw/AIM3_remote/opt`.

## Access and layout

| what | where |
|---|---|
| SSH | `ssh -F ~/.ssh/config_server server` → `luck@192.168.111.172` |
| repository checkout | `/home/luck/rxw/AIM3_remote/AIM3` (synced from the WSL working copy with `rsync -az --delete --exclude results/server/`) |
| tools built for this batch | `/home/luck/rxw/AIM3_remote/opt/{bin,lib,include,envs,src,logs}` |
| CSV output of a job (`$AIM3_RESULTS`) | `/home/luck/rxw/AIM3_remote/out/<job>/` — deliberately **outside** the checkout, so a server run never rewrites a committed Tier-A table |
| scratch matrices (`$AIM3_SCRATCH`) | `/home/luck/rxw/AIM3_remote/scratch/<job>/` — exported Macaulay matrices, never committed |
| job logs | `results/server/<YYYYMMDD>_<EXP>/<name>.{log,time,meta}` inside the checkout |

Nothing is written outside `/home/luck/rxw/AIM3_remote/`; the sibling directories under
`/home/luck/rxw/` belong to unrelated projects.

## Hardware

```
host    luck-Precision-7920-Tower
os      Ubuntu 22.04.5 LTS, kernel 6.8.0-138-generic
cpu     Intel(R) Xeon(R) Gold 6230R @ 2.10GHz, 104 threads (2 sockets x 26 cores x 2)
memory  251 GB
disk    /  11 TB, 8.5 TB free
gpu     1x NVIDIA GeForce RTX 3080 (10 GB)   [unused by this batch]
```

## Software

| tool | version | provenance |
|---|---|---|
| **Magma** | **V2.24-2** | pre-existing, `/home/luck/bin/magma` → `/home/luck/magma-2018`, licence file present |
| **SageMath** | **10.7** (2025-08-09) | conda-forge, `$AIM3_OPT/envs/sage10/bin/sage`; the *system* Sage is 9.5, below the `server/README.md` requirement |
| **msolve** | **0.10.1** | source build from `github.com/algebraic-solving/msolve` @ `1e3af01`, `$AIM3_OPT/bin/msolve` |
| **GNU parallel** | 20260822 | source build, `$AIM3_OPT/bin/parallel` |
| FLINT | 3.3.1 (conda, in the sage10 env) | msolve HEAD needs FLINT ≥ 3 |
| M4RI | libm4ri 0.0.20200125 (system) and the copy inside Sage 10.7 | used through Sage's dense `GF(2)` matrices |
| Python | 3.10.12 (system), 3.11 (sage10 env) | the experiment scripts are zero-dependency and run under the system Python |
| gcc / make / cmake / tmux / git | 11.4.0 / GNU / 3.22 / 3.2a / 2.54.0 | system |
| `jq` | **not installed** | log parsing is done in Python |

### Magma notes that cost time to discover

1. `magma` is only on the **login-shell** `PATH`; a plain `ssh server '<cmd>'` cannot find it.
   Always `ssh -F ~/.ssh/config_server server 'bash -lc "..."'`.
2. A non-interactive run must close stdin — `magma -b script.m < /dev/null` — or it waits forever.
3. `GetVersion()` returns **three** values; `printf "%o.%o-%o", GetVersion()` raises
   *"Not enough arguments for % conversions"*.  Assign first: `a,b,c := GetVersion();`.
4. `load` accepts a **string literal only**, so a system file cannot be named on the command line.
   `run_f4_degrees.py` and `run_resultants.py` generate a three-line wrapper `.m` instead.
5. `GetMemoryLimit()` is `0` (unlimited) — the per-tier memory caps of `server/README.md` are
   enforced by `server/jobs/run_job.sh` with `ulimit -v`, not by Magma.
6. `SetNthreads(k)` works and **imposes no upper limit of its own** — it accepts 128 as readily as
   1, above the machine's 104 hardware threads (measured: 1, 2, 4, 8, 16, 26, 32, 52, 64, 104, 128
   all accepted, `GetNthreads()` echoing each).  So the ceiling is the hardware, not Magma, and a
   job that asks for more than 104 will oversubscribe silently.  Magma's F4 uses the threads for
   linear algebra only; the jobs in this batch ran single-threaded Magma with 8 msolve threads,
   because the batch is a farm of many small jobs rather than one large one.

### Sage notes

* `matrix(GF(2), r, c)` is M4RI; `M.echelonize()` then `M.pivots()` gives rank **and** the pivot
  columns in one pass, which is what `sage/rank_backend.sage` needs.
* `.sage` files are pre-parsed: integer literals become Sage `Integer`s.  `sys.exit(0)` therefore
  exits with status 1, and `round(x, 3)` returns a `RealDoubleElement` that `json` refuses to
  serialise.  Both are wrapped in `int()` / `float()` in the backend.
* **`^` is exponentiation, not XOR.**  The preparser rewrites `a ^ b` to `a ** b` and `a ^= b` to
  `a **= b` inside `.sage` files, so every bitwise XOR must be written `_xor(a, b)`
  (`from operator import xor as _xor`).  This bit `sage/ore_scale.sage`: a Möbius transform written
  with `arr[i] ^= arr[i ^ step]` raised `MemoryError` (it was computing a tower of powers).  A
  `MemoryError` is the *lucky* case — `s ^= f0` on small values preparses to `s **= f0` and returns
  a wrong number silently.  Grep every new `.sage` file for `^` before trusting its output.
* `ulimit -v 32G` does **not** break Sage start-up (checked), so the memory caps are usable.
* `OrePolynomialRing(K, K.frobenius_endomorphism())` works (it is marked experimental and warns), and
  `f.left_xlcm(g)` returns `(l, u, v)` with `u*f == v*g == l`.  It is **50–1000× slower per `lclm`**
  than the pure-Python `src/aim3/ore.py`: use it as a correctness reference, not as a workhorse.
  The same ring over the rational function field `K(T)` (twist `T ↦ T²`, built with
  `LF.hom([LF(T)^2], base_map=K.frobenius_endomorphism())`) is accepted but does not finish an
  `lclm` of two degree-14 elements in 6 h at n = 8 — the degree in `T` doubles at every Euclidean
  step.

### msolve notes

* Default behaviour computes a rational parametrisation, which **re-runs F4 about 20 times** with
  different random choices.  The step degrees are identical every time, but the step count and the
  wall clock then describe 20 runs, not one — `run_f4_degrees.py` passes `-g 2` (reduced Gröbner
  basis only).
* `-v 2` prints a fixed-width table (`deg sel pairs mat density new-data time`); the summary block
  further down contains a `max. matrix data R x C` line that must not be parsed as a step.

## Verification performed

```
$ bash scripts/smoke_test.sh                                   # see below
$ python3 tests/run_tests.py                                   # 6/6 before this batch, 8/8 after
$ python3 server/jobs/verify_backend.py                        # 76/76, BACKEND ALIGNED
```

The two added tests are `test_export_matches_boolmac` (the exported file describes the same matrix
`boolmac.py` streams, and the Sage backend's column reversal reproduces its `rank`/`dim_low`) and
`test_legacy_field_table_extended` (the added moduli are irreducible and the shared entries are
unchanged).

`verify_backend.py` is the gate for task 2: it compares the pure-Python engine and the Sage/M4RI
engine on the same instances (E03 n=11 e=3 and n=13 e=7, E12 n=11 `--D4`, E06 n=7 D=4), for
`rows, cols, rank, syz, dim(span ∩ deg≤D−1)`, the derived `falls_D4` / `cross_excess`, **and** the
point-set ideal dimension `dim I(V)_{≤3}`; and it re-checks both against the committed Tier-A CSV
rows.  All 76 comparisons agreed, so `--backend sage` is cleared for larger n.

At n = 7, ℓ = 2 (`S_quad`) msolve and Magma independently report the same F4 profile —
`d_solv = 3`, 5 steps — which is the cross-tool validation asked for in task 6.

One extension was needed to reach n > 16 at all: `legacy/aim3_lab.CONWAY` stops at n = 16, so
`src/aim3/legacy_fields.py` copies the n = 17..24 moduli from `src/aim3/gf2n.TOY_MODULI` into it
**after verifying each one is irreducible** (Rabin's test).  `legacy/` itself is untouched
(CLAUDE.md rule 3), and the moduli shared with the legacy table are asserted to be identical, so
Tier-A and server instances describe the same fields.

## Reproducing

```bash
ssh -F ~/.ssh/config_server server 'bash -lc "
  cd /home/luck/rxw/AIM3_remote/AIM3 && bash server/setup_env.sh"'
```

Per-job environment (set by `server/jobs/batch1.sh`):

```bash
export AIM3_SAGE=/home/luck/rxw/AIM3_remote/opt/envs/sage10/bin/sage
export AIM3_RESULTS=/home/luck/rxw/AIM3_remote/out/<job>
export AIM3_SCRATCH=/home/luck/rxw/AIM3_remote/scratch/<job>
export PATH=/home/luck/rxw/AIM3_remote/opt/bin:$PATH     # msolve, parallel
```
