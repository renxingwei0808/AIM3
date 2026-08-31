# Server plan (2× Intel Xeon Gold 6230R, 26 cores/52 threads each, 2.1 GHz, 256 GB RAM) — driven by local Claude Code over SSH

## Software to install / verify
* Magma (licensed, already available) — check `SetNthreads(52)`; F4 uses threads for linear algebra.
* SageMath ≥ 10 (conda-forge `sage` or the Sage binary) — Ore polynomials, msolve bindings, M4RI via `sage.matrix`.
* msolve (https://msolve.lip6.fr) — fast F4/FGLM over F₂ with `-t <threads>`, prints step degrees with `-v 2`.
* M4RI/M4RIE (dense GF(2) linear algebra) — used through Sage or a small C driver for Macaulay ranks up to ~10⁶ × 10⁶ bits.
* GNU parallel, tmux, Python 3.11+ (numpy optional).

## Job conventions
* One instance = one (experiment, n, exps, seed, distribution) tuple; one process each; farm with `jobs/run_matrix.sh`.
* Every job writes `results/server/<date>_<exp>/<params>.json` + `.log`; a watchdog kills jobs over the RAM cap (`ulimit -v`).
* Memory caps per job by tier: A 4 GB, B 32 GB, C up to 200 GB (one job at a time).
* Never run Magma jobs that need > 1 day without a checkpoint of the intermediate basis (`SetVerbose` logs are the checkpoint for degree data).

## What runs here (priority order)
1. E03 D=4 Macaulay at n = 16..24 and D=5 at n = 11..14 (M4RI): confirms the 2n-fall theorem beyond Python reach.
2. E04 F4 step-degree logs for toy AIM3, n = 8..16, ℓ = 1,2,3, 20 seeds × 3 matrix distributions (Magma + msolve cross-check).
3. E06 / E10 block-order F4 on the kernel model vs S_quad (Magma).
4. E05 resultants of the two-valued-inverse reduced system, n = 8..24 (Magma).
5. E07 skew-polynomial LCLM / approximant-basis degree growth, n = 8..20 (Sage).
6. E09 weak-iv statistics on the real parameter sets (Python, many cores).

See `SIZING.md` for matrix-size and memory estimates.
