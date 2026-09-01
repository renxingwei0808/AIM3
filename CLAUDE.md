# Instructions for Claude (Cowork / Claude Code) working in this repository

## What this project is
Independent academic cryptanalysis of AIM3 (AIMer v3.0). Read `docs/ROADMAP.md` first, then the README of the experiment you
work on. The roadmap fixes the six technique tracks T1–T6 and the experiments E00–E12; do not invent new top-level directories.

## Ground rules
1. **Never overclaim.** Use the status tags [PROVED]/[VERIFIED n≤…]/[CONJ]/[OPEN]. A computation is [VERIFIED] only for the
   parameters actually run. "Attack" means a complete procedure with an accounted cost against B1/B2/B3 (see
   `docs/COMPLEXITY_ACCOUNTING.md`); anything else is a "model" or an "observation".
2. **Reproducibility.** Every script takes its parameters on the command line, prints them, and appends a CSV row to `results/`
   with `tier, host, commit, seconds`. Fix seeds. Keep control groups (random systems of the same shape; unit-LU vs GL_n).
3. **Do not edit `legacy/`.** Extend `src/aim3/` or add scripts inside the experiment directories. Import via `experiments/_common.py`.
4. **Toy parameters are not the real ones.** Extrapolate only under the rule in COMPLEXITY_ACCOUNTING §4. Never report a
   timed-out run as a value; report it as a lower bound.
5. **Spec references.** When stating a fact about AIM3 cite the page of the v3.0 specification (see `docs/FACT_SHEET.md`).
   Third-party complexity constants must be checked against the original PDFs before being quoted.
6. **Encoding conventions** of the reference implementation are unpinned until E00 is done; do not quote full-size iv-dependent numbers before that.

## Tier A (Cowork, this laptop / web)
Pure Python, n ≤ 16, minutes. Suitable: E01, E02, E03 (D≤4, n≤13), E05 statistics, E06 (n≤9), E08, E11, E12 (deg 3, n≤16),
all writing tasks (paper/, docs/). Run `bash scripts/smoke_test.sh` before committing.

## Tiers B/C (Claude Code on the local WSL machine → remote server over SSH)
* Server: 2× Xeon Gold 6230R (52 cores/104 threads), 256 GB RAM; Magma licensed. See `server/README.md` and `server/SIZING.md`.
* Work in `tmux`; launch long jobs with `nohup … > results/server/<date>_<exp>/<params>.log 2>&1 &`; poll, do not block.
* Respect memory caps (`ulimit -v`); one tier-C job at a time above 50 GB.
* Magma/Sage files under `experiments/*/magma|sage` are **untested skeletons**: validate on n ≤ 13 against the Python CSVs
  (`results/E03_sbox_degree_fall.csv` etc.) before trusting larger-n output.
* Copy back only summary CSVs and logs; large matrices stay on the server.

## Writing
Paper in English (`paper/`), roadmap and internal notes may be in Chinese. Every paper claim gets a row in `paper/claims_ledger.md`.

---

## Environment — verified 2026-08-31 (checked by Claude Code, WSL laptop)

### Language of interaction
**All final deliverables, summaries and questions addressed to the user are written in Chinese (中文).**
Code, identifiers, commit messages, CSV headers and the paper in `paper/` stay in English.

### GitHub (upload target)
* Repo `renxingwei0808/AIM3` — **private**, default branch `main`.
* Permission of the authenticated account (`gh` CLI, user `renxingwei0808`):
  `admin=true, maintain=true, push=true, pull=true, triage=true` → **full admin, push allowed**.
* Token scopes: `repo`, `workflow`, `read:org`, `gist`. Protocol: HTTPS via `gh`.
* Local working copy: `/home/star/AIM3_local` (clone of `origin/main`).
* **Upload policy: never push on my own initiative.** Commit locally in `/home/star/AIM3_local` and push
  only when the task prompt explicitly asks for it. When it does: branch + PR for anything larger than a
  doc tweak; never force-push `main`.

### Remote compute server
* Access: `ssh -F ~/.ssh/config_server server` → `luck@192.168.111.172`
  (host `luck-Precision-7920-Tower`, key `~/.ssh/server_ed25519`). Connection verified.
* Ubuntu 22.04.5 LTS, **104 threads**, **251 GB RAM**, `/` = 11 TB with 8.5 TB free, 1× RTX 3080 (10 GB).
* **All code and runs live in `/home/luck/rxw/AIM3_remote/`.** Do not write outside it
  (`/home/luck/rxw/` holds unrelated projects: `yuh_*`, `research_MDS`, `benchmarks`, …).
* Available: `python3` 3.10.12, `sage` 9.5, `gap`, `gp`, `git`, `tmux` 3.2a, `gcc/make/cmake`, `mpirun`,
  system `libm4ri` and `libflint`.
* **Missing** (install before the jobs in `server/README.md` that need them): `msolve`, GNU `parallel`,
  Sage is 9.5 not ≥ 10.

### Magma on the server — available ✅
* `/home/luck/bin/magma` (wrapper) → `/home/luck/magma-2018`, **Magma V2.24-2**, licence file present.
* Verified working: `GF(2,127)`, `PolynomialRing(..., "grevlex")`, `GroebnerBasis`, `HilbertSeries`,
  `SetNthreads(4)` / `GetNthreads()` (multi-threaded F4 linear algebra is usable — scale up to 52 per `server/README.md`).
* **Three traps found while testing — always follow these:**
  1. `magma` is only on the **login-shell** `PATH` (`~/.bashrc`). A plain `ssh server '<cmd>'` will not find it.
     Always wrap: `ssh -F ~/.ssh/config_server server 'bash -lc "..."'`.
  2. Non-interactive runs must close stdin or Magma hangs forever:
     `magma -b script.m < /dev/null` (use `-b` for batch, no banner).
  3. `GetVersion()` returns three values — `a,b,c := GetVersion();` before printing, otherwise
     `printf` raises "Not enough arguments for % conversions".
* `GetMemoryLimit()` returns `0` (unlimited) → the memory caps in `server/README.md` must be enforced by us
  with `ulimit -v` in the launch wrapper, not by Magma.

### Standard remote job recipe
```bash
ssh -F ~/.ssh/config_server server 'bash -lc "
  cd /home/luck/rxw/AIM3_remote/<exp> &&
  mkdir -p results &&
  ulimit -v <kb> &&
  nohup magma -b <script>.m </dev/null > results/<params>.log 2>&1 &
  echo started \$!"'
```
Poll the log; never block on a long job. Copy back only summary CSV/JSON + logs.
