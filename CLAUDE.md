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
