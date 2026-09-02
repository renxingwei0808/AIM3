# Experiments still running on the server at the close of batch 2

Written 2026-09-01 18:20 PDT, when `claude/server-batch-2` was pushed.  Four jobs were left running,
in two groups.  Neither group changes any conclusion already in `REPORT_server_batch2.md`; each
either **resolves a stated open question** or **closes a stated coverage gap**.  This file exists so
that whoever receives the repository can decide whether to let them finish, and can do so without
re-deriving what they are for.

Everything runs in `tmux` on `luck@192.168.111.172` under `/home/luck/rxw/AIM3_remote/AIM3`, detached
from any interactive session.  They survive logout and will keep running until they finish, hit their
`--budget`, or are killed.

---

## Group A — why `falls_D5 = n(4n+7)` fails at n = 18

**Sessions** `B2_d5_n18_e7`, `B2_d5_n18_e11`.  Started 2026-09-01 18:07 PDT.
**Command** `python3 server/jobs/e03_d5.py 18 --exps {7,11} --backend sage --budget 86400`
**Expected cost** ≈ 13 h wall each (the e = 5 run took 13 h 15 m), ≈ 25 GB peak RSS each, one core
each for the Python phases and M4RI's own threading for the echelon.  They run concurrently.

### The question

`REPORT_server_batch2.md` §2.2 establishes `falls_D5 = (#cubics)·2n − 18n + net` at n = 14 and 16,
which reproduces batch 1's closed forms `n(4n+7)` (generic exponent) and `n(4n−6)` ({−3} class).
The n = 18 run then measured **1152 = 64n**, where `4n+7` predicts **79n** — and the per-n count
*decreased* from 71 at n = 16, which no `4n + c` permits.  The run is complete, not censored, and its
random control behaves correctly (`falls_D5 = 0`).

Two readings are open, and one exponent cannot separate them:

1. **`n(4n+7)` is a fit that does not extend past n = 16.**  Then batch 1's closed forms are
   [VERIFIED n ≤ 16] and nothing more, and any use of them at larger n has to be withdrawn.
2. **(18, 5) is degenerate in a way `aim3.structures.degenerate_class` does not detect.**  n = 18 is
   the first size where that function returns an empty class for *every* admissible exponent, which
   is itself unusual and makes this reading worth taking seriously.

### What each outcome means

| e = 7 and e = 11 give | reading | consequence |
|---|---|---|
| both ≈ 64n (same as e = 5) | (1) — the form is wrong past n = 16 | narrow the ledger row to n ≤ 16; the D = 5 fall count needs a new form, fitted on ≥ 4 sizes |
| both ≈ 79n | (2) — e = 5 is special at n = 18 | keep the form; **`degenerate_class` has a missing class** and must be extended — that is a correctness bug in `src/aim3/structures.py` affecting every experiment that selects "generic" exponents |
| the two disagree with each other | a third class structure at n = 18 | neither reading; needs the remaining exponent (e = 13) and probably n = 20 |

Every outcome is informative, which is why both were launched rather than one.

### If they are stopped

Nothing already reported becomes wrong.  §2.2 already states the closed form as
`[VERIFIED n = 14, 16]` and explicitly not beyond, and §7 lists this as unfinished.  The cost of
stopping is that the ledger keeps a **known-unexplained** data point instead of a resolved one, and
that a possible bug in `degenerate_class` — which would matter well beyond E03 — stays unexamined.

---

## Group B — the n = 9 cell that task 3.1 asked for

**Sessions** `B2_multblk_n9_U`, `B2_multblk_n9_W`.  Started 2026-09-01 18:19 PDT.
**Command** `python3 experiments/E06_kernel_bilinear_model/kernel_bilinear_toy.py 9 5 --backend sage
--mult-block {U,W} --budget 86400`
**Expected cost** ≈ 1–3 h each, well under 60 GB.

### The question

The batch prompt's task 3.1 asks for block-restricted XL at **n = 7, 9, 11**.  The first pass
measured n = 7 (to D = 5) and n = 11 (to D = 4); n = 9 was missed.  These two jobs fill it.

### What the outcome means

Almost certainly nothing new: the finding is that the corank stays six orders of magnitude away from
the corank-1 linearisation threshold (909378 at n = 7 D = 5, 731600 at n = 11 D = 4), and U and W
track each other to within a handful of dimensions.  n = 9 landing between them is the expected
result.  The value is **coverage**, not discovery — a ledger row that says "n = 7, 9, 11 as asked"
rather than one that needs a footnote.

### If they are stopped

§3.1 already carries an explicit *Coverage shortfall* paragraph naming the missing cell, and the
ledger row in §8 is written for the sizes actually measured.  So stopping costs only the footnote.

---

## Operating them

```bash
ssh -F ~/.ssh/config_server server 'tmux ls'                       # what is still alive
ssh -F ~/.ssh/config_server server \
  'tail -5 /home/luck/rxw/AIM3_remote/AIM3/results/server/20260901_E0*/B2_*.log'
ssh -F ~/.ssh/config_server server 'tmux kill-session -t B2_d5_n18_e7'   # stop one
```

To take the results into the repository once they finish:

```bash
python3 server/jobs/collect.py                                     # pulls CSVs + logs back
python3 server/jobs/summarise_batch2.py --inject results/REPORT_server_batch2.md
```

then update `REPORT_server_batch2.md` §2.2 (group A) or §3.1 (group B) and the affected row of §8.
`collect.py` already has the merge entries for both — group A lands in `results/E03_server_D5.csv`
(dedupe key `(n, e)`), group B in `results/E06_server_multblock_{U,W}.csv`.

**Do not `pkill -f` on a pattern that also matches the ssh command line.**  It matches the remote
shell running the `pkill` itself and kills the session before the command reaches the target; this
cost ≈ 20 minutes of a batch-1 job and is recorded in `EXTRA_EXPERIMENTS.md` §10.  Use
`tmux kill-session -t <name>`.
