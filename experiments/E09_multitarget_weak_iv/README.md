# E09 — Multi-target, weak-iv and reused-key scenarios (track T6)

AIM3 explicitly targets multi-target one-wayness (spec §3.1, §6.1). Questions:
1. **Weak iv classes.** Is there a property Π of the generated (A_j, b_j) (e.g. A₁A₀⁻¹ Frobenius-sparse, common invariant
   subspace, low-rank part of the twisted-bilinear operator) that makes any track-T2/T3 attack cheaper by a factor 2^s, with
   density ρ(Π) among iv's? Amortised cost over N targets: min over Π of 2^{-s}·cost / (N·ρ(Π)) once N·ρ ≥ 1.
2. **Reused key with different iv** (misuse model, as in the Magic Pot paper): q pairs (iv_j, ct_j) with the same pt give
   3nq equations in n + 2nq variables (P1). Compute d_solv for q = 2,3 on toy n (E04 machinery) — does the extra pt-coupling help?
3. **Zero-input rejection leakage**: valid keys avoid ℓ+1 affine hyperplanes' points; quantify (negligible, but state it).
Tier A/B. Deliverable: a scenario table with amortised formulas; results in the single-key and multi-target models are reported separately.

---

## Reused key, q ≤ 6 (Tier-A batch 3): `reused_key_toy.py`

```
python3 reused_key_toy.py 9  -l 2 -q Q --seeds 1 --exps 2,4,5   # Q = 2..6
python3 reused_key_toy.py 11 -l 2 -q Q --seeds 1 --exps 2,3,4   # Q = 2..5
python3 reused_key_toy.py 11 -l 2 -q Q --seeds 1,2 --exps 3,4,10
python3 reused_key_toy.py 9  -l 2 --single-only --exps 4,5,8    # exponent-set comparison at q=1
```
Rows go to `results/E09_theory_reused_key.csv`.  **Every run below has saturation ≤ 0.29**, so the D=3 numbers are
informative (the batch-3 brief's threshold was 0.9); the n = 7 rows already in the CSV are at 0.75–0.99 and are excluded
from all statements here.  A D=3 pass that hits its `--budget` now returns "NOT MEASURED" instead of a row: a censored rank
is a lower bound, and a lower bound on the rank says nothing about a syzygy count (CLAUDE.md ground rule 4).

Notation: `syz(q)` is the D=3 syzygy count of the q-instance system, `syz(1)` that of one instance.  "Excess" means
`syz(q) − q·syz(1)`, i.e. the coupling beyond a direct sum of q independent instances.

### The C(q,2)·n law holds — for exponent sets whose single instance is generic

| n | exps | syz(1) | (ℓ+1)·8n | q | syz(q) | excess | C(q,2)·n | saturation |
|---|------|--------|----------|---|--------|--------|----------|-----------|
| 9 | 2,4,5 | 216 | 216 | 2 | 441 | 9 | 9 | 0.918 |
| 9 | 2,4,5 | 216 | 216 | 3 | 675 | 27 | 27 | 0.147 |
| 9 | 2,4,5 | 216 | 216 | 4 | 918 | 54 | 54 | 0.133 |
| 9 | 2,4,5 | 216 | 216 | 5 | 1170 | 90 | 90 | 0.120 |
| 9 | 2,4,5 | 216 | 216 | 6 | 1431 | 135 | 135 | 0.110 |
| 11 | 2,3,4 | 264 | 264 | 2 | 539 | 11 | 11 | 0.136 |
| 11 | 2,3,4 | 264 | 264 | 3 | 825 | 33 | 33 | 0.110 |
| 11 | 2,3,4 | 264 | 264 | 4 | 1122 | 66 | 66 | 0.096 |
| 11 | 2,3,4 | 264 | 264 | 5 | 1430 | 110 | 110 | 0.089 |

**[VERIFIED n ∈ {9,11}, ℓ=2, q ≤ 6]** `syz(q) = q·syz(1) + C(q,2)·n` — one extra syzygy per *pair* of instances per
coordinate, exactly what "the q instances share pt and nothing else" predicts.  The count is seed-independent (n = 11,
seeds 1 and 2 give identical numbers).

### The law is not universal: e = n−1 changes it

Running the same experiment with an all-generic exponent set (`COMPLEXITY_ACCOUNTING.md` rule 5) at n = 11 turned up a
different regime, and the cause is a single exponent:

| n | exps | contains e = n−1 | syz(1) | (ℓ+1)·8n | new elements at D=3, q=1 |
|---|------|------------------|--------|----------|--------------------------|
| 9 | 2,4,5 / 2,5,7 | no | 216 | 216 | 1 |
| 9 | 2,4,8 / 4,5,8 / 5,7,8 | yes | 217 | 216 | 27 |
| 11 | 2,3,4 / 3,4,5 | no | 264 | 264 | 0 / 1 |
| 11 | 3,4,10 (any order) | yes | 265 | 264 | 32 |
| 13 | 2,3,4 / 3,5,8 | no | 312 | 312 | 0 |
| 13 | 3,5,12 | yes | 313 | 312 | 38 |

**[VERIFIED n ∈ {9,11,13}, ℓ=2]** putting e = n−1 anywhere in the exponent set adds exactly **one** D=3 syzygy to the
single instance and exactly **3n−1** new degree-≤2 elements; every set without it sits on the (ℓ+1)·8n baseline.
This is a *different* statistic from the single-S-box D=4 fall count that `structures.degenerate_class` classifies — that
function calls e = n−1 generic, and it is right to, so this is a new class, for the full ℓ=2 system at D=3.  Note the
recommended AIM3 parameters never use e = n−1 (spec Table 2: n=128 → 3,7,11; n=192 → 11,23,47; n=256 → 3,7,11,19), so
this is a statement about toy instances, not about AIM3 itself; but it is a trap for anyone choosing toy exponents.

With e = n−1 present, the coupling excess is larger and follows a different law:

| n | exps | q | syz(q) | excess = syz(q) − q·syz(1) | n·[5·C(q,2) + (q−1)] |
|---|------|---|--------|----------------------------|----------------------|
| 11 | 3,4,10 | 2 | 596 | 66 | 66 |
| 11 | 3,4,10 | 3 | 982 | 187 | 187 |
| 11 | 3,4,10 | 4 | 1423 | 363 | 363 |
| 11 | 3,4,10 | 5 | 1919 | 594 | 594 |
| 11 | 3,4,10 | 6 | 2470 | 880 | 880 |
| 9 | 4,5,8 | 2 | 488 | 54 | 54 |
| 9 | 4,5,8 | 3 | 804 | 153 | 153 |
| 9 | 4,5,8 | 4 | 1165 | 297 | 297 |

**[VERIFIED n ∈ {9,11}, ℓ=2, q ≤ 6, 2 seeds at n=11]** the excess is `n·[5·C(q,2) + (q−1)]`, five times the generic
pair term plus a linear term.  **[OPEN]** why 5, and whether the coefficient depends on n (only two values of n are
available at Tier A; the server can reach n = 13, 15, 17).

### What this does and does not say

Reused-key coupling exists but is **quadratically small**: C(q,2)·n extra syzygies against q·(ℓ+1)·8n present ones — a
relative gain of O(q/n), and the new elements at D=3 stay far below the 3nq equations added.  Nothing here is an attack;
d_solv is not measured (that needs D ≥ 4 and is a Tier-B/C job).  The saturation column is the guard: a run above ~0.9
tells you only that the D=3 row space filled the ambient space.
