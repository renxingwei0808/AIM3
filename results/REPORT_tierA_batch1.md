# Tier-A batch 1 — report

**Date** 2026-08-31 · **Branch** `claude/tier-a-batch-1` · **Tier** A only (pure Python standard
library, one 4-core VM, 15 GB RAM) · **Nothing here needed Magma, Sage or the server.**

All numbers below come from CSV/MD files under `results/`, each row carrying
`tier, host, commit, seconds, censored`. **No run in this batch was censored**: every experiment
finished inside its budget, so no number in this report is a lower bound.

---

## 1. Anchor reproduction

Run before any change was made to the repository (host `vm`, commit `7249219`):

| anchor | required | measured | |
|---|---|---|---|
| E02 official Table 4 | 262.6 / 381.0 / 651.5 | 262.6 / 381.0 / 651.5 | ✅ |
| E02 Boolean series | 218.4 / 312.5 / 529.1 | 218.4 / 312.5 / 529.1 | ✅ |
| E03 (n=11, e=3) syz₃ | 88 | 88 | ✅ |
| E03 (n=11, e=3) falls_D4 | 22 | 22 | ✅ |
| E03 (n=11, e=3) cubics outside the ideal | 0 | 0 | ✅ |
| E03 (n=12, e=5) falls_D4 | 25 | 25 | ✅ |
| E06 (n=7, D=3) kernel syz | 0 | 0 | ✅ |
| E06 (n=7, D=3) dim(span ∩ deg≤2) | 84 | 84 | ✅ |
| E12 (n=14) cubic relations on (t₀,t₁) | 0 | 0 | ✅ |
| E12 (n=14) cubic relations on (pt,t₀,t₁) | 3444 | 3444 | ✅ |

All ten reproduce exactly. One pre-existing defect was found and fixed on the way:
`scripts/smoke_test.sh` exited 1 before printing `smoke test OK`, because `set -o pipefail` turned
the SIGPIPE of `… | head -2` into a script failure.

## 2. Infrastructure

* `experiments/_common.py::record()` / `record_table()` — append a CSV row and add
  `tier, host, commit, seconds, censored`; refuse to append to a file whose header does not match
  (so a schema change cannot silently corrupt a table). E02/E03/E04/E05/E06/E08/E12 all write
  through it.
* `tests/run_tests.py` — runs every `test_*` of `tests/test_core.py` without pytest; called at the
  end of `scripts/smoke_test.sh`. **6/6 pass** (two tests are new: `test_boolmac_matches_legacy`
  and `test_kat_conventions_pinned`).
* `src/aim3/boolmac.py` — streamed one-pass Boolean Macaulay statistics. It never stores the
  matrix and gets the rank and the top-degree rank from a single elimination (the column order puts
  degree-D monomials in the highest bit positions, so the pivots below that region are exactly a
  basis of rowspace ∩ deg≤D−1). It is numerically identical to `legacy/aim3_macaulay.macaulay`,
  which is asserted by a test. This is what made the E03 sweep cost minutes instead of hours
  (E03 at n=14 went from ~40 s to 15 s; the memory saving is what makes n=16 and the full system at
  n=11 possible at all). It also carries an optional `budget` that raises `MacaulayTimeout` with
  the partial counts as explicit lower bounds.
* `legacy/` was **not** modified. `docs/ROADMAP.md` gained only an appended status section.

---

## 3. Results, task by task

### 3.1 E03 — single-S-box degree falls, full table

**Run:** every admissible e (gcd(e,n)=1, 2 ≤ e ≤ n−2) for **n = 11, 12, 13, 14, 15, 16** — 36 runs,
0 censored, ≤ 81 s each. The batch asked for n = 11..14; n = 15 and 16 were added because they were
cheap and because they are what pins the rule below.
**Files:** `results/E03_sbox_degree_fall.csv`, **`results/E03_summary.md`** (full table + analysis).

1. **syz₃ = 8n in all 36 runs**, no exception. [VERIFIED n=11..16, every admissible e]
2. **`outside_ideal` = 0 in 33 of 36 runs** — the D=4 Macaulay matrix exhausts the cubic part of
   I(V_e), so there are no cubic generators outside ⟨P1,P2,P3⟩. The exceptions are exactly
   **e = 2** (n = 11, 13, 15), where **2n** cubics of the ideal are not reached at D=4.
3. **falls_D4 depends only on the unordered pair {e, e⁻¹ mod n}** (verified on every inverse pair
   present) and is 2n except on four classes:

   | class of e | falls_D4 | rows |
   |---|---|---|
   | {2, 2⁻¹} | **7n** | (11,2),(11,6),(13,2),(13,7),(15,2),(15,8) |
   | {−2, (−2)⁻¹} | 4n | (11,5),(11,9),(13,6),(13,11),(15,7),(15,13) |
   | {1+n/2} (n even) | 17n/2 | (12,7),(16,9) |
   | {−3, (−3)⁻¹} | 2n+1 | (11,7),(11,8),(13,4),(13,10),(14,9),(14,11),(16,5),(16,13) |
   | {n/2−1} (n even) | 2n+1 | (12,5),(16,7) |
   | generic | 2n | the other 12 |

   This reproduces **all 36 measured counts exactly**. It is [CONJ]: a fit to n ≤ 16, not a theorem.

**H2 — "do all non-degenerate e give 2n (or 2n+1)?"** Under the *previous* `degenerate_class`
(kmax = 4) the answer is **no in both directions**, and here are all the exceptions:
* 3 untagged e are not 2n/2n+1: **(11,9), (13,11), (15,13)**, all with 4n falls;
* 11 tagged e are exactly 2n: (11,3), (11,4), (13,3), (13,5), (13,8), (13,9), (14,3), (14,5),
  (15,4), (15,11), (16,11).

Under the calibrated rule the answer is yes for all 36 by construction — which is a restatement of
the fit, not independent evidence.

**The 2n vs 2n+1 discriminator — found, not [OPEN].** The extra fall appears iff
`3e ≡ −1 (mod n)` **or** `3e⁻¹ ≡ −1 (mod n)` **or** `e = n/2 − 1`. The **sign matters**:
`3e ≡ +1` gives the generic 2n — (11,4), (13,9), (14,5), (16,11) all satisfy it and give exactly 2n.
All the families the batch asked to check (k·e ≡ ±1, ±2 for k ≤ 6, for e and for e⁻¹; e vs n/2;
e² ≡ ±1) are tabulated in `results/E03_summary.md`; none of the others separates the two classes.
The batch's own list is reproduced: (12,5), (14,9), (16,5), (16,7) → 2n+1 and (11,3), (13,3),
(14,3), (14,5), (16,3) → 2n. ✅

**H3 — the degenerate classes: the earlier guess is corrected by the data.**
e = 2 gives **7n, not 9n**. "2e ≡ ±1 → 4n" is wrong as stated: 2e ≡ −1 gives 4n, but 2e ≡ +1 gives
**7n** (2e ≡ 1 means e = 2⁻¹, i.e. the class of e = 2). e = 1+n/2 gives exactly **17n/2 = 8.5n**
(the earlier "≈ 8.5n" is exact).

**Change to `aim3.structures.degenerate_class`, and why.** The kmax = 4 rule is wrong in both
directions (H2 above), so it was replaced by the class rule; `degenerate_class(n, e, legacy_kmax=4)`
still reproduces the old behaviour, and `predicted_falls_D4(n, e)` returns the predicted count.
Reason recorded in `experiments/E03_syzygy_and_degree_fall/README.md`. **All ten recommended AIM3
exponents remain generic (2n falls)** — now asserted inside `exponent_audit.py` and by a unit test.

### 3.2 E04 — the full system at D=4

**Run:** `full_system_falls.py` (new), n = 7, 9, 11, ℓ = 1 and 2, seeds 1 and 2, systems
`S_quad` / `P1only` / same-shape random quadratic control, plus a `--per-block` direct-sum control.
**File:** `results/E04_full_system_falls.csv`.

1. **D=3 is an exact direct sum.** syz₃(S_quad) = (ℓ+1)·8n on every non-saturated run:
   168 = 3·8·7, 216 = 3·8·9, 264 = 3·8·11 (ℓ=2); 144 = 2·8·9, 176 = 2·8·11 (ℓ=1). The random
   control has 0 syzygies at the same shape. [VERIFIED n=7,9,11 (ℓ=2) and n=9,11 (ℓ=1)]
2. **D=4 is not measurable at Tier A — the full system saturates.** S_quad has 3n(ℓ+1) equations in
   only (ℓ+1)n variables, so its D=4 Macaulay row space already fills the degree-≤3 space:
   `dim(span ∩ deg≤3)/C(nv,≤3)` is 0.9936–0.9997 in **every** run, and the **same-shape random
   control saturates identically** (n=11, ℓ=2: S_quad 6016/6018, random 6017/6018). `falls_D4` then
   measures the size of the monomial basis, not the algebra. Every such run is recorded with
   `verdict = indeterminate(saturated)`.
3. **The naive (ℓ+1)·2n criterion is not a direct-sum test even without saturation**, because the
   joint Macaulay matrix also multiplies by cross-block monomials (pt·t₁, t₀·t₁, …) absent from
   E03's 2n-variable ring. The `--per-block` control replaces it with
   `cross_excess = falls(all) − Σᵢ falls(S-box i alone in the same variables)`.
4. **A clean consistency check falls out of the per-block control at n=11** (the only per-block runs
   that are themselves non-saturated): the falls of a single S-box embedded in 3n variables are
   **77, 22, 22** for exponents 2, 3, 4 — i.e. exactly 7n, 2n, 2n, the E03 single-S-box values for
   those exponents. The extra multipliers add nothing for one S-box.

### 3.3 E12 — first layer vs the point-set ideal ⇒ **this is where Gate 1 is decided**

**Run:** `first_layer_ideal.py` (new) at n = 9, 11 (with D=4) and n = 14 (D=3); `cross_relations.py`
at n = 11, 14, 16; `explain_n16.py` (new).
**Files:** `results/E12_first_layer_ideal.csv`, `results/E12_cross_relations.csv`,
`results/E12_cross_relations_explained.csv`.

1. **n = 14, the ideal-theoretic comparison the batch asked for:**
   Mac₃ rank of the first-layer system (6n quadratics, 42 variables) = **3388**, syz₃ = 224 = 2·8n;
   dim of all cubics vanishing on {(pt,t₀,t₁)} = **3444**;
   **excess = 3444 − 3388 = 56 = 4n exactly, difference 0.** [VERIFIED n=14]
2. **n = 11, the D=4 per-S-box control:** falls(both) = 44, falls(S₀) = falls(S₁) = 22,
   **cross_excess = 0 → direct sum**, in a non-saturated regime (dim_low/C(33,≤3) = 0.35), with both
   exponents in the generic class. [VERIFIED n=11]
3. **n = 9** gives cross_excess = 1 — a single element, not a growing term; both exponents there are
   degenerate ({−2,(−2)⁻¹} and {2,2⁻¹}), which is why the per-S-box counts are 5n and 7n.
4. **n = 16 reports dim = 1 for the cubics on (t₀,t₁), and it is not a cross relation.**
   `explain_n16.py` shows the relation space is supported on the **t₁ variables alone**, and that
   Tr(t₁^{−r₁}) — the image criterion of S_{e₁} — vanishes on every graph point with Boolean degree
   wt(−r₁) = n − (e⁻¹ mod n) = 16 − 13 = **3**. At n = 14 the same criteria have degree 9 and 11,
   which is why the dimension there is 0. Sharpened claim: **the cubic relation space on (t₀,t₁)
   consists exactly of the univariate image criteria; its dimension is #{i : wt(−rᵢ) ≤ 3}; there is
   no cross-S-box cubic.** [VERIFIED n=14, 16]

### 3.4 E06 — kernel bilinear model and its random control

**Run:** n = 5, 7 at D ≤ 4 and n = 9 at D ≤ 3, with the new `random_bilinear_system(n)` control
(6n equations, blocks U = 3n / W = 3n, random U×W bilinear + affine, planted root) and a random
quadratic control of the S_quad shape. **File:** `results/E06_kernel_model_falls.csv`.

| | n=5 | n=7 | n=9 |
|---|---:|---:|---:|
| kernel, D=3, dim(span ∩ deg≤2) | 59 | **84** | **108** |
| 2·#eqs = 12n | 60 | 84 | 108 |
| kernel, D=3, syz | 1 | 0 | 0 |
| kernel, D=4, new falls | 937 | 1924 | — |
| kernel, D=4, saturation | 0.41 | 0.30 | — |
| **random bilinear, D=3 falls** | **0** | **0** | **0** |
| **random bilinear, D=4 falls** | **0** | **0** | — |

1. *Is the D=3 fall count always 2·#eqs = 12n?* **Yes at n = 7 and n = 9** (84 and 108 exactly).
   At n = 5 it is 59 = 12n − 1, together with a single syzygy — so the identity is
   [VERIFIED n=7,9] and **fails by one at n=5**.
2. *Does D=4 produce new falls?* **Yes**, a lot (937 and 1924), and in a **non-saturated** regime,
   unlike every full-system D=4 run.
3. *What does the random bilinear control do?* **Nothing**: at D=3 its row space contains only the
   6n equations themselves, and at D=4 `dim_low` equals rank_{D=3} exactly — zero falls at both
   degrees, for n = 5, 7, 9. So the kernel model's falls are structural, not a generic consequence
   of the bipartite bilinear shape. This is the single most useful control in the batch.

### 3.5 E05 — the two-valued-inverse reduced system, explicit

**Run:** `reduced_system_toy.py` (new), n = 8..14, seeds 1 and 2.
**File:** `results/E05_reduced_system_stats.csv`.

All four construction steps verified for every (n, seed): H_e as a linearised polynomial matches
`SboxStructures.H` on 100 random trace-zero points; {G_{e,0}(t), G_{e,1}(t)} = {x, x+a} on 100
random x ≠ 1; the true key is a root of both equations for the correct branch.

| n | Eq1 terms | Eq2 terms | deg_{t₀}(Eq1) | deg_{t₁}(Eq1) | Bézout deg·deg | 2^{2n} | roots (true + spurious) |
|---:|---:|---:|---:|---:|---:|---:|---|
| 8 | 64 | 678 | 254 | 254 | 64 516 | 65 536 | 2 = 1+1 |
| 10 | 100 | 1476 | 1022 | 1022 | 1 044 484 | 1 048 576 | 1 = 1+0 |
| 12 | 145 | 2689 | 4094 | 4094 | 16 760 836 | 16 777 216 | 3 = 1+2 |
| 14 | 197 | 4446 | 16382 | 16382 | 268 369 924 | 268 435 456 | 1 = 1+0 |

Eq1 has ≈ n² terms (n+1 at each Boolean degree 0..n−1), Eq2 ≈ 1.3–1.6 n³. The Bézout bound is
≈ 2^{2n} — exactly the input the T3.2 negative result needs. Root counts over both seeds and all n:
1 to 3, i.e. 0 to 2 spurious solutions, so pseudo-solution filtering is cheap.

**Two corrections to the ROADMAP §3.3 formulation (see §5 below).**

### 3.6 E08 — parameter audit

`exponent_audit.py` for n = 128, 192, 256 → `results/E08_audit_{128,192,256}.csv`, now carrying
`falls_D4_pred` and `falls_over_n` from the E03-calibrated rule. Degenerate exponents are rare at
the real sizes — 4 of 62 at n=128, 2 of 62 at n=192, 4 of 126 at n=256 — and **all ten recommended
exponents are generic**, asserted by the script itself. Note that at even n the two most damaging
classes ({2,2⁻¹} → 7n and {−2,(−2)⁻¹} → 4n) are unreachable because gcd(2,n) ≠ 1, so the
recommended field sizes are structurally protected from them.

### 3.7 E00 — KAT pinning: **done, not blocked**

`https://aimer-signature.org` was reachable and links the reference implementation at
`github.com/samsungsds-opensource/AIMer`; the repository holds `AIMer.zip` (commit `734cc76`,
*"update AIMer V3.0 reference implementation"*) with both `Reference_Implementation/` and `KAT/`.
`aim3.c`/`aim3.h` carry the AIM3 exponents, so it is the v3.0 primitive.

Reading the C fixed one convention axis `Conventions` did not have — the reference stores each
squeezed field element as a **column** and forms **A = L·U** — so `Conventions` now has four axes
(`tape_bit_order`, `matrix_apply`, `tape_matrix_order`, `lu_order`) and `ALL_CONVENTIONS`
enumerates the 16 combinations. `kat_pin.py` searched all 16 against the first three KAT records of
every parameter set:

| set | n | matching conventions |
|---|---:|---|
| AIM3-I | 128 | `lsb / A / col / LU` and `lsb / AT / row / LU` |
| AIM3-III | 192 | same |
| AIM3-V | 256 | same |

The two survivors are the same matrix (transposing both the tape reading and the application
cancels), so the convention is **unique** and is now the default in `src/aim3/aim3_ref.py`.
`docs/FACT_SHEET.md` §0.2 updated. The first three records of each set are committed to
`experiments/E00_reference_and_facts/kat_vectors.json`, so `kat_pin.py --vectors` reproduces the
result offline and `tests/test_core.py::test_kat_conventions_pinned` guards it.
**Consequence: full-size iv-dependent numbers may now be quoted** (CLAUDE.md ground rule 6 is
discharged).

---

## 4. **Gate 1**

> *Gate 1 (ROADMAP §7): if E04 finds cross-S-box degree falls → T2 becomes top priority; otherwise
> T2 and T3 run in parallel.*

**Verdict: no cross-S-box degree falls at D=4 — the falls are a direct sum, on the evidence
available at Tier A. T2 and T3 proceed in parallel; T2 is not promoted.**

The evidence, in decreasing order of strength:

| # | measurement | result | regime |
|---|---|---|---|
| 1 | E12, n=14: dim I(V)_{≤3}(pt,t₀,t₁) − rank Mac₃ (first layer) | **56 = 4n, difference 0** | ideal-theoretic, 16382 points vs 12384 monomials |
| 2 | E12, n=11: `cross_excess` at D=4 with the per-S-box control | **0** | non-saturated (0.35), generic exponents |
| 3 | E04, n=7,9,11, ℓ=1,2: syz₃(S_quad) | **(ℓ+1)·8n exactly** | direct sum at D=3 |
| 4 | E12, n=14 and 16: cubic relations on (t₀,t₁) | **only univariate image criteria** | no cross-S-box cubic exists to feed a cross fall |
| 5 | E04 full system at D=4 | **indeterminate** | saturated, so it neither supports nor contradicts |

The one non-zero residual is E12 at n=9 (`cross_excess = 1`), where both exponents are degenerate;
a single element that does not grow with n is not a cross-fall mechanism.

**Caveats, stated plainly.** This is a Tier-A verdict: (a) it covers ℓ = 1 (two S-boxes) and the
*first layer* of ℓ = 2, not the full three-S-box system, whose D=4 behaviour is saturated at every
size reachable in pure Python; (b) the strongest evidence is at a single n for each observable
(n=14 ideal-theoretic, n=11 per-block), which is below the ≥ 4 values of n that
`COMPLEXITY_ACCOUNTING §4` requires before extrapolating; (c) it says nothing about D ≥ 5.
**The server run in §7 is what would turn this into a [VERIFIED] statement across n.** Until then
the honest tag is **[VERIFIED n=11 (D=4 falls), n=14 (cubic ideal); CONJ in general]**.

**Consequence for priorities.** T1 keeps its closed-form target — the direct-sum structure is what
lets the low-degree Hilbert function of S_quad be written down. T2 gains a *different* boost from
this batch: the kernel model shows large, non-saturated D=4 falls exactly where the random bilinear
control shows **none** (§3.4), so the promising anomaly is in the kernel model's variable choice,
not in cross-S-box coupling. T3 is now concrete and cheaper than assumed (2 branch bits, §5).

---

## 5. Findings that contradict, or correct, what was expected

1. **E03 H3 was wrong in two places.** e = 2 gives **7n**, not 9n. "2e ≡ ±1 → 4n" conflates two
   different classes: 2e ≡ −1 gives 4n, 2e ≡ +1 gives 7n.
2. **`degenerate_class(kmax=4)` was wrong in both directions** — 11 false positives and 3 false
   negatives in 36 runs. Replaced (details in §3.1); this is a **decision for review** (§8).
3. **The D=4 fall count of the full system is not a usable observable at Tier A.** Every run
   saturates, *including the random control*. Any earlier reading of a full-system `falls_D4` at
   n ≤ 11 as structural should be discarded.
4. **Comparing a joint fall count with (ℓ+1)·2n is not a direct-sum test at all**, saturation
   aside, because of the cross-block multipliers. The batch's stated criterion is replaced by the
   per-S-box control; both numbers are recorded so the old criterion can still be read off.
5. **`cross_relations.py 16 3` returns 1, not 0.** It is the univariate image criterion
   Tr(t₁^{−r₁}) of Boolean degree wt(−r₁) = 3, not a cross relation (§3.3.4). Claim K survives and
   is sharpened; but a future run at deg = 4 must subtract these criteria before claiming anything.
6. **The ROADMAP §3.3 reduced system cannot be written as a sparse polynomial.** Its second
   equation composes a non-linearised sparse polynomial (G_{e₂,ε₂}, exponents ≈ 2ⁿ) with another
   one; the expansion has ≈ 2ⁿ terms. The equivalent mixed form — invert the *first* S-box, run the
   *last* one forward through its own P1 relation — is sparse, and this is what E05 builds.
7. **The reduced system has 2 branch bits, not 3.** The forward form of the last S-box has no
   branch ambiguity, so the attacker tries 4 branches, not 8. A factor 2 in T3's favour.
8. **E06's "2·#eqs falls at D=3" fails at n = 5** (59 instead of 60, with one syzygy). It holds
   exactly at n = 7 and n = 9.
9. **E03's `outside_ideal` is not symmetric under e ↔ e⁻¹**, although `falls_D4` is: e = 2 leaves
   2n cubics of the ideal unreachable at D=4, e = 2⁻¹ leaves none. **[OPEN]** — no explanation yet.
10. `scripts/smoke_test.sh` had never printed `smoke test OK` (pipefail/SIGPIPE); fixed.

## 6. Not done, and why

* **E05 resultant / Kronecker timings, polynomial decomposition** — explicitly out of scope for
  Tier A (Magma); the input format for the server is written up in the E05 README.
* **E04 at n = 13, 15 and ℓ = 3** — the first non-saturated full-system sizes, but out of reach in
  pure Python (n=13 at D=4 is ~9·10⁴ rows × 9·10⁴ columns). Queued in §7.
* **E12 deg = 4 on (t₀,t₁) at n = 16** — needs a 65k × 41k rank over F₂ (M4RI). Queued.
* **E06 at n = 11..16** — queued.
* **A uniform-GL_n matrix distribution as a second control** (E08 next-step (b)) — not run; every
  toy instance in this batch uses the spec's unit-LU.
* **The bit cost of one AIM3 evaluation** (COMPLEXITY_ACCOUNTING §1, spec Table 5 addition chain) —
  still the 2^14 placeholder.
* **No attack** is claimed anywhere in this batch. Everything here is structure and modelling.

## 7. Suggested first server batch (exact commands)

Run from the repository root on the server, in `tmux`, one Tier-C job above 50 GB at a time.
Validate every Magma/Sage skeleton against the Python CSVs at n ≤ 13 **before** trusting larger n
(`results/E03_sbox_degree_fall.csv` is the reference table).

```bash
# (1) Gate 1 at the first non-saturated full-system sizes — the highest-value job.
#     n=13 is the smallest n where the D=4 Macaulay matrix is not saturated (rows/cols ~ 0.99).
for n in 13 15; do
  nohup python3 experiments/E04_full_system_hilbert/full_system_falls.py $n -l 2 \
      --seeds 1,2,3 --systems S_quad,P1only,random --per-block --budget 86400 \
      > results/server/$(date +%Y%m%d)_E04/n${n}_l2.log 2>&1 &
done
# and the same for l=1 and l=3 (l=3 needs a 4-exponent instance = the AIM3-V shape).

# (2) E03 outside n<=16, to test the calibrated fall-count rule where it was NOT fitted.
bash experiments/E03_syzygy_and_degree_fall/sweep.sh 17 18 19 20
python3 experiments/E03_syzygy_and_degree_fall/summarise.py --md   # must still report 0 mismatches

# (3) E12: the first-layer system at D=4 for n=13,14, and deg=4 cross relations at n=16.
python3 experiments/E12_cross_sbox_relations/first_layer_ideal.py 13 --D4 --budget 86400
python3 experiments/E12_cross_sbox_relations/first_layer_ideal.py 14 --D4 --budget 172800
python3 experiments/E12_cross_sbox_relations/cross_relations.py 16 4     # M4RI: 65534 x 41449
python3 experiments/E12_cross_sbox_relations/explain_n16.py 16           # subtract the image criteria

# (4) E06: the kernel model's fall curve where it is still non-saturated.
for n in 11 13 15; do
  python3 experiments/E06_kernel_bilinear_model/kernel_bilinear_toy.py $n 4 --budget 86400
done

# (5) E05 (Tier C, Magma): resultants of the reduced system. Regenerate the instance with the same
#     (n, seed) as the Python run, emit both equations, and eliminate t1.
#     Input format: experiments/E05_two_valued_inverse/README.md, section "Input format for the server".
#     Fit log2(time) vs n for n = 8..24; the T3 termination criterion is slope >= 2.

# (6) E04 with SetVerbose("Faugere",2) / msolve -v 2 to log F4 step degrees — the observable that
#     saturation does not destroy, and the one that actually answers Q2 (d_solv vs n).
```

## 8. Open decisions for a human

1. **Keep the recalibrated `degenerate_class`?** It changes the public meaning of the tags (they are
   now the classes {2,2⁻¹}, {−2,(−2)⁻¹}, {1+n/2}, {−3,(−3)⁻¹}, {n/2−1} instead of `ke ≡ ±1,±2`) and
   the `degenerate` column of `results/E08_audit_*.csv`. The old rule is still reachable via
   `legacy_kmax=4`. Recommended: keep it — it is the only version that matches the measurements.
2. **Is the Gate-1 verdict strong enough to act on?** It rests on one n per observable, below the
   extrapolation discipline of COMPLEXITY_ACCOUNTING §4. The alternative is to hold Gate 1 open
   until the server job §7(1) reports n = 13 and 15.
3. **Should the paper state the fall-count rule as [CONJ]?** It is a fit to 36 runs at n ≤ 16 with
   no proof and no test outside that range; §7(2) is the cheapest way to falsify it.
4. **E05's mixed-form reduced system** replaces the ROADMAP §3.3 formulation. The ROADMAP text
   should be updated (it currently describes a system that cannot be expanded); this batch only
   appended a status note, it did not rewrite §3.3.
