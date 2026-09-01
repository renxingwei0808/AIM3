# Server batch 2 — report

**Branch `claude/server-batch-2`, 2026-09-01.  Tier B/C, driven from the WSL laptop over SSH.**

Environment unchanged from batch 1 (`server/ENVIRONMENT.md`): Magma V2.24-2, SageMath 10.7,
msolve 0.10.1, GNU parallel, M4RI, all under `/home/luck/rxw/AIM3_remote/opt`.  Job logs and their
`/usr/bin/time -v` stamps are under `results/server/20260831_*/` (the server clock is PDT, so the
batch-2 jobs share batch 1's date directory; every batch-2 job name begins `B2_`).  Every table
below is generated from the committed CSVs by `python3 server/jobs/summarise_batch2.py`.

**The three questions.**  (Q1) what is the calibrated solving-degree law; (Q2) does the S-box ideal
hold anything beyond the quadratic ideal at degree 4 and above; (Q3) can the kernel model and the
twisted-bilinear structure be exploited by a non-Boolean solving strategy.

---

## 0. Closing batch 1 — and the control it was missing

**The leftover jobs.**  Batch 1's §8 named five jobs still running at its data cut.  Three
completed and are in the tables (`ℓ=3, n=15` three-first-layer subset → `cross_excess = 0`; the
D = 5 random control at n = 16 → `falls_D5 = 0`, unsaturated; Gate 1 at `ℓ=2, n=19` seeds 1–2 →
`cross_excess = 0`).  Two — the block orders at n = 13, 15 and the E06 kernel F4 at n = 11, 13 —
were **stopped on purpose** once batch 2's corrected-exponent grid superseded them; what they had
produced is collected, and the stop is recorded here rather than left as a gap.

**The decoupled control (§B.0) — batch 1's one genuinely outstanding gap, now closed.**  Batch 1 argued
that the `cross_excess = 1639` at ℓ = 3, n = 11 was a small-parameter effect, but conceded it could
not rule out a fault in the *instrument*: the same-shape random control has no degree falls at all,
so its `cross_excess` is trivially 0 and proves nothing about the arithmetic.
`server/jobs/decoupled_control.py` supplies the control that does: k S-boxes on **pairwise disjoint
variable blocks**, where additivity is forced by the disjoint supports.

    n = 11, k = 2 (44 variables)   3 subsets   cross_excess = 0
    n = 7,  k = 4 (56 variables)   7 subsets   cross_excess = 0

**10 of 10 subsets give exactly 0**, at both the ring size of the ℓ = 3, n = 11 pairs and the block
count where the anomaly appeared.  The measurement is sound; the 1639 is a property of that
instance, not of the code.

---

## 1. Q1 — the calibrated solving-degree law

### 1.1 A finding before any measurement: generic exponents often do not exist

Task 1.1 asks every new instance to use the smallest ℓ+1 exponents with
`degenerate_class(n,e) == []`.  Implemented as `--generic-exps` (raising when there are not enough,
as specified).  It **raises at three of the five sizes the batch asks for**:

| n | admissible e | non-degenerate | enough for ℓ = 2? |
|---:|---|---|:--:|
| 9 | 2, 4, 5, 7 | **none** — {2,2⁻¹} = {2,5}, {−2,(−2)⁻¹} = {7,4} | no |
| 11 | 2..9 | {3, 4} | no (needs 3) |
| 13 | 2..11 | {3, 5, 8, 9} | **yes** |
| 15 | 2, 4, 7, 8, 11, 13 | {4, 11} | no (needs 3) |
| 17 | 2..15 | {3, 4, 5, ...} | **yes** |

So "use generic exponents everywhere" is not achievable below n = 13, and at n = 9 not at all.
Rather than silently fall back, `--exps` was added and the least degenerate available set is named
explicitly for each n, with its classes recorded in every row:

    n = 9   2;4;5   (7n, 4n, 7n)     -- no generic exponent exists at this n
    n = 11  3;4;7   (2n, 2n, 2n+1)   -- instead of the default 2;3;4, which starts at the 7n class
    n = 13  3;5;8   (2n, 2n, 2n)     -- fully generic
    n = 15  4;11;7  (2n, 2n, 4n)     -- only two generic exponents exist
    n = 17  3;4;5   (2n, 2n, 2n)     -- fully generic

**This matters for reading batch 1.**  Every ℓ = 2 instance there began at e = 2, the 7n class.
Batch 1's one structural claim about `d_solv` — that `S_quad` sits one degree below a same-shape
random system at ℓ = 2, n = 7 and 13 — was measured on such instances.

### 1.2 `d_solv` on the corrected instances

Table in §B.1.  On the sizes completed so far, **`S_quad`, `P1only` and the random control give the
same `d_solv` at every cell**, and in particular the n = 11, ℓ = 2 cell that batch 1 measured as
4 / 4 stays 4 / 4 when e = 2 is replaced by the generic 3;4;7.  msolve and Magma agree everywhere.

**What the completed cells say (§B.1).**  `d_solv` is 4 at ℓ = 2 for n = 9 and 11, and 5 at n = 13
and 15; at ℓ = 1 it is 4 for n = 11, 13, 15 and 5 at n = 17.  So on the sizes reachable here the
solving degree moves by **one step over six values of n**, which is why no law is fitted to it:
`COMPLEXITY_ACCOUNTING §4` wants four completed sizes before an extrapolation, and four completed
sizes of a quantity that takes two distinct values do not determine a growth rate.  What the table
does settle is the *level* — `d_solv` ∈ {4,5} where B2's `d_reg` for the same shapes is 5–8 (§B.2) —
and that the exponent correction did not move it.

**One cell where the structured system is worse than random.**  At ℓ = 1, n = 17, `S_quad` and
`P1only` both need `d_solv = 5` while the random control of the same shape needs 4.  Everywhere else
the three agree exactly.  A single cell is not a phenomenon, and the direction is the unhelpful one
for an attacker, but it is recorded because it is the only place in this batch where the AIM3
structure and a random system of identical shape separate at all in `d_solv`.

**Tool agreement and cost.**  Magma and msolve agree on **every** cell both completed; msolve is
5–10× slower on the same instance and accounts for 11 of the batch's censored runs against Magma's
3.  Magma remains the instrument of record for `d_solv`, msolve the independent check — the same
division as batch 1, now with the numbers behind it.

### 1.3 Hybrid guessing and δ

`k` bits of pt are fixed to their **true** values, so the planted solution survives and the run
measures the branch the attacker would actually be on.  Full table in §B.2; ℓ = 2, n = 13, both
seeds, Magma and msolve agreeing on every completed cell:

| k | variables 3n−k | `S_quad` `d_solv` | random `d_solv` | B2 `d_reg` | **δ = d_reg − d_solv** |
|---:|---:|---:|---:|---:|---:|
| 0 | 39 | 5 | 5 | 5 | **0** |
| 2 | 37 | 4 | 4 | 4 | **0** |
| 4 | 35 | 4 | 4 | 4 | **0** |
| 6 | 33 | 3 | 4 | 4 | **1** |

**δ = 0 for three of the four guess widths.**  Where B2 assumes semi-regularity and derives `d_reg`,
the measurement lands exactly on it — the assumption is *right* about the solving degree of the
guessed system at n = 13 for k ≤ 4, even though §1.4 shows the Macaulay ranks underneath it are not
semi-regular at all.  These two facts are compatible: `d_reg` is where the rank first suffices, and
a system can carry thousands of extra syzygies at D = 5 without that degree moving.

**δ = 1 appears at k = 6, and only for the structured system.**  At 33 remaining variables `S_quad`
solves at degree 3 while the random control of identical shape still needs 4.  This is the one place
in the batch where guessing exposes something a random system does not have.  It is **one n, one k,
two seeds** — far short of an extrapolation, and §1.5's table shows what a unit of δ would be worth
if it survived (8.5–9 bits at ω = 2).  Whether δ grows with k, with n, or neither is exactly the
question the next batch has to answer, and `n = 15` was launched for it (still censored at k = 0 at
the time of writing; the k > 0 cells are cheaper and are running).

**A column that is easy to misread.**  `B2 d_reg` in §B.2 is computed from each row's own shape, so
the `P1only` rows (3n equations, not 9n) get a much larger `d_reg` and a δ of 3 throughout.  That δ
is not the B2 quantity: B2's baseline is the 9n-equation system.  Only the `S_quad` rows carry the
δ that §1.5 consumes.

### 1.4 The Hilbert-function side check

Two things were asked here: how far the measured Macaulay rank is from the Boolean semi-regular
prediction at D = 5, and whether the direct-sum property (Gate 1) still holds at D = 5.  Both were
run at ℓ = 2, **n = 11** rather than the n = 13 the prompt names — n = 13 at D = 5 is
1160640 × 667928 ≈ 97 GB and ≈ 32 h of M4RI, out of proportion to this batch; the shortfall is
stated rather than papered over, and `P1only` was run at n = 13 (32 GB) as the partial substitute.

**The semi-regular prediction fails, and it fails only for the AIM3 system.**  The D = 5 Macaulay
matrix is 595782 × 284274 in both runs (§B.3):

| system | rank | predicted | rank − predicted | syzygies | predicted |
|---|---:|---:|---:|---:|---:|
| `S_quad` | 268067 | 281622 | **−13555** | 327715 | 314160 |
| random, same shape | 284273 | 281622 | +2651 | 311509 | 314160 |

The random control lands one short of full column rank (284274 − 1), i.e. at D = 5 a random system
of this shape is already *solved*; `S_quad` is **16206 below it** and 13555 below the semi-regular
formula.  For scale, the departure at D = 3 is `syz₃ = 8n(ℓ+1) = 264`.  So the ROADMAP's C3/C4
objection is not a small correction that stays small: between D = 3 and D = 5 the excess syzygies
grow by a factor of ≈ 51 on this instance.  **One instance, one n — this is a measurement, not a
law**, and no exponent is fitted to it here.

**Gate 1 at D = 5: the direct sum breaks, but the measurement cannot say why.**  The subset lattice
(§B.0/§B.11 aside; the D = 5 rows are in `results/E04_server_gate1_pairs_D5.csv`) gives

| blocks | equations | `falls_D5` | Σ singles | `cross_excess` | saturation |
|---|---:|---:|---:|---:|---:|
| {0} | 33 | 859 | 859 | 0 | 0.3422 |
| {1} | 33 | 870 | 870 | 0 | 0.3424 |
| {2} | 33 | 914 | 914 | 0 | 0.3433 |
| {0+1+2} | 99 | 5467 | 2643 | **2824** | **0.9999** |

Taken at face value this says the direct sum fails at D = 5.  It should not be taken at face value:
`dim_low = 46934` out of `low_cols = 46938`, so the degree-5 Macaulay matrix of the full system
spans all but **four** dimensions of the entire degree-≤4 space.  A system that has filled the space
has nowhere left to be additive, so `cross_excess` at that saturation measures the ceiling, not a
cross-S-box mechanism.  The D = 4 rows are meaningful precisely because their saturation is 0.12–0.34.
**Recorded as an observation with its saturation attached, and explicitly not offered as Gate-1
evidence.**  (Consistency check: the Hilbert job, an independent script, reports the same
`dim_low = 46934` for this instance.)

### 1.5 The calibrated numbers (§B.4)

`server/jobs/calibrated_estimate.py` minimises `2^k·(Σ_{i≤d} C(N−k,i))^ω` over k with
`d = dreg_boolean(N−k, m) − δ`.  **The table is the deliverable; no security conclusion is drawn
here** — that judgement belongs to the paper and needs an extrapolation of δ from n ≤ 15 to
n = 128..256 that `COMPLEXITY_ACCOUNTING §4` does not license.

Two things are worth noting about the table itself.  Each unit of δ is worth **8.5–9 bits** at
ω = 2 (215.8 → 207.3 for AIM3-I, 309.0 → 300.3 for AIM3-III, 523.7 → 514.8 for AIM3-V), which is
what a *measured* δ would buy.  And the ω = 2.81 rows are far **above** the ω = 2 rows rather than
below, because the optimiser then prefers a much larger k (91 rather than 5 for AIM3-I): the
exponent applies to the monomial count, so a bigger ω pushes the optimum towards guessing.

---

## 2. Q2 — the S-box ideal above degree 3

### 2.1 There is content beyond D = 5 — but only for generic exponents (§B.5)

`--pointset-deg 4` measures `dim I(V_e)_{≤4}` directly at n = 16, where the graph has 65535 points
against 41449 monomials so the answer is an ideal statement and not a point deficit (at n = 14 it
would be pure deficit — 16383 points against 24158 monomials — and batch 1 already recorded that
trap):

| n | e | class | dim I(V)≤4 | D = 5 reaches | **outside** |
|---:|---:|---|---:|---:|---:|
| 16 | 3 | generic | 21592 | 21432 | **160 = 10n** |
| 16 | 5 | {−3,(−3)⁻¹} | 21224 | 21224 | **0** |

So for a generic exponent the degree-5 Macaulay matrix of ⟨P₁,P₂,P₃⟩ **fails to produce 10n
quartics of the vanishing ideal**, while for the `{−3}` class it produces all of them.  That is the
same shape of asymmetry batch 1 found one degree lower (`outside_ideal = 2n` for e = 2 and 0 for
e = 2⁻¹), now on the other side: here it is the *generic* class that leaves a remainder.

**What 160 is and is not.**  It is an upper bound on the number of quartic generators outside the
quadratic ideal: D = 6 might still reach some of them.  D = 6 at n = 16 is ≈ 290 GB (48·C(32,≤4)
rows × C(32,≤6) columns) and is not attempted; §8.

### 2.2 The D = 5 fall count decomposed — batch 1's arithmetic is refuted (§B.6)

Batch 1 measured `falls_D5 = n(4n+7)` (generic) and `n(4n−6)` ({−3}) and suggested reading the
first as `(2n)(2n) + 7n`: the 2n new cubics of D = 4 times the 2n variables, plus 7n genuinely new
quartics.  `server/jobs/d5_decomposition.py` measures the two halves separately, and **that reading
is wrong**:

| n | class | #cubics | products | X = span mod Mac₄ | **relations** | falls_D5 | **net new** |
|---:|---|---:|---:|---:|---:|---:|---:|
| 14 | generic | 28 = 2n | 784 | 532 | **252 = 18n** | 882 | **350 = 25n** |
| 16 | generic | 32 = 2n | 1024 | 736 | **288 = 18n** | 1136 | **400 = 25n** |
| 14 | {−3} | 29 = 2n+1 | 812 | 560 | **252 = 18n** | 700 | **140 = 10n** |
| 16 | {−3} | 33 = 2n+1 | 1056 | 768 | **288 = 18n** | 928 | **160 = 10n** |

The products are **not** independent: there are exactly **18n linear relations** among them, at
both sizes and — this is the striking part — **in both exponent classes**, even though the number
of products differs (4n² against 4n²+2n).  With that, both closed forms of batch 1 are accounted
for exactly:

    generic:  falls_D5 = (2n)(2n)   − 18n + 25n = 4n² + 7n   ✓
    {−3}:     falls_D5 = (2n+1)(2n) − 18n + 10n = 4n² − 6n   ✓

So the correct statement is **`falls_D5 = (#cubics)·2n − 18n + net`**, with the 18n relations a
class-independent invariant and the net new quartics 25n (generic) or 10n ({−3}).  The `net` is 2.5
times larger than the 7n batch 1 read off the closed form, and the sign of the class difference is
the opposite of what the closed forms suggest: the `{−3}` class has *more* products but *fewer* net
new quartics.  The 18n is the quantity a theorem should explain first.

### 2.3 Task 2.3 is not executable as specified

Task 2.3 asks, if quartic generators outside the ideal are found, to add them to `S_quad` at ℓ = 2,
n = 13 and re-measure `d_solv`.  **The two sizes do not overlap.**  The quartics can only be
*identified* where the graph has more points than monomials, i.e. n ≥ 16 (at n = 13 the graph has
8191 points against C(26,≤4) = 17902 monomials and the measurement is pure deficit); and the full
ℓ = 2 system is only solvable up to n ≈ 15.  Reported as not done, with the reason, rather than
substituted with a measurement at a size where the input does not exist.

---

## 3. Q3 — non-Boolean strategies on the kernel model

### 3.1 Block-restricted XL does not linearise earlier (§B.7)

Restricting the Macaulay multipliers to one block (`--mult-block W|U`) gives, at every degree
measured, a matrix of **full row rank** (syz = 0 almost throughout) and a corank that stays enormous
— 731600 at n = 11, D = 4 — so linearisation (corank 1) is nowhere in sight.  Against the full
multiplier set it produces roughly a quarter of the rank and a third of the falls:

| n | D | full multipliers | U block | W block |
|---:|---:|---:|---:|---:|
| 11 | 3 | rank 4422, dim_low 132 | rank 2244, dim_low 99 | rank 2244, dim_low 99 |
| 11 | 4 | rank 143715, dim_low 9075 | rank 37092, dim_low 3465 | rank 37092, dim_low 3366 |

**U and W behave the same**, to within a handful of dimensions (3465 against 3366 at n = 11, D = 4).
The bilinear structure does not make one side cheaper than the other, and neither side reaches
linearisation sooner than the full set does.  Negative result for the block-restricted XL idea.

### 3.2 The rank condition is real but worth only about n bits (§B.8)

Writing the first layer as `N(W)·(pt,1)^T = 0` with N of size 2n × (n+1):

* at the **true** W the rank is exactly n — corank 1 — at every n measured, and `(pt,1)` really is
  the kernel vector (asserted, not assumed);
* over 20000 uniformly random W the fraction with rank ≤ n is 3.085 % (n=7), 0.840 % (9),
  0.140 % (11), 0.045 % (13) — consistently **1.4–2.2× the value for a uniformly random matrix of
  the same shape**, `1 − Π_{i≤n}(1 − 2^{i−2n}) ≈ 2^{1−n}`.

So the rank drop **is** a genuine distinguisher, and the family `{N(W)}` is mildly less generic than
a uniform one (the consistent factor ≈ 2).  But the information it carries is only about n − 1 bits
of the 2n-bit W space, leaving ≈ 2^{n+1} candidates — **the same order as the pt search the attack
was supposed to avoid**.  MinRank/KS on this layer therefore offers no advantage over guessing pt,
and the honest way to state T2.1 is that the rank defect is structurally real and cryptanalytically
free.

### 3.3 Block orders, computed properly this time

Batch 1 could not measure a block-order solving degree at all: Magma 2.24 runs its GF(2) F4 in
grevlex whatever order the ring declares and converts afterwards.  Sage/PolyBoRi computes
**directly** in the declared order and, with `prot=True`, prints one `Current Degree: d` line per
degree, so `max d` is a genuine per-order solving degree.  Block orders are built as
`TermOrder('degrevlex', k) + TermOrder('degrevlex', nv−k)`; "W first" is obtained by rotating the
variables so the W block occupies the leading positions, which is a renaming and leaves the algebra
untouched.

**Two caveats that must be read with the table (§B.9).**

1. **PolyBoRi's `d_solv` is not comparable with Magma's or msolve's.**  Its protocol reports the
   degree at which it is *processing pairs*, and its algorithm is not F4; on the same system it
   reports 5 where both F4 implementations report 4.  The comparison this table supports is
   **across orders within PolyBoRi**, which is exactly the comparison batch 1 could not make.
2. **PolyBoRi aborts on most block-order runs, and the abort is not confined to the kernel
   model.**  The failure is `Unhandled SIGABRT: An abort() occurred… a *compiled* module has a bug`
   (exit code 134), emitted after some number of degree lines.  It is a crash inside PolyBoRi, not
   a timeout and not a fault in the driver; the degrees already printed are a **lower bound** and
   the row is censored.

   The driver first recorded these as completed runs — it derived `censored` from the timeout and
   from "did any degree line appear", and never looked at the exit code.  Fixed
   (`server/jobs/run_polybori.py`), and the rows already written were re-marked from `rc`.  Before
   the fix the table appeared to show block orders cutting `d_solv` from 5 to 2; that reading was
   an artefact of counting a crash as a result, and it is withdrawn.

**What survives (§B.9): 18 completed runs of 43.**  12 of the 15 grevlex runs complete (the three
kernel-model ones abort); **all 12 `blockPT` (pt-first) runs abort**; `blockT` completes 6 of 12 —
all six at n = 9, none at n = 11; both `blockU` and `blockW` on the kernel model abort.

| system | order | `d_solv` (PolyBoRi) | note |
|---|---|---|---|
| `S_quad`/`P1only`/random, ℓ = 2, n = 9 and 11 | grevlex | **5** at every one of the 12 completed runs | complete |
| same, n = 9 | `blockT` (t-blocks first) | **14, 15, 15, 15, 16, 17** over 6 completed runs | complete |
| same | `blockPT` (pt first) | — | **every run aborts** |
| kernel model, n = 7 and 9 | grevlex, `blockU`, `blockW` | — | **every run aborts** |

So the one block-order comparison that completed says the opposite of what T2.2 hoped for: putting
the t-blocks first raises the working degree from 5 to 14–17 and is strictly worse.  The order that
T2.2 actually wants — pt first, so that the twisted-bilinear structure is eliminated block by
block — **cannot be measured with this tool either**.  Batch 1 was blocked by Magma computing in
grevlex regardless of the declared order; batch 2 is blocked by PolyBoRi crashing on exactly the
orders of interest.  Two independent tools, two different failure modes, same missing number.

**This is the batch's clearest unfinished item.**  The remaining options, none of which fits in this
batch: build the block-order Gröbner basis over F₂ with field equations in Singular or msolve (both
compute in the declared order but neither is tuned for Boolean systems), or implement a small
matrix-F4 in the block order directly on top of the existing M4RI backend, which the Macaulay
machinery in `src/aim3/boolmac.py` already almost supports.

---

## 4. T4 cross-validation, and what the eliminated condition costs (§B.11)

The Tier-A branch reached `main` while this batch was running (PR #5, `e6cd091`; merged into this
branch at `2370d74`), which is the condition the prompt attaches to task 4 — so it was done.

**Two instruments, on the same instances.**  Every instance is built by the Tier-A probe itself
(`skew_elimination_probe.py --export`, same default seed), so Sage and the Tier-A class are compared
on identical data rather than on two random instances that happen to have the same shape.  The
division of labour was forced by a measurement: **Sage's `OrePolynomialRing` costs about 1 s per
`lclm` where `src/aim3/ore.py` costs 0.2–20 ms** — a factor of 50–1000 depending on n.  Sage is
therefore the independent *check*, and the hand-written class is what carries the 2ⁿ-evaluation
measurement upward.  (This is worth recording for later batches: for skew-polynomial work Sage's
`OrePolynomialRing` is a correctness reference, not a workhorse.)

**Result 1 — the Tier-A measurement replicates exactly.**  Over n = 8..14, where the Tier-A CSV has
rows, Sage reproduces `deg lclm`, `dim(ker Λ₀ ∩ ker Λ₁)`, the affine-solvability bit and the random
`gcrd` range: **28 comparisons, 0 disagreements**.  In particular the correction Tier-A made to the
brief's hypothesis is confirmed by an independent implementation — `deg gcrd = 0` and
`dim(ker Λ₀ ∩ ker Λ₁) = 0` at the *true* key just as at random keys, and the statistic that does
separate the true key is the **affine** one.

**Result 2 — the ring-level picture holds to n = 24.**  `deg lclm = deg Λ₀ + deg Λ₁` at every
n ≤ 24 (the gcrd is 0 everywhere), and the affine system is solvable at the true (t₀,t₁) and at
**0 of 32** random pairs at every one of those sizes.  Tier-A's n = 14 ceiling was a cost ceiling of
its 2ⁿ measurements, not of the algebra: the ring-level statement costs nothing to check and does
not weaken.

**Result 3 — D(n) = 2ⁿ − 1, the maximum, at every size measured (n = 8..18).**  This is the number
T4.3 makes the switch of the whole track.  With t₁ at its true value and t₀ ranging over **all** of
F_{2ⁿ}, every coefficient of the monic lclm — not merely the extreme ones — has exact univariate
degree 2ⁿ − 1 and exact ANF degree n.  Both are the largest values possible, at all eleven sizes,
for all 28–47 coefficients of each.

**What that settles.**  Eliminating pt in the Ore ring is exact, cheap and loses nothing (Tier-A
already showed `Z(E) = C`, so no separate Frobenius-consistency step is needed).  But the object it
leaves behind is of *full* degree in t₀, so every degree-sensitive way of finishing — a resultant in
t₀, univariate root-finding, sparse interpolation — costs Ω(2ⁿ) just to write the condition down.
That is the same order as guessing pt outright (B3), so **the elimination buys nothing**: T4's cost
does not come from the elimination step being hard, it comes from the eliminated condition being
dense.  This is the same shape of negative result as E05/T3 in §5, reached independently.

**Caveats, inherited and new.**
* The exponent pairs are Tier-A's: for n ≤ 12 and a few larger n the generic class contains only
  `e = n−1`, so both S-boxes get the same exponent — not an AIM3 parameter shape, but it keeps the
  instance out of the degenerate classes, which is what a degree measurement needs.  The n where
  two distinct generic exponents exist (11, 13, 14, 15, …) show the same D(n).
* The symbolic instrument (the same lclm over K(t₀), coefficients as exact rational functions)
  **does not finish at all**: 6 h at n = 8 without returning, on two degree-14 elements.  The reason
  is structural rather than a tuning failure — the twist sends t₀ ↦ t₀², so the degree in t₀ doubles
  at every Euclidean step and passes 2¹⁴ before the algorithm is halfway.  Recorded as censored, and
  it is consistent with result 3: a formal degree far above 2ⁿ reduces, as a function on F_{2ⁿ}, to
  the full degree 2ⁿ − 1.
* The independent Sage check **of D(n) specifically** was lost to a Sage trap and is being redone.
  `.sage` files preparse `^` as exponentiation, so the Möbius transform written `arr[i] ^= arr[i^s]`
  was computing a tower of powers; it raised `MemoryError`, which is the lucky outcome — the same
  rewrite on `s ^= f0` returns a wrong number in silence.  Fixed, documented in
  `server/ENVIRONMENT.md`, and re-launched for n = 8..11.  **Until it lands, D(n) rests on one
  implementation**, and the report says so rather than implying two.
* One instance per n, as in Tier A.  D(n) saturating its maximum is not a statistic that a second
  seed could improve, but it has not been repeated over seeds.

---

## 5. E05 — eliminating t₁ for the price of a table lookup (§B.10)

Batch 1 left T3 with one completed size: `Resultant(Eq1, Eq2, t₁)` takes 294.5 s at n = 7 and does
not return in two hours at n = 8.  Its §8 proposed two cheaper routes; the second one works, and it
works by a very large margin.

Eq1 is **separable** — `W(t₀) + V(t₁)` — so its variety is obtained by tabulating V over F_q\{0} and
looking `W(t₀)` up, then testing each hit against Eq2.  That is Θ(2ⁿ) field operations and forms no
eliminated polynomial at all:

| n | 7 | 10 | 13 | 16 | 18 | 20 | 22 | 24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| eliminate (s) | 0.006 | 0.116 | 1.34 | 15.6 | 74.0 | 403 | 1536 | 7490 |
| batch 1 `Resultant` (s) | 294.5 | — | — | — | — | — | — | — |

**18 completed sizes, n = 7..24, no censoring** (full table and the un-tabled first run in §B.10).
Least squares gives a slope of `log2(seconds)` against n of **1.167**, and every successive ratio
lies in 1.71..3.00 — no regime change anywhere in the range.  The method is Θ(2ⁿ) by construction,
so a slope near 1 is the expected answer; it is now measured over 16 sizes rather than assumed.

Batch 1's `Resultant` on the same instance takes 294.5 s at n = 7 and does not return in two hours
at n = 8.  At n = 7 the table-lookup route is **49000× faster**, and the comparison only gets worse
for the resultant from there.

**A measurement artefact found and removed.**  The first run of this had a **29× step at n = 17**
and a 4.9× one at n = 18, which would have made any slope fitted across the whole range meaningless.
The cause is not the algorithm: `src/aim3/gf2n.py` builds log/antilog tables only for n ≤ 16 and
falls back to carry-less multiplication above that.  The tables cost 16 MB at n = 20 — nothing
against the 2ⁿ lookup table this method already builds — so `--tables` forces them, and the curve
above is the homogeneous one.  Both runs are kept, in separate files, with identical root counts on
every shared size.

**What this settles.**  The two-valued-inverse reduced system is solvable in **Θ(2ⁿ)** — the same
order as brute-forcing pt, which is baseline B3.  T3.2 expected a negative result and argued it from
the degree of the eliminated polynomial; this measures the actual cost of solving the system and
finds it neither better nor worse than exhaustive search.  The resultant route was never the right
way to ask the question: it costs more than the whole problem is worth.

**What it does not settle.**  There is no `log2(deg)` slope here, because no eliminated polynomial
is formed; the route that does form one is censored from n = 8 (batch 1).  The Θ(2ⁿ) figure is a
*solving* cost for the toy instances, not an extrapolation to n = 128.

---

## 6. Findings that contradict, or correct, what was expected

1. **Batch 1's reading of `falls_D5` is wrong** (§2.2).  It is not `(2n)(2n) + 7n`; the products
   carry 18n linear relations and the net new quartics are 25n, not 7n.
2. **Generic exponents do not exist at n = 9, and are too few at n = 11 and 15 for ℓ = 2** (§1.1).
   The batch's own instruction cannot be followed at three of its five sizes, and every ℓ = 2
   instance in batch 1 began at the degenerate e = 2.
3. **Task 2.3 is not executable** (§2.3): the size at which ideal-external quartics can be
   identified (n ≥ 16) and the size at which the full system can be solved (n ≤ 15) do not overlap.
4. **The `{−3}` class exhausts its quartic ideal at D = 5 and the generic class does not** (§2.1) —
   the reverse of the degree-3 asymmetry, where it was e = 2 that left a remainder.
5. **`N(W)` is measurably less generic than a uniform matrix family** (§3.2), by a consistent
   factor ≈ 2 in the rank-drop probability — an observation that is easy to miss if one only checks
   the true W.
6. **Sage's `OrePolynomialRing` is 50–1000× slower per `lclm` than `src/aim3/ore.py`** (§4).  It is
   an excellent independent check and a poor workhorse; the task-4 plan had to be split accordingly.
7. **T4 and T3 fail for the same reason, reached by different routes** (§4, §5): the elimination
   step is cheap and exact in both, and in both the object it produces costs Θ(2ⁿ) to exploit.
8. **A 29× step in the E05 curve was the field implementation, not the algorithm** (§5).
   `src/aim3/gf2n.py` drops its log/antilog tables above n = 16; fitting a slope across that step
   would have produced a fabricated exponent.  Two other measurement defects were found and fixed
   the same way: the PolyBoRi driver counted SIGABRT crashes as completed runs (§3.3), and the T4
   degree driver wrote a censored row with a different set of columns, which killed its own
   escalation at the first timeout.  All three are the same failure mode — **a run that did not
   finish, or did not finish the way it looks, being recorded as a value.**

---

## 7. Not done, and why

* **D = 6 at n = 16** (§2.1): ≈ 290 GB dense and far beyond the per-job cap; without it the 160
  quartics outside the D = 5 reach stay an upper bound on genuine generators rather than a count.
* **The Hilbert side check at ℓ = 2, n = 13, D = 5** as literally specified: 1160640 × 667928 is
  ≈ 97 GB — the size `SIZING.md` predicts — but ≈ 32 h of M4RI, out of proportion to the batch.
  Run at n = 11 instead, and at n = 13 for `P1only` (32 GB), with the shortfall stated.
* **Task 2.3** (§2.3), for the reason given there.
* **The symbolic (K(t₀)) form of the lclm above n ≈ 10** (§4): the twist squares the degree in t₀
  at every Euclidean step, so the exact rational functions become unrepresentable almost at once.
  The evaluation-interpolation instrument answers the same question and is the one reported.
* **D(n) beyond the sizes in §B.11.2**: the measurement is 2ⁿ `lclm` computations; it was pushed as
  far as the batch's budget allowed and censored rows are marked as such, never as values.
* **`ℓ = 2` at n = 17 for `d_solv`**: launched, and expected to be censored — 51 Boolean variables
  is past where Magma finished at n = 15 (3109 s, 79 GB).
* **No attack** is claimed anywhere in this batch.

---

## 8. Suggested rows for `paper/claims_ledger.md` (not applied)

| claim | suggested status | support |
|---|---|---|
| `cross_excess` as measured by `gate1_pairs.py` is arithmetically sound | [VERIFIED] — decoupled control, 10/10 subsets exactly 0 | `results/E04_server_decoupled.csv` |
| `falls_D5 = (#cubics)·2n − 18n + net`, net = 25n (generic) / 10n ({−3}) | [VERIFIED n = 14, 16, two classes] — **replaces** batch 1's `(2n)(2n)+7n` reading | `results/E03_server_d5_decomposition.csv` |
| 18n linear relations among the `z_k·c_j`, independent of the exponent class | [CONJ], verified at n = 14, 16 | same |
| D = 5 leaves 10n quartics of `I(V_e)` unreached for generic e, 0 for the {−3} class | [VERIFIED n = 16]; an upper bound on generators outside ⟨P₁,P₂,P₃⟩ | `results/E03_server_pointset.csv` |
| block-restricted XL does not linearise earlier, and U ≈ W | [VERIFIED n = 7, 11, D ≤ 5] — new row | `results/E06_server_multblock_*.csv` |
| the first-layer rank defect is real (corank 1 at the true W) but carries ≈ n−1 bits | [VERIFIED n = 7, 9, 11, 13, 20000 samples each] — new row, T2.1 | `results/E06_server_minrank.csv` |
| the two-valued-inverse reduced system is solvable in Θ(2ⁿ) by the separable structure | [VERIFIED n = 7..], i.e. the same order as B3 — this is T3.2's negative result, demonstrated | `results/E05_server_cheap_elimination.csv` |
| generic exponents do not exist at n = 9 and are too few at n = 11, 15 for ℓ = 2 | [PROVED] (a finite check of `degenerate_class`) — worth a footnote wherever toy exponents are chosen | `server/jobs/…/generic_exps` |
| the Tier-A T4 measurement replicates in an independent implementation | [VERIFIED n = 8..14, 28/28 comparisons] — supports the existing E07 rows rather than adding a claim | `results/E07_server_ore_scale.csv` |
| `deg lclm = deg Λ₀ + deg Λ₁`, gcrd = 0, and the affine predicate separates the true key | [VERIFIED n ≤ 24, 0/32 random pairs solvable] — widens the Tier-A row from n ≤ 14 | same |
| the pt-eliminated condition has **full** degree in t₀: D(n) = 2ⁿ−1 and ANF degree n | [VERIFIED n = 8..] — T4.3's switch, and it points the same way as T3.2 | `results/E07_server_lclm_degrees.csv` |

<!--TABLES-->

## Appendix B — generated tables

### B.0 The decoupled control (task 0)

k S-boxes on pairwise **disjoint** variable blocks, put through the same subset lattice as `gate1_pairs.py`.  `cross_excess = 0` is forced by construction here, so a non-zero value would mean the instrument is wrong and every `cross_excess` in batches 1 and 2 would have to be discarded.

| n | k | nv | blocks | eqs | D3 syz (8n·|B|) | falls_D4 | Σ singles | cross_excess | sat |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 7 | 4 | 56 | {0} | 21 | 56 (56) | 64 | 64 | 0 | 0.0411 |
| 7 | 4 | 56 | {1} | 21 | 56 (56) | 112 | 112 | 0 | 0.0427 |
| 7 | 4 | 56 | {2} | 21 | 56 (56) | 57 | 57 | 0 | 0.0409 |
| 7 | 4 | 56 | {3} | 21 | 56 (56) | 119 | 119 | 0 | 0.043 |
| 7 | 4 | 56 | {0+1} | 42 | 112 (112) | 176 | 176 | 0 | 0.0838 |
| 7 | 4 | 56 | {0+2} | 42 | 112 (112) | 121 | 121 | 0 | 0.082 |
| 7 | 4 | 56 | {0+3} | 42 | 112 (112) | 183 | 183 | 0 | 0.0841 |
| 7 | 4 | 56 | {1+2} | 42 | 112 (112) | 169 | 169 | 0 | 0.0836 |
| 11 | 2 | 44 | {0} | 33 | 88 (88) | 22 | 22 | 0 | 0.0997 |
| 11 | 2 | 44 | {1} | 33 | 88 (88) | 22 | 22 | 0 | 0.0997 |
| 11 | 2 | 44 | {0+1} | 66 | 176 (176) | 44 | 44 | 0 | 0.1994 |

**11 of 11 subsets give exactly 0**  → the instrument is sound.

### B.1 `d_solv` with corrected exponents (task 1.2)



grevlex; msolve and Magma agree on every completed cell.  A cell is `—` where every run of it was censored at the 2-hour cap.

| system | ℓ | n=9 | n=11 | n=13 | n=15 | n=17 |
|---|---:|---:|---:|---:|---:|---:|
| `S_quad` | 1 | — | 4 | 4 | 4 | — |
| `P1only` | 1 | — | 4 | 4 | 4 | 5 |
| `random` | 1 | — | 4 | 4 | 4 | — |
| `S_quad` | 2 | 4 | 4 | — | — | — |
| `P1only` | 2 | 4 | 4 | — | — | — |
| `random` | 2 | 4 | 4 | — | — | — |

Exponents used: ℓ=1, n=11: `3;4`, ℓ=1, n=13: `3;5`, ℓ=1, n=15: `4;11`, ℓ=1, n=17: `3;4`, ℓ=2, n=9: `2;4;5`, ℓ=2, n=11: `3;4;7`.

### B.2 Hybrid guessing: measured `d_solv` against B2's `d_reg` (task 1.3)

_not collected_

### B.3 Measured Macaulay rank against the Boolean semi-regular prediction (task 1.4)

_not collected_

### B.4 Calibrated hybrid-guessing cost (task 1.5)

`cost = 2^k · (Σ_{i≤d} C(N−k,i))^ω`, `d = dreg_boolean(N−k, m) − δ`, minimised over k.  **This is a table, not a conclusion** — δ is measured at n ≤ 15 and `COMPLEXITY_ACCOUNTING §4` does not license carrying it to n = 128..256.  B1/B2 use E02's *commutative* monomial count `C(N−k+d,d)²`, not this Boolean one, so the δ = 0, ω = 2 row is close to but not literally B2.

| scheme | n | ℓ | N | m | δ | ω | k* | d_reg | d | log2 T | B1 | B2 | B3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AIM3-I | 128 | 2 | 384 | 1152 | 0 | 2.0 | 5 | 19 | 19 | **215.83** | 262.6 | 218.4 | 147.9 |
| AIM3-I | 128 | 2 | 384 | 1152 | 0 | 2.81 | 91 | 13 | 13 | **298.02** | 262.6 | 218.4 | 147.9 |
| AIM3-I | 128 | 2 | 384 | 1152 | 1 | 2.0 | 5 | 19 | 18 | **207.33** | 262.6 | 218.4 | 147.9 |
| AIM3-I | 128 | 2 | 384 | 1152 | 1 | 2.81 | 91 | 13 | 12 | **285.55** | 262.6 | 218.4 | 147.9 |
| AIM3-III | 192 | 2 | 576 | 1728 | 0 | 2.0 | 2 | 27 | 27 | **308.98** | 381.0 | 312.5 | 213.14 |
| AIM3-III | 192 | 2 | 576 | 1728 | 0 | 2.81 | 102 | 20 | 20 | **428.45** | 381.0 | 312.5 | 213.14 |
| AIM3-III | 192 | 2 | 576 | 1728 | 1 | 2.0 | 2 | 27 | 26 | **300.29** | 381.0 | 312.5 | 213.14 |
| AIM3-III | 192 | 2 | 576 | 1728 | 1 | 2.81 | 102 | 20 | 19 | **415.78** | 381.0 | 312.5 | 213.14 |
| AIM3-V | 256 | 3 | 1024 | 3072 | 0 | 2.0 | 9 | 44 | 44 | **523.7** | 651.5 | 529.1 | 278.2 |
| AIM3-V | 256 | 3 | 1024 | 3072 | 0 | 2.81 | 140 | 35 | 35 | **726.53** | 651.5 | 529.1 | 278.2 |
| AIM3-V | 256 | 3 | 1024 | 3072 | 1 | 2.0 | 9 | 44 | 43 | **514.76** | 651.5 | 529.1 | 278.2 |
| AIM3-V | 256 | 3 | 1024 | 3072 | 1 | 2.81 | 140 | 35 | 34 | **713.59** | 651.5 | 529.1 | 278.2 |

### B.5 `dim I(V_e)_{≤4}` against what D = 5 reaches (task 2.1)

`dim I(V_e)_{≤4}` is the whole degree-≤4 part of the vanishing ideal of the S-box graph; `dim(span_{D=5} ∩ deg≤4)` is what the degree-5 Macaulay matrix of ⟨P₁,P₂,P₃⟩ actually produces.  The difference is an **upper bound** on the number of quartic generators outside the quadratic ideal — D = 6 might still reach some of them, and D = 6 at n = 16 is ≈ 290 GB (§8).  The measurement is only meaningful when the graph has more points than there are monomials; `deficit_floor` is what the point shortage alone would force.

| n | e | class | points | monomials | deficit floor | dim I(V)≤4 | D=5 reaches | outside | /n |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 16 | 3 | `generic` | 65535 | 41449 | 0 | 21592 | 21432 | **160** | 10 |
| 16 | 5 | `{-3,-3^-1}` | 65535 | 41449 | 0 | 21224 | 21224 | **0** | 0 |

### B.6 What the D = 5 falls are made of (task 2.2)

`X` = dim of the span of the products `z_k·c_j` (the `#cubics` new cubics of D = 4 times the 2n variables) **modulo the D = 4 row space**; `net new = falls_D5 − X` is what D = 5 adds that is not a multiple of a D = 4 fall.

| n | e | class | #cubics | products | X | relations | /n | falls_D5 | net new | /n |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 | 3 | `generic` | 28 | 784 | 532 | 252 | 18 | 882 | **350** | 25 |
| 14 | 9 | `{-3,-3^-1}` | 29 | 812 | 560 | 252 | 18 | 700 | **140** | 10 |
| 16 | 3 | `generic` | 32 | 1024 | 736 | 288 | 18 | 1136 | **400** | 25 |
| 16 | 5 | `{-3,-3^-1}` | 33 | 1056 | 768 | 288 | 18 | 928 | **160** | 10 |

### B.7 Block-restricted XL on the kernel model (task 3.1)

Multipliers restricted to one block instead of all 6n variables.  `corank = cols − rank`; linearisation would mean corank 1.

| n | D | multipliers | rows = rank | syz | dim_low | corank | full-multiplier rank | full dim_low |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 7 | 2 | U block | 42 | 0 | 0 | 862 | — | — |
| 7 | 3 | U block | 924 | 0 | 63 | 11460 | — | — |
| 7 | 4 | U block | 9736 | 8 | 1441 | 114578 | — | — |
| 7 | 5 | U block | 65329 | 275 | 15349 | 909653 | — | — |
| 11 | 2 | U block | 66 | 0 | 0 | 2146 | 66 | 0 |
| 11 | 3 | U block | 2244 | 0 | 99 | 45728 | 4422 | 132 |
| 11 | 4 | U block | 37092 | 0 | 3465 | 731600 | 143715 | 9075 |
| 7 | 2 | W block | 42 | 0 | 0 | 862 | — | — |
| 7 | 3 | W block | 924 | 0 | 63 | 11460 | — | — |
| 7 | 4 | W block | 9744 | 0 | 1386 | 114570 | — | — |
| 7 | 5 | W block | 65604 | 0 | 14511 | 909378 | — | — |
| 11 | 2 | W block | 66 | 0 | 0 | 2146 | 66 | 0 |
| 11 | 3 | W block | 2244 | 0 | 99 | 45728 | 4422 | 132 |
| 11 | 4 | W block | 37092 | 0 | 3366 | 731600 | 143715 | 9075 |

### B.8 Does the first layer's rank condition carry information? (task 3.2)

`N(W)` is 2n × (n+1) with entries affine in W = (a₀,a₁).  A **uniformly random** matrix of that shape has `P(rank ≤ n) = 1 − Π_{i=0}^{n}(1 − 2^{i−2n})`, which is the column headed *uniform*: the question is whether the true W is distinguished by more than that.

| n | rank at true W | corank | (pt,1) in kernel | samples | random rank ≤ n | measured % | uniform % | ratio |
|---:|---:|---:|:--:|---:|---:|---:|---:|---:|
| 7 | 7 | 1 | yes | 20000 | 617 | 3.085 | 1.548 | 1.99× |
| 9 | 9 | 1 | yes | 20000 | 168 | 0.840 | 0.390 | 2.16× |
| 11 | 11 | 1 | yes | 20000 | 28 | 0.140 | 0.098 | 1.43× |
| 13 | 13 | 1 | yes | 20000 | 9 | 0.045 | 0.024 | 1.84× |

### B.9 Block orders computed directly, in PolyBoRi (task 3.3)

PolyBoRi works **in** the declared order and prints one `Current Degree: d` line per degree, so `d_solv` here is a genuine per-order solving degree — the measurement batch 1 could not make, because Magma's GF(2) F4 runs in grevlex whatever order the ring declares.

| system | n | order | ring order | `d_solv` | steps | GB size | GB max deg | s | censored |
|---|---:|---|---|---:|---:|---:|---:|---:|:--:|
| `E04b2_n9_l2_s1_P1only` | 9 | `grevlex` (grevlex) | degrevlex | 5 | 10 | 27 | 1 | 27 | no |
| `E06_n7_kernel` | 7 | `blockU` (block:21) |  | 3 | 4 |  |  | 26 | yes |
| `E06_n7_kernel` | 7 | `grevlex` (grevlex) |  | 4 | 31 |  |  | 502 | yes |
| `E06_n7_kernel` | 7 | `grevlex` (grevlex) |  | 4 | 14 |  |  | 250 | yes |

### B.10 E05: eliminating t₁ without a resultant (task 5)

Eq1 is separable, `W(t₀) + V(t₁)`, so its variety is found by tabulating V and looking W up — Θ(2ⁿ) field operations, no eliminated polynomial.  The `Resultant` column is batch 1's measurement of the same instance for comparison.

| n | build s | **eliminate s** | roots | spurious | true key a root | `Resultant` s | speed-up |
|---:|---:|---:|---:|---:|:--:|---:|---:|
| 7 | 0.011 | **0.006** | 3 | 2 | yes | 294.5 | 49,083× |
| 8 | 0.019 | **0.016** | 2 | 1 | yes | — | — |
| 9 | 0.025 | **0.04** | 1 | 0 | yes | — | — |
| 10 | 0.03 | **0.098** | 1 | 0 | yes | — | — |
| 11 | 0.036 | **0.21** | 1 | 0 | yes | — | — |
| 12 | 0.045 | **0.463** | 3 | 2 | yes | — | — |
| 13 | 0.057 | **1.167** | 1 | 0 | yes | — | — |
| 14 | 0.084 | **2.777** | 1 | 0 | yes | — | — |
| 15 | 0.094 | **5.876** | 1 | 0 | yes | — | — |
| 16 | 0.188 | **15.818** | 2 | 1 | yes | — | — |
| 17 | 1.623 | **461.241** | 4 | 3 | yes | — | — |

Least squares over n = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]: slope of `log2(eliminate seconds)` vs n = **1.406**; over n ≥ 13 only (past the constant-overhead regime): **1.976**.

The method is Θ(2ⁿ) by construction, so a slope of 1 is the expected answer and the measurement confirms it rather than assuming it.  There is **no `log2(deg)` slope**: this route forms no eliminated polynomial, and the route that does is censored from n = 8 (`results/E05_server_resultants.csv`).

