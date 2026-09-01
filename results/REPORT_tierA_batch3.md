# Tier-A batch 3 — report

Branch `claude/tier-a-batch-3`, base `main` at `11ec8ab` (server batch 1 merged).  Pure Python, standard library only;
`bash scripts/smoke_test.sh` and `python3 tests/run_tests.py` (17/17) pass on the final tree.

Write scope respected: `docs/`, `paper/`, `src/aim3/` (one new file + one backward-compatible method), `tests/`,
`results/README.md`, `experiments/E07_ore_skew_elimination/`, `experiments/E09_multitarget_weak_iv/`,
`results/E07_theory_*.csv`, `results/E09_theory_*.csv`, this report.  **Not touched:** `legacy/`, `server/`,
`results/server/`, `src/aim3/boolmac.py`, and the E03/E04/E06/E10/E12 scripts (the server branch owns them).

---

## 1. Sync tidy (task 1)

* `paper/claims_ledger.md`: the §9 suggestions of `REPORT_server_batch1.md` are in — G and H widened to n ≤ 20, H3 turned
  to [OPEN] (reproduced at n = 17, 19), K/L/L2/L3 widened, M2 restated as T1.3 with the ℓ = 1,2,3 subset lattice at
  n = 13..19, M3 extended to n = 13; new rows M5 (the ℓ=3 n=11 `cross_excess` excess, [OPEN]), U (measured `d_solv`
  with the e = 2 caveat) and V (block orders as a measurement limitation).  Batch 3 then added Q2, W, W2, W3 and
  restated Q (below).  40 rows now.
* `results/README.md`: the five Tier-A batch-2 `E0*_theory_*.csv`, all twelve server batch-1 `E*_server*.csv`, both batch
  reports and (batch 3) `E07_theory_skew_degrees.csv` are documented.  Every file in `results/` now has a row.
* `docs/ROADMAP.md` §3.1: the Thm T1.2 candidate mechanism is marked **refuted** and the target restated per
  `THEORY_NOTES.md` §2.6–2.7 — the new cubics split n + n, both halves need weight-2 multipliers on P₂/P₃
  (`x^{2^i}y^{2^j}·P₂` and `y^{2^i+2^j}·P₃`), and the single extra element of the {−3} and {n/2−1} classes is a constraint
  any description must explain.
* `docs/ROADMAP.md` §3.2: the two server facts added — Magma's GF(2) F4 computes in grevlex whatever order is declared
  (T2.2 needs PolyBoRi), and the kernel model is ~1000× slower than `S_quad` in F4 but far faster than a same-shape
  random bilinear system.
* `docs/COMPLEXITY_ACCOUNTING.md`: new rule 5 — toy exponents must come from the generic class
  (`degenerate_class(n,e) == []`) unless the degenerate class is the point of the experiment; and server batch 1's
  `d_solv` table is unusable for extrapolation because the default exponent list starts at e = 2.
  (Old rules 5/6/7 renumbered 6/7/8.)

Commit `4de9131`.

## 2. The Ore prototype (task 2.1)

`src/aim3/ore.py` — `SkewPoly` over R = F_{2^n}⟨φ⟩ with φ·a = a²·φ: add, multiply, `apply` (the linearised action),
`scale`/`monic`, `root_space`, `divmod_right`, `gcrd`, `lclm`, `lclm_cofactors`, plus `solve_gf2n`, `linearise_values`
(Moore-matrix solve) and `from_bit_matrix` (legacy `apply_mat` bit matrix → SkewPoly).

The module docstring keeps two things apart that are easy to conflate and that the measurement depends on: **the ring**,
where φⁿ ≠ 1, degrees add and the Euclidean algorithm lives; and **the action** on F_{2^n}, which factors through
R/(φⁿ−1) — so a degree-d element has a root space of dimension ≤ min(d, n), not d.

`tests/test_ore.py` — the four required tests all pass, plus three more:

| test | result |
|---|---|
| (a) `(f*g).apply(v) == f.apply(g.apply(v))` and associativity | pass, n = 5,7,8,11, 20 random triples each |
| (b) `lclm(f,g)` right-divisible by f and by g | pass, n = 5,7,11, 25 pairs each |
| (c) `root_space(gcrd(f,g))` = intersection of the root spaces, **exhaustive over F_{2^n}**, n ≤ 10 | pass, n = 5,7,9,10; half the cases have a planted common right factor so the test is not vacuous (asserted) |
| (d) `deg lclm = deg f + deg g − deg gcrd` | pass, n = 5,7,11,13 |
| `f = q*g + r` with deg r < deg g | pass |
| `from_bit_matrix` vs legacy `apply_mat` on a unit-LU matrix | pass, n = 5,7,9,11 |
| `linearise_values` interpolates an arbitrary F_2-linear map | pass |

`tests/run_tests.py` now discovers every `tests/test_*.py` module instead of only `test_core.py`.  Commit `07fb752`.

## 3. T4: what is left after pt is eliminated (task 2.2) — the main result

`experiments/E07_ore_skew_elimination/skew_elimination_probe.py`, output `results/E07_theory_skew_degrees.csv`
(78 rows), README with the full method.  Total runtime ≈ 2 min.

Set-up: `Λ_i(t_i) = (φ^{e_i} + t_i)·A_i`, `Eq_i: Λ_i(t_i)(pt) = c_i(t_i) = 1 + b_i^{2^{e_i}} + t_i b_i`.
`Λ_i(pt) = c_i` is asserted on every instance built.

### 3.1 A correction to the brief's hypothesis

The brief conjectured: *at the true (t₀,t₁) the GCRD has degree ≥ 1 because pt and pt+a lie in the intersection of the two
root spaces; at random (t₀,t₁) the GCRD is 0.*  **That is false**, and the reason matters:
`ker(φ^e + t) = {0, t^r}` is one-dimensional, so `ker Λ_i = A_i^{-1}({0, a_i})` is a 1-dimensional F_2 space, and two
1-dimensional subspaces of an n-dimensional space meet only in 0 unless they coincide (probability ≈ 2^{-n}).  Nothing
about the true key forces it, because **the true key satisfies affine equations, not homogeneous ones**.

Measured (`--probe`, n = 6..14, 32 random (t₀,t₁) per n): `deg gcrd = 0` and `dim(ker Λ₀ ∩ ker Λ₁) = 0` at the true key
in every single run — indistinguishable from random.  The statistic that *does* separate the true key is the affine one:
the system {Λ₀(v) = c₀, Λ₁(v) = c₁} is solvable at the true (t₀,t₁) for every n, and for only 1 of 288 random pairs
(one accidental hit at n = 9, consistent with the ≈ 2·2^{-n} chance).

### 3.2 The degree measurement — exact, not interpolated

t₁ is fixed at its true value and t₀ runs over **all** of F_{2^n}, so nothing here is an interpolation lower bound:
the Boolean (ANF) degree comes from a Möbius transform over all 2^n values, the univariate degree from
`c_j = Σ_{a≠0} E(a)a^{-j}` (every map F_{2^n}→F_{2^n} is a polynomial of degree ≤ 2^n−1).  `--selftest` checks both
instruments against closed-form degrees before any measurement is trusted.

| n | e₀,e₁ | deg Λ | deg lclm (in R) | **Boolean deg E** | **univariate deg E** | Boolean deg C | \|C\| | \|Z(E)\| |
|---|-------|-------|-----------------|-------------------|----------------------|---------------|-----|--------|
| 6 | 5,5 | 10,10 | 20 | 5 | 62 (of 63) | 5 | 2 | 2 |
| 7 | 6,6 | 12,12 | 24 | 7 | 127 | 6 | 2 | 2 |
| 8 | 7,7 | 14,14 | 28 | 8 | 255 | 7 | 2 | 2 |
| 9 | 8,8 | 16,16 | 32 | 9 | 511 | 8 | 2 | 2 |
| 10 | 9,9 | 18,18 | 36 | 10 | 1023 | 9 | 2 | 2 |
| 11 | 3,4 | 13,14 | 27 | 11 | 2047 | 10 | 2 | 2 |
| 12 | 11,11 | 22,22 | 44 | 12 | 4095 | 11 | 2 | 2 |
| 13 | 3,5 | 15,17 | 32 | 13 | 8191 | 12 | 2 | 2 |
| 14 | 3,5 | 16,18 | 34 | 14 | 16383 | 13 | 2 | 2 |

E is the LCLM-eliminated condition `u(c₀(t₀)) + v(c₁(t₁))`; C is the solvability predicate.  Every entry is an
**exact value**.  `Z(E) = C` in every run (mutual inclusion checked both ways), so the condition derived in R — where
φⁿ ≠ 1 — loses nothing against the action on F_{2^n}: **no separate Frobenius-consistency step was needed** (question 3
of the brief).  Exactly one t₀ per instance makes `gcrd(Λ₀,Λ₁)` nontrivial in R; it is never the true one.

**The brief's n = 12 target was met and exceeded: the measurement runs to n = 14** (and would run further; the cost is
2^n LCLM computations, 61 s at n ≤ 14).

Controls:
* `--control` (wrong t₁): |C| ∈ {0, 2} — 0 at n = 7, 9, 13 and 2 elsewhere — and the true t₀ is **never** in C.
* `--joint` (Boolean degree in all 2n bits of (t₀,t₁), n = 6..11): |C| is 2ⁿ−2 or 2ⁿ−3 out of 2^{2n}, Boolean degree
  2n or 2n−1 out of a possible 2n.

### 3.3 Conclusion: **exponential**, and it is proved, not just measured

> **Lemma.** Fix t₁ ≠ 0 with gcd(e₁,n) = 1 and A₁ invertible.  The homogeneous part of Eq₁ has a 1-dimensional kernel,
> so its solution set in pt is empty or a coset {pt, pt+w} of size 2 with A₁(w) = a₁ = t₁^{r₁}.  Each candidate pt gives
> `x₀ = A₀(pt)+b₀` and hence at most one `t₀ = S[e₀](x₀)`.  Therefore **|C| ≤ 2**.  The indicator of {α,β} ⊆ F_2^n with
> α ≠ β has ANF degree exactly n−1 (the weight-(n−1) coefficient omitting bit i is [α_i ≠ β_i]).  Hence
> **deg_ANF C = n−1 for every n**.

The experiment confirms the hypotheses hold in practice and that the prediction is sharp: |C| = 2 at every n = 6..14 and
the second element is exactly the predicted `S[e₀](A₀(pt+w)+b₀)`; with a wrong t₁ the count is 0 or 2, as the lemma allows.

So the answer to **T4.3's trichotomy is exponential** for the LCLM route: D(n) = 2ⁿ−1 univariate, n Boolean.  T4.1's
conjecture "total degree ≈ n in the conjugates" is numerically right and operationally empty — n is the *largest* Boolean
degree available on n bits.  The skew resultant of the first layer is the indicator of a 2-element set: informationally
equivalent to guessing t₀ and testing it, which is what one wants to avoid.

**Status tags.** [PROVED] |C| ≤ 2 and deg_ANF C = n−1, for all n.  [VERIFIED n ≤ 14] the exact degrees of E, |C| = 2,
the identity of the second element, Z(E) = C, and the trivial GCRD at the true key.  [OPEN] whether another elimination
order escapes the collapse — see §6.  **No extrapolation is made from these curves**; every number is an exhaustive count
at that n.  Note also that ROADMAP §3.4's kill criterion is literally about **approximant-basis row degrees**, which
Tier A did not measure; what is measured here is the LCLM elimination, which is one specific approximant.

Commit `a045220`.

## 4. Reused key, q ≤ 6 (task 3)

`experiments/E09_multitarget_weak_iv/reused_key_toy.py` (new `--exps`, `--single-only`; a budget-censored D=3 pass is now
reported as NOT MEASURED instead of recording a lower bound as a value).  All runs below have saturation ≤ 0.29.

**The law holds** for exponent sets whose single instance sits on the (ℓ+1)·8n baseline —
`syz(q) = q·syz(1) + C(q,2)·n`, seed-independent:

| n | exps | syz(1) | q = 2 | 3 | 4 | 5 | 6 |
|---|------|--------|-------|---|---|---|---|
| 9 | 2,4,5 | 216 | (9) † | 27 | 54 | 90 | 135 |
| 11 | 2,3,4 | 264 | 11 | 33 | 66 | 110 | — |
| | | | C(q,2)·n → | 27 / 33 | 54 / 66 | 90 / 110 | 135 |

† the n = 9, q = 2 run is at saturation 0.918 (batch-2 data); shown for continuity, not counted.

**And it is not universal.** Following the new generic-exponent rule turned up a second regime whose cause is one
exponent: **e = n−1**.  Any exponent set containing it (in any position) gives the single ℓ=2 instance exactly one extra
D=3 syzygy and exactly 3n−1 extra new degree-≤2 elements — verified at n = 9, 11, 13, against sets without it that sit
exactly on the baseline.  Then the coupling excess becomes `n·[5·C(q,2) + (q−1)]`:

| n | exps | q=2 | 3 | 4 | 5 | 6 |
|---|------|-----|---|---|---|---|
| 11 | 3,4,10 | 66 | 187 | 363 | 594 | 880 |
| 9 | 4,5,8 | 54 | 153 | 297 | — | — |
| | predicted | ✓ | ✓ | ✓ | ✓ | ✓ |

This is a **different statistic** from the single-S-box D=4 fall count that `structures.degenerate_class` classifies —
that function calls e = n−1 generic, correctly, for its own statistic.  No recommended AIM3 parameter set uses e = n−1
(spec Table 2), so this is a toy-parameter trap, not a property of AIM3.  **[OPEN]** where the coefficient 5 comes from
and whether it depends on n — only two values of n are available at Tier A.

Interpretation, unchanged in substance from batch 2: reused-key coupling is real but **quadratically small** — C(q,2)·n
extra syzygies against q·(ℓ+1)·8n present ones, a relative gain of O(q/n).  `d_solv` is not measured; that needs D ≥ 4.

Commit `1eefcd4`.

## 5. Not done / limitations

* `d_solv` for the reused-key system (needs D ≥ 4 — Tier B/C).
* Approximant / weak-Popov bases over F_{2^n}[t]⟨φ⟩ with t kept symbolic: the literal object of ROADMAP §3.4's kill
  criterion.  Tier A measured the LCLM elimination instead; the two agree on this instance shape but the row-degree
  question is formally still open.
* The last-layer relation `Λ_ℓ(pt)(t₀,t₁) = c_ℓ(pt)` is not part of the probe; only the two first-layer relations are.
* n ≤ 10 and n = 12 have only one generic exponent (n−1), so both S-boxes there share it — not an AIM3 parameter shape.
  n = 11, 13, 14 use two distinct generic exponents and give the same degrees, so the conclusion does not rest on it.
  (These E07 instances do contain e = n−1 at n ≤ 10, 12; §4's anomaly is about a different statistic — the ℓ=2 D=3
  syzygy count — and the E07 degrees at n = 11, 13, 14, where e ≠ n−1, are identical in form.)
* Commits are one per task except that the E09 run data landed in the E07 commit `a045220` (a single `git add -A`);
  the E09 write-up is `1eefcd4`.

## 6. For the server / next batch

1. **Sage cross-validation of E07.**  `experiments/E07_ore_skew_elimination/instance_n11.json` (n = 11, e = (3,4)) plus
   `sage/skew_probe.sage`.  The JSON carries `n`, `modulus`, `e`, `A_bitmatrix`, `A_skew_coeffs`, `b`, `x`, `t`, `c`,
   `pt` with the bit convention documented in the file.  Sage must reproduce: `Λ_i(pt) = c_i`;
   `deg right_gcd(Λ₀,Λ₁) = 0`; `deg left_lcm = deg Λ₀ + deg Λ₁`; and exactly two solvable t₀, one of them `t[0]`.
   Regenerate any n with `--export FILE.json --n N --no-csv`.
2. **The open T4 question.**  Compute an approximant / weak-Popov basis over `F_{2^n}[t]⟨φ⟩` (t symbolic throughout) for
   n = 8..20 and measure its row degrees — that is the literal kill criterion.  Tier A's result predicts they blow up,
   but it does not prove it.
3. **Re-run the `d_solv` table** with generic exponents (COMPLEXITY_ACCOUNTING rule 5) — already on the server's list.
4. **e = n−1 at larger n.**  `reused_key_toy.py N -l 2 --single-only --exps ...` at n = 15, 17, 19 decides whether the
   +1 syzygy / 3n−1 new elements and the coefficient 5 are n-independent.
5. **Reused key at D = 4** for n = 9, 11, q = 2, 3 (the `d_solv` question the ROADMAP §3.6 actually asks).
