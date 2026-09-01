# Server batch 1 — experiments run beyond the batch prompt

**Branch `claude/server-batch-1` · 2026-08-31 · for review before the PR is merged**

This document lists everything in this branch that the batch prompt did **not** ask for, why it was
run, what it cost, and what would have been wrong or missing without it.  The prompt-specified work
is described in `results/REPORT_server_batch1.md`; this file is only about the additions.

## 0. Why there were additions at all

The prompt fixes the seven tasks and the decision rules precisely.  Three situations make "just
follow the instructions" produce a wrong conclusion, and those are the only three that triggered
extra work:

* **The rule's own data requirement cannot be met with the sizes listed.**  Gate 1 demands
  `cross_excess = 0` on ≥ 4 values of n; half the listed sizes are saturated and inadmissible.
* **A measurement could be an artefact.**  A non-zero number with no control cannot distinguish
  "structure" from "more equations reach further".
* **The tool may not do what the plan assumes.**  Whether Magma actually runs F4 in a declared
  block order, for instance.  Without checking, one reports a measurement that never happened.

Each item below gives: what, why, cost, result, and **what would have been wrong without it**.

---

## 1. The subset lattice (`server/jobs/gate1_pairs.py`) — the most important addition

**Prompt.**  `full_system_falls.py --per-block`, i.e. singletons and the whole set only, with
pairwise localisation *only if* `cross_excess` grows.

**Added.**  A general tool measuring **every subset** of the S-box blocks — singletons, pairs,
triples, the whole set — and putting a same-shape random system, cut into blocks at the same
positions, through the identical lattice.

**Where.**  ℓ = 2 at n = 13, 15; ℓ = 3 at n = 11 (default and generic exponent lists), n = 13,
n = 15 (three-first-layer subset only).

**Cost.**  7 jobs, ≈ 8.6 h wall clock, peak 15.4 GB (ℓ = 3, n = 15 triple).

**Result.**

1. It upgrades "direct sum" from *matches at the top* to **exactly additive at every level of the
   lattice** — at ℓ = 2, n = 13 all three singletons, all three pairs and the triple satisfy
   `falls(B) = Σ_{i∈B} falls({i})`.
2. **It overturned a conclusion this batch had already written down.**  At ℓ = 3, n = 11 the
   `cross_excess` of 1639 localises onto the *three first-layer S-boxes* (144) — exactly the
   structure AIM3-V has and AIM3-I/III does not.  The report was written that way.  The **n = 13
   run then gave exactly 0, and n = 15 gave 0**, which is what forced the reading down to a
   small-parameter effect.

**Without it.**  The report would have claimed a cross-S-box degree-fall mechanism specific to the
AIM3-V shape and recommended promoting T2 on that basis.  That would have been wrong.

**Reviewer's note.**  This is the addition most worth keeping.  Its value is not extra data points;
it is the only instrument that can localise a non-zero `cross_excess` onto specific S-box subsets
and thereby expose its instability.

---

## 2. Extra Gate-1 sizes

**Prompt.**  ℓ = 2 at n = 13, 15, 17; ℓ = 1 at n = 13, 15; ℓ = 3 at n = 11, 13.

**Added.**  ℓ = 1 at **n = 17, 19** (3 seeds each) and ℓ = 2 at **n = 19**.

**Why.**  The prompt's own rule admits a run only at saturation < 0.9.  Measured:

| ℓ | n = 13 | 15 | 17 | 19 |
|---:|---:|---:|---:|---:|
| 1 | 0.998 ✗ | 0.978 ✗ | 0.525 ✓ | 0.470 ✓ |
| 2 | 0.455 ✓ | 0.398 ✓ | 0.349 ✓ | 0.313 ✓ |

**Both of the ℓ = 1 sizes the prompt lists are saturated** — not one of them can enter the verdict.
ℓ = 2 needed n = 19 to reach four admissible sizes (13, 15, 17, 19).

**Cost.**  ≈ 1 h for ℓ = 1, plus the ℓ = 2 n = 19 seeds (≈ 1 h each; the D = 4 matrix is
282834 × 425924 ≈ 15 GB).

**Result.**  19 of 22 non-saturated runs give `cross_excess` exactly 0, and ℓ = 2 reaches the four
sizes `COMPLEXITY_ACCOUNTING §4` requires.

**Without it.**  ℓ = 1 would contribute nothing and ℓ = 2 would have three sizes — below the
extrapolation threshold, so Gate 1 would still be [CONJ].

---

## 3. Generic-exponent re-run of ℓ = 3 at n = 11

**Added.**  The same lattice with `--exps 3,4,7,8` (all generic or `{−3}` class) against the
default `2,3,4,5`.

**Why.**  `default_exps` always starts at e = 2, which is the degenerate `{2, 2⁻¹}` class (7n
falls).  The anomaly could have been an exponent effect.

**Cost.**  ≈ 67 min.

**Result.**  Under the default exponents the pair `{0,1}` carries an excess of 11 = n; with generic
exponents **every pair is exactly 0**, and the first-layer triple drops from 144 to 77 = 7n.  A
quantity that changes with the exponent list is not a cross-S-box mechanism — this is the second
independent reason the n = 11 row is not believed.

---

## 4. Subtracting the *ideal* of the image criteria (`server/jobs/e12_criteria_subtract.py`)

**Prompt.**  Run `explain_n16.py 16 --deg 4`; if the script has no `--deg`, "subtract the univariate
image criteria the way its README prescribes".  `--deg` was added, so the letter of the task was
already met.

**Added.**  A separate computation that builds each criterion from its ANF by Möbius transform (so
its Boolean degree is **measured**, not assumed), spans the entire degree-≤ 4 part of the ideal
those criteria generate, and subtracts that.

**Why.**  `explain_n16.py` decides "univariate" by asking whether a relation lies in the span of
monomials touching only t₀ or only t₁.  That is right at degree 3 and **wrong at degree 4**: if
g = Tr(t₁^{−r₁}) is a cubic relation, every product `g·t₀ᵢ` is a quartic relation touching both
blocks which is not a cross relation at all, only a multiple of a univariate one.

**Cost.**  Three runs, ≈ 8 min total.

**Result.**  At n = 16 the quartic relation space has dimension 33.  Read with the script's own
test: 17 on t₁ only, 0 on t₀ only → **"16 genuine cross-S-box quartics"**.  t₀ has exactly 16
variables, and those 16 are precisely `g·t₀ᵢ`.  Done properly: `dim R = 33`, `dim C = 33` →
**0 genuine cross-S-box relations**.  Claim K strengthens from degree 3 to degree 4.

**Without it.**  The report would have claimed 16 genuine cross-S-box quartic relations — wrong, and
wrong in the direction that would be misread as support for T2.

**Also run.**  n = 16 at deg = 3, to reproduce the known Tier-A answer (1 − 1 = 0); and n = 14 at
deg = 4, which showed the measurement is **meaningless at that size** (16382 graph points against
24158 monomials forces ≥ 7776 relations for any instance whatsoever).  No positive result, but it
documents a trap.

---

## 5. E03 at D = 5: a random control and a second exponent class

**Prompt.**  Run D = 5 at n = 14 and 16 and report whether new falls to degree ≤ 4 appear.

**Added.**  At each n, (a) a same-shape random quadratic control, and (b) a second exponent from
the `{−3}` class rather than only a generic one.

**Cost.**  ≈ 3 h (inside the 54 min of `E03_D5_n14` and the 4 h 59 min of `E03_D5_n16`), peak
8.7 GB.

**Result.**

| n | e | class | falls_D5 | / n | saturation |
|---:|---:|---|---:|---:|---:|
| 14 | 3 | generic | 882 | **63 = 4n+7** | 0.583 |
| 16 | 3 | generic | 1136 | **71 = 4n+7** | 0.517 |
| 14 | 9 | {−3} | 700 | **50 = 4n−6** | 0.576 |
| 16 | 5 | {−3} | 928 | **58 = 4n−6** | 0.512 |
| 14 | — | random control | 7966 | 569 | **1.000 (saturated, useless)** |
| 16 | — | random control | **0** | 0 | 0.584 |

1. **The control only becomes usable at n = 16.**  At n = 14 the random system is *fully saturated*
   at D = 5, so its 569n is an artefact and no baseline.  At n = 16 it is not saturated (0.584,
   comparable to the S-box's 0.517) and its `falls_D5` is **exactly 0 for both exponents** — which
   is what makes the D = 5 falls structural rather than a space-filling effect.
2. **A candidate closed form**: all four runs land exactly on `n(4n+7)` for generic e and `n(4n−6)`
   for the `{−3}` class.  Flagged as a **candidate, not a rule** (two sizes per class), and
   deliberately **not** written into `src/aim3/structures.py`.
3. **A sign flip worth a theorist's attention**: the `{−3}` class gets *one more* cubic than generic
   at D = 4 and *13n fewer* quartics at D = 5, at both sizes.

**Without it.**  D = 5 would be two isolated numbers with no control to say whether they mean
anything, and neither the closed form nor the sign flip would be visible.

---

## 6. E05 at n = 7

**Prompt.**  n = 8..20.

**Added.**  After n = 8 hit the 2-hour cap and the escalation stopped as the rules require, n = 7
was run.  n = 8 was also re-run under a 2-minute cap purely to record the system's shape on the
censored row (otherwise that row says only "censored" and the reader cannot tell what was attempted).

**Why.**  Stopping as instructed left E05 with **zero completed sizes**, so the slope task 7 asks
for had nothing at all to stand on.  A smaller size gives at least a one-step lower bound.

**Cost.**  13 min 50 s + 2 min.

**Result.**  n = 7 completes: resultant 294.5 s, degree 22204 (98 % of the Bézout-style bound
22676 — **no degree collapse**), 3 roots including the true key (2 spurious), and it
**factors into 11 irreducibles** of degrees 1, 1, 1, 10, 145, 172, 1250, 3931, 4745, 5454, 6494.
That is the **T3.3 lead**: it does factor, but the largest factor is still 6494 ≈ 2^{1.86n} and the
factorisation costs 534 s, nearly twice the resultant.

One step from n = 7 to n = 8 multiplies the time by at least 24, so `log2(time)` rises by ≥ 4.6 over
one unit of n — clearing the ROADMAP §3.3 slope-2 criterion **as a one-step lower bound**.  The
report therefore quotes **no slope**.

**Without it.**  E05 would be a single blank censored row: no factorisation result, no growth
information at all.

---

## 7. The E10 `genref` baseline and two Magma behaviour probes

**Prompt.**  Block orders `elim` (pt block first / t block first) and, for the kernel model, U > W
and W > U.

**Added (a): a `genref` baseline** — grevlex on the **same general polynomial ring**.  Magma 2.24's
`BooleanPolynomialRing` accepts an order **name only**; an `elim` order raises
`Bad argument types`.  Block orders therefore have to run on
`PolynomialRing(GF(2), n, "elim", S)` with the field equations added by hand — a **slower engine**.
Comparing a block-order timing against the Boolean-ring grevlex would show block orders roughly
1000× slower for the wrong reason: that is the ring, not the order.

**Added (b): two probe scripts** (seconds each):
* whether `BooleanPolynomialRing` accepts an elim order → **no**;
* whether Magma honours a declared elimination order → **the returned basis genuinely is in the
  elim order** (the same toy system gives 30 elements under grevlex and 15 under elim), **but
  `SetVerbose("Faugere",2)` reports `Order: Graded Reverse Lexicographical` in both cases** —
  Magma's GF(2) F4 phase always runs in grevlex and converts afterwards.

**Cost.**  `E10_blockorder_E06_kernel` ≈ 2 h; the probes seconds.

**Result.**  **The "F4 step degrees under a block order" that task 6 / E10 asks for cannot be
obtained from Magma 2.24 at all** — such a log is a grevlex log.  This is a tool limitation, not an
implementation problem.  What remains comparable is total time and the shape of the resulting
basis, and the E10 section was rewritten accordingly.

The kernel model is the one informative instance (it has no unique solution, so its basis *is*
order-dependent): `genref` 131.5 s / basis 4432, `blockU` 485.3 s / 2357, `blockW` 1127.0 s / 2772 —
**block orders give a smaller basis for 3.7× and 8.6× the time**, all of it in the order conversion.

**Without it.**  The report would have presented a block-order step-degree comparison Magma never
performed, and would have read a ring difference as an order effect.

---

## 8. The `--no-points` fast path for E03 at n = 19, 20

**Added.**  A `--no-points` option on `sbox_degree_fall.py`, used to re-sweep n = 19 and 20.

**Why.**  The point-set matrix is O(2ⁿ) — at n = 20 it is 1048575 × 10701, ≈ 5.7 GB and ≈ 9 min to
export *per exponent* — and **the deliverable of task 4, the fall-rule mismatch count, does not need
it at all**.  At full speed n = 19 (16 exponents) and n = 20 (6 exponents) could not finish inside
the batch.

**Cost.**  n = 19 ≈ 23 min, n = 20 ≈ 14 min (against 10 min *per exponent* the other way).

**Result.**  **All 40 admissible (n, e) at n = 17..20 completed, 0 mismatches.**  The price is that
those rows have empty `IV_deg3` and `outside_ideal`; `collect.py` keeps the more complete row when
a pair is duplicated, and `e03_summary_merged.sh` excludes the unparsable rows from `summarise.py`
and reports how many it excluded.

**Without it.**  n = 19 would be 5/16 and n = 20 2/6 — task 4's "every admissible e" unmet.

---

## 9. Two unit tests and the field-modulus extension

* **`src/aim3/legacy_fields.py`** — `legacy/aim3_lab.CONWAY` stops at n = 16 and every script raises
  `KeyError` above it.  This adds the n = 17..24 moduli from `gf2n.TOY_MODULI` **after verifying
  each is irreducible by Rabin's test** and asserting that entries shared with the legacy table are
  unchanged.  `legacy/` itself is untouched.  Not an "extra experiment" — **nothing above n = 16
  runs without it**.
* **Two unit tests** (6/6 → 8/8).  `test_export_matches_boolmac` checks row for row that the
  exported matrix is the one `boolmac.py` streams and reproduces the Sage backend's column
  reversal; `test_legacy_field_table_extended` checks the added moduli.  Reason: **a wrong column
  order does not raise, it returns a plausible number** — the one place the new backend can fail
  silently.

---

## 10. Additional work that produced nothing (stated plainly)

* **The `random` sub-lattice inside the subset-lattice jobs** (inside ℓ = 2 n = 13/15 and ℓ = 3
  n = 11, ≈ 2–3 h): a random system has **no degree falls at all**, so its `cross_excess` is
  trivially 0 and it **cannot discriminate**.  Said so in report §2.3, and the control that would
  actually settle the ℓ = 3 question — four S-boxes on *disjoint* variable blocks — is listed as not
  done.  This compute was wasted.
* **`E12_criteria_n14_deg4`** (1 min): shows the measurement is meaningless at that size.  No
  positive result; records a trap.
* **`E05_resultants_n8shape`** (2 min): pure bookkeeping, to fill the censored row's shape columns.
* **One operator error.**  Stopping a superseded job with `pkill -f "gate1_pairs.py 13 -l 3"` also
  matched a target job that had just been launched and killed it; it was relaunched, **losing about
  20 min**.  The drivers now kill by process group (`run_f4_degrees.py`, `run_resultants.py`), but
  that particular kill was a hand-typed command, not a script fault.
* **Prompt-specified jobs stopped on purpose.**  E10's block orders at n = 13, 15 (once the tool
  limitation was established they could not answer the question they were launched for, and each
  held ≈ 9 GB the Gate-1 jobs needed); and the full ℓ = 3 n = 13 lattice (superseded by a targeted
  triple job).  Both are in report §8.

---

## 11. Totals, and what to keep

**Additional compute:** ≈ 14–15 h wall clock (mostly in parallel), peak 15.4 GB, all inside the
per-tier caps of `server/README.md`.

**The additions that changed a conclusion** — deleting any of these makes the report wrong or leaves
it short of the prompt's own decision threshold:

| addition | consequence of not running it |
|---|---|
| subset lattice at ℓ = 3, n = 13 and 15 | would claim a cross-S-box mechanism in AIM3-V — **wrong** |
| Gate 1 at ℓ = 1 n = 17/19 and ℓ = 2 n = 19 | ℓ = 1 contributes nothing; ℓ = 2 has 3 sizes, below the ≥ 4 rule |
| ideal-subtraction of the image criteria | would claim 16 genuine cross-S-box quartics at n = 16 — **wrong** |
| D = 5 random control at n = 16 | no way to tell whether the D = 5 falls are structural |
| E05 at n = 7 | zero completed sizes, no factorisation result |
| E10 `genref` baseline + the two probes | would report a block-order comparison Magma never performed |
| E03 `--no-points` | task 4's "every admissible e" unmet at n = 19, 20 |

**Recommendation.**  Keep all seven — each either corrected a wrong conclusion or was necessary to
reach a threshold the prompt itself sets.  The items in §10 may be dropped or marked
"for the record" at the reviewer's discretion; leaving their CSV rows in place is harmless, because
the report already states that they carry no discriminating power.

**Still for a human to decide** — unchanged, and listed in report §9: whether to chase the ℓ = 3
n = 11 outlier, whether the 0.9 saturation threshold is strict enough, how the `d_solv` negative
result should be written up, whether E10 should be re-scoped, and which `paper/claims_ledger.md`
rows to apply.
