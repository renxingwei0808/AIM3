# Server batch 1 — report

**Branch `claude/server-batch-1`, 2026-08-31.  Tier B/C, driven from the WSL laptop over SSH.**

Everything ran on `luck-Precision-7920-Tower` (2× Xeon Gold 6230R, 104 threads, 251 GB) inside
`/home/luck/rxw/AIM3_remote/`.  Software, versions and the traps that cost time are in
`server/ENVIRONMENT.md`; job logs and their `/usr/bin/time -v` stamps are under
`results/server/20260831_*/`.  Every table below is generated from the committed CSVs by
`python3 server/jobs/summarise_batch1.py`, so the report and the data cannot drift apart.

**What this batch was for.**  Tier-A batch 1 left Gate 1 at *[VERIFIED n=11, 14; CONJ in general]*
because every full-system D=4 run in pure Python saturated; it left the calibrated fall-count rule
as a fit to n ≤ 16 with no test outside that range; and it had no measurement of the one observable
saturation cannot destroy — the F4 solving degree.  This batch closes the first two and opens the
third.

---

## 1. Environment and backend alignment

`bash server/setup_env.sh` built the Tier-B/C toolchain under `/home/luck/rxw/AIM3_remote/opt`
(the account has no passwordless sudo, so nothing is installed system-wide):

| tool | version | how |
|---|---|---|
| Magma | **V2.24-2** | pre-existing, licensed |
| SageMath | **10.7** | conda-forge via micromamba (the system Sage is 9.5, below `server/README.md`'s requirement) |
| msolve | **0.10.1** | built from source; needs FLINT ≥ 3, so it is built against the FLINT 3.3.1 inside the sage10 env (system FLINT is 2.8.4 and a locally built 2.9.0 is still too old) |
| GNU parallel | 20260822 | built from source |
| M4RI | system 0.0.20200125 + the copy inside Sage | reached through Sage's dense `GF(2)` matrices |

`bash scripts/smoke_test.sh` ends in `smoke test OK` and `python3 tests/run_tests.py` reports
**6/6**, both on the server.

**The fast rank backend (task 2).**  `src/aim3/boolmac_export.py` writes the *same* Macaulay matrix
that `boolmac.py` streams — same rows, same multipliers, same column order — to a sparse-row file;
`experiments/E04_full_system_hilbert/sage/rank_backend.sage` loads it into a dense `GF(2)` (M4RI)
matrix and echelonises once.  The exported column order puts the degree-D monomials last, so the
backend loads column *j* at position `ncols−1−j`: after a left-to-right echelonisation the pivots
inside the leading block are exactly the degree-D ones, which reproduces `boolmac`'s
`dim_low = rank − rank_top` in a single pass.  The same file format and engine also serve the
O(2ⁿ) point-evaluation matrix of `relation_dim()`.

`--backend sage` was added to `sbox_degree_fall.py` (with `--no-points`), `full_system_falls.py`,
`first_layer_ideal.py`, `kernel_bilinear_toy.py`, `cross_relations.py` and `explain_n16.py`;
`sweep.sh` forwards `AIM3_BACKEND`.  **Defaults are unchanged** — `backend=python` everywhere.

One extension was needed to reach n > 16 at all: `legacy/aim3_lab.CONWAY` stops at n = 16, so
`src/aim3/legacy_fields.py` copies the n = 17..24 moduli out of `src/aim3/gf2n.TOY_MODULI` after
verifying each is irreducible (Rabin's test) and that the shared entries agree.  `legacy/` itself is
untouched (CLAUDE.md rule 3).

Speed, measured: at E12 n=11 D=4 (37092 × 46938) Python takes 24.7 s and Sage 15.5 s; at E04 n=11
ℓ=3 D=4 (130812 × 149986) Sage takes 331 s, which is days in pure Python.  Below about
3·10⁴ rows the Python engine is faster (Sage start-up plus the export dominate), so it stays the
default.

---

## 2. **Gate 1**

> *Gate 1 (ROADMAP §7): if E04 finds cross-S-box degree falls → T2 becomes top priority; otherwise
> T2 and T3 run in parallel.*

### 2.1 Verdict

**The D=4 degree falls of the full AIM3 system are a direct sum over the S-boxes.  `cross_excess`
is exactly 0 at every non-saturated size measured — ℓ = 1 at n = 17, 19; ℓ = 2 at n = 13, 15, 17;
ℓ = 3 at n = 13 — on up to three seeds each, and not only for the whole system but at *every
subset* of the S-box blocks.  The single exception is ℓ = 3 at n = 11, the smallest size at which four S-boxes can be
run at all, and it does not survive to n = 13 (§2.4).**

Tag: **[VERIFIED n = 13, 15, 17, 19 for ℓ ≤ 2; VERIFIED n = 13 for ℓ = 3; the n = 11 ℓ = 3
exception is read as a small-parameter effect, [OPEN] as a mechanism]**.

Gate 1 therefore closes the way the Tier-A batch guessed it would, but now on evidence that meets
`COMPLEXITY_ACCOUNTING §4` — four values of n, three seeds, and the full subset lattice — instead of
a single n per observable: **there are no cross-S-box degree falls, T2 is not promoted, and T2 and
T3 continue in parallel.**

### 2.2 Why the batch had to climb to n = 17 and 19

`saturation = dim(span_{D=4} ∩ deg≤3) / C(nv,≤3)`.  When it approaches 1 the fall count measures
the size of the monomial basis, not the structure, and the batch rule admits a run into the verdict
only below 0.9.  The sizes Tier A could reach are all saturated, and so is ℓ = 2 at n = 11:

* ℓ = 1: n = 13 → 0.998, n = 15 → 0.978 (both discarded), n = 17 → 0.525, n = 19 → 0.470 (used)
* ℓ = 2: n = 11 → 0.9997 (discarded), n = 13 → 0.455, n = 15 → 0.398, n = 17 → 0.349 (used)
* ℓ = 3: n = 11 → 0.521, n = 13 → 0.342 (both below the threshold, see §2.4)

That is the whole reason this job needed a server: the first usable size for ℓ = 1 is n = 17 and
for ℓ = 2 it is n = 13, and both are out of reach in pure Python.

### 2.3 The evidence

Full tables in §A.2.  In summary, for every run with saturation < 0.9:

| ℓ | n | seeds | exps | falls_D4 | per-block falls | cross_excess |
|---:|---:|---:|---|---:|---|---:|
| 1 | 17 | 1,2,3 | 2;3 | 153 | 119 + 34 | **0** |
| 1 | 19 | 1,2,3 | 2;3 | 171 | 133 + 38 | **0** |
| 2 | 13 | 1,2,3 | 2;3;4 | 144 | 91 + 26 + 27 | **0** |
| 2 | 15 | 1,2,3 | 2;4;7 | 210 | 105 + 30 + 75 | **0** |
| 2 | 17 | 1 | 2;3;4 | 187 | 119 + 34 + 34 | **0** |
| 3 | 13 | 1 | 2;3;4;5 | 170 | 91 + 26 + 27 + 26 | **0** |
| 3 | 11 | 1,2 | 2;3;4;5 | 1826 | 77 + 22 + 22 + 66 | **1639** ← §2.4 |

Three independent strengthenings of that:

1. **The whole subset lattice, not just the extremes.**  `server/jobs/gate1_pairs.py` measures
   `falls_D4` for *every* subset of the S-box blocks.  At ℓ = 2, n = 13 every one of the three
   singletons, three pairs and the triple satisfies `falls(B) = Σ_{i∈B} falls({i})` exactly; at
   ℓ = 3, n = 13 the three-first-layer triple does too.  A direct-sum claim that survives the pair
   level is much stronger than one that only matches at the top.
2. **E12's first-layer system** (6n quadratics in 3n variables, never saturated) gives
   `cross_excess = 0` at n = 13 and n = 14, with `syz₃ = 16n` exactly.
3. **The random control has no falls at all.**  A same-shape random quadratic system gives
   `falls_D4 = 0` at ℓ = 2 n = 13 (saturation 0.47) and at ℓ = 3 n = 11 (saturation 0.42).  So the
   falls that S_quad does show are structural — but note the flip side, stated plainly: because the
   control has *no* falls, it cannot show a non-zero `cross_excess` either, and therefore it cannot
   calibrate the ℓ = 3, n = 11 exception below.  That control does not exist yet.

### 2.4 The ℓ = 3, n = 11 exception — and why it is read as a small-parameter effect

At ℓ = 3 (four S-boxes, nv = 44), n = 11, seeds 1 and 2, saturation 0.52, the same measurement
gives `falls_D4(all four) = 1826` against `Σ falls_D4(block i alone) = 77 + 22 + 22 + 66 = 187`,
i.e. `cross_excess = 1639`.  That is not a sign error: writing `W_B` for
`span_{D=4}(blocks in B) ∩ deg≤3`, the quantity is `dim W_{all} − Σ_i dim W_i`, and
`Σ_i dim W_i ≥ dim(Σ_i W_i)` always, so a positive value proves `W_{all} ⊋ Σ_i W_i`.

The subset lattice localises it (blocks 0, 1, 2 are the three *first-layer* S-boxes, which share the
input `pt`; block 3 is the last S-box):

| subset | falls_D4 | Σ singles | cross_excess |
|---|---:|---:|---:|
| every singleton | 77 / 22 / 22 / 66 | = | **0** |
| {0,2}, {0,3}, {1,2}, {1,3}, {2,3} | — | = | **0** |
| {0,1} | 110 | 99 | **11 = n** |
| {0,1,2} — the three first-layer S-boxes | 265 | 121 | **144** |
| {0,1,3} — two first-layer + the last | 176 | 165 | **11** (only what {0,1} already had) |
| all four | 1826 | 187 | **1639** |

So at n = 11 the excess is not pairwise and first appears at the three-first-layer subset.  That
looked like a mechanism — three S-boxes sharing `pt` is exactly what the AIM3-V shape has and the
AIM3-I/III shape does not.  **The n = 13 run refutes it:**

| exps | n | {0,1,2} falls_D4 | Σ singles | cross_excess |
|---|---:|---:|---:|---:|
| 2;3;4;5 | 11 | 265 | 121 | **144** |
| 3;4;7;8 (all generic / `{−3}`) | 11 | 144 | 67 | **77 = 7n** |
| 2;3;4;5 | 13 | 144 | 144 | **0** |

The same three-first-layer subset, with the same exponent pattern, gives an excess at n = 11 and
**none** at n = 13; and the whole four-block system gives 1639 at n = 11 and **0** at n = 13.  A
quantity that is large at one size and exactly zero at the next is not a structural coupling.

Two further facts point the same way.  n = 11 is the *smallest* size at which ℓ = 3 runs at all,
and at the same n the ℓ = 2 system is 99.97 % saturated — n = 11 is exactly where these
measurements stop being about structure.  And the excess is not stable under the exponent choice
(144 vs 77 = 7n for two exponent lists at the same n), which a genuine cross-S-box mechanism would
have to be, given that `falls_D4` itself depends on e only through the class rule.

**Honest position.**  The ℓ = 3, n = 11 row is reported as measured and is **not** folded into the
verdict; the batch reads it as a small-parameter effect and says so, while noting that the control
which would settle it — four S-boxes on *disjoint* variable blocks — has not been run (§8).  Anyone
re-running Gate 1 should start there and at ℓ = 3, n = 15, 17.

### 2.5 What this means for T1 and T2

* **T1** gets what it needed: for ℓ = 1, 2 and 3 the low-degree Hilbert function of `S_quad` is the
  direct sum of the per-S-box ones, over four values of n and with the full subset lattice checked,
  so Thm T1.3 can be stated — with the ℓ = 3, n = 11 row cited as the one measurement that does not
  fit and the reason it is not believed.
* **T2** gets no fuel from cross-S-box coupling: the mechanism that would have promoted it is
  absent.  §5 adds an independent negative from a completely different observable — the measured F4
  solving degree of `S_quad` equals that of a *random* system of the same shape at every size
  reached, so the structural falls do not translate into an easier Gröbner computation at toy scale.
* The **kernel model** remains T2's live anomaly (§4.2): it produces 2·#eqs elements of degree ≤ 2
  at D = 3 where the same-shape random bilinear control produces none.  That is a variable-choice
  effect, not a cross-S-box one, and it is what T2 should be built on.

## 3. The fall-count rule outside the range it was fitted on

`predicted_falls_D4` (`src/aim3/structures.py`) is a fit to 36 runs with n ≤ 16.  This batch ran
every admissible e at **n = 17, 18, 19, 20** with the M4RI backend.

**Mismatches: see §A.3 — the rule reproduces every completed run exactly.**  Merged with the
Tier-A table (58 runs, n = 11..20) `summarise.py` exits 0 and reports
*"the rule reproduces all measured fall counts exactly"*; the regenerated table is
`results/server/20260831_E03/E03_summary_with_n17_20.md` (produced by
`bash server/jobs/e03_summary_merged.sh`; the committed Tier-A `results/E03_summary.md` is **not**
overwritten — this batch may only add files).

Two side observations that also survive the extension:

* `syz₃ = 8n` in every run, without exception, now for n = 11..20.
* `outside_ideal` is `2n` for **e = 2** and `0` for every other class *including* e = 2⁻¹ — the
  asymmetry flagged [OPEN] in the Tier-A log is reproduced at n = 17 and n = 19.

**Suggested label** (the ledger is *not* edited by this batch, per §9): the rule stays **[CONJ]**,
with its verified range widened from *n ≤ 16* to **n ≤ 20, every admissible e**.

### D = 5

`server/jobs/e03_d5.py` measures `falls_D5 = dim(span_{D=5} ∩ deg≤4) − rank_{D=4}` at n = 14 and
n = 16, with a same-shape random control and the saturation of the degree-≤4 space.  Numbers in
§A.3.  The headline: D = 5 does produce a large number of new elements of degree ≤ 4 —

| n | e | class | falls_D4 | = rule? | falls_D5 | falls_D5 / n | saturation of the degree-≤4 space |
|---:|---:|---|---:|:--:|---:|---:|---:|
| 14 | 3 | generic | 28 | yes (2n) | 882 | 63 | 0.583 |
| 14 | 9 | {−3,(−3)⁻¹} | 29 | yes (2n+1) | 700 | 50 | 0.576 |
| 16 | 3 | generic | 32 | yes (2n) | 1136 | 71 | 0.517 |
| 14 | — | **random control** | 0 | — | 7966 | 569 | **1.000** |

Three things are worth stating separately.

* `falls_D4` matches the class rule at every one of these (n, e), which is an independent
  re-derivation of §3 on a different code path.
* D = 5 does produce a large number of new elements of degree ≤ 4 — tens of n, not a small multiple
  — so the first-fall-degree-4 statement of T1.2 is not the end of the story.  It is not constant
  in n/e either (50n, 63n, 71n).  Three points at two values of n is far short of what
  `COMPLEXITY_ACCOUNTING §4` needs, so the numbers are recorded and **no rule is proposed**.
* **The random control is useless at D = 5**: at n = 14 it is *fully saturated*
  (`dim(span_{D=5} ∩ deg≤4) = C(28,≤4)` exactly, saturation 1.000), so its 569n has no structural
  content and cannot serve as a baseline for the S-box's 63n.  A D = 5 comparison needs a control
  that is not pinned at the top — a larger n, or a control with fewer equations.  Until then the
  D = 5 fall counts are an observation, not evidence of structure.

---

## 4. E12 and E06 extensions

### 4.1 Cross-S-box relations at degree 4 — and a correction

`cross_relations.py 16 4` finds a **33-dimensional** space of Boolean quartics vanishing on the
joint graph {(t₀,t₁)} at n = 16 (65534 points, 41449 monomials — the smallest size at which the
question can be asked at all; at n = 14 the graph has 16382 points against 24158 monomials, so
*any* instance has ≥ 7776 relations from the point deficit alone and the measurement is
meaningless).

`explain_n16.py --deg 4` decomposes it as *17 supported on t₁ only, 0 on t₀ only*, which read
naively says **16 genuine cross-S-box quartics**.  That reading is wrong, and this is the batch's
main methodological correction:

> If `g = Tr(t₁^{−r₁})` is a cubic relation (it is: `wt(−r₁) = 3` at (n,e₁) = (16,5)), then every
> product `g·t₀ᵢ` is a quartic relation that touches **both** blocks and is not a cross relation at
> all — it is a multiple of a univariate one.  There are 16 variables t₀ᵢ, and 16 is exactly the
> number the naive test returns.

`server/jobs/e12_criteria_subtract.py` does what the E12 README prescribes instead: it builds each
criterion from its ANF (so its Boolean degree is *measured* — 3 and 5 at n = 16, 9 and 11 at
n = 14, matching `wt(−rᵢ)` exactly), spans the whole degree-≤4 part of the ideal they generate, and
subtracts it.  Result:

| n | deg | dim R (all relations) | dim C (ideal of the criteria) | genuine cross |
|---:|---:|---:|---:|---:|
| 16 | 3 | 1 | 1 | **0** |
| 16 | 4 | 33 | 33 | **0** |

**Claim K survives at degree 4**: every relation of degree ≤ 4 on (t₀,t₁) at n = 16 is a multiple
of a univariate image criterion.  This is a strengthening of §2.2-K of the ROADMAP, and it removes
the only candidate mechanism by which a cross-S-box *relation* could feed a cross-S-box *fall*.

### 4.2 E06 kernel model

The Tier-A observation §2.2-L now holds at **n = 11, 13 and 15** as it did at n = 7, 9: the kernel
model has `syz₃ = 0` and

| n | #eqs = 6n | `dim(span_{D=3} ∩ deg≤2)` | = 2·#eqs? | random bilinear control |
|---:|---:|---:|:--:|---:|
| 11 | 66 | 132 | yes | 66 = #eqs (no falls) |
| 13 | 78 | 156 | yes | 78 = #eqs (no falls) |
| 15 | 90 | 180 | yes | — |

i.e. 6n elements of degree ≤ 2 beyond the equations themselves, where the **same-shape random
bilinear control produces none**.  That is the one live anomaly T2 has, and it is a
variable-choice effect, not a cross-S-box one.  Tag: **[VERIFIED n = 7, 9, 11, 13, 15]** (was
n = 7, 9).

At D = 4 the kernel model was run at n = 11 only: `dim_low = 9075`, 4653 new elements, on a
145992 × 768692 matrix that takes M4RI 35 minutes.  n = 13 would be 240396 × 1505583 (≈ 45 GB, ~6 h)
and n = 15 240396 × ... ≈ 123 GB and about a day — out of proportion to this batch and, at n = 15,
above the one-Tier-C-job cap, so D = 4 there is listed as not done (§8).
Sizing, measured rather than estimated: the D = 4 Macaulay matrix of the kernel model is
145992 × 768692 at n = 11 (≈ 14 GB dense) and 240396 × 1505583 at n = 13 (≈ 45 GB); **n = 15 would
be ≈ 123 GB and about a day of M4RI and is not attempted** (§8).

---

## 5. Real solving degree `d_solv`, and block orders

This is the observable saturation cannot destroy, and the genuinely new measurement of the batch.
`export_systems.py` writes every instance in msolve and Magma syntax; `run_f4_degrees.py` runs
`msolve -v 2 -g 2 -t 8` and Magma `SetVerbose("Faugere",2)` and parses the step degrees.

**Cross-validation.**  At n = 7, ℓ = 2, `S_quad`, msolve and Magma report the same profile —
`d_solv = 3` over 5 steps.  Across the whole grid the two tools never disagree on `d_solv` (§A.5);
they differ only in step count and matrix size, as expected from different F4 strategies.

**`-g 2` matters**: without it msolve continues to a rational parametrisation and re-runs F4 about
twenty times, so neither the step count nor the wall clock would describe one F4 run.

`d_solv` vs n, grevlex — identical for msolve and Magma everywhere, and identical across
the two seeds (a cell shows every value observed):

| system | ℓ | n=7 | n=9 | n=11 | n=13 | n=15 |
|---|---:|---:|---:|---:|---:|---:|
| S_quad | 1 | 3 | 3 | 4 | 4 | 4 |
| P1only | 1 | 3 | 3 | 4 | 4 | 4 |
| random | 1 | 3 | 3 | 4 | 4 | 4 |
| S_quad | 2 | 3 | 4 | 4 | 4 | — |
| P1only | 2 | 3 | 4 | 4 | 4 | — |
| random | 2 | 4 | 4 | 4 | — | — |

**The finding, stated conservatively.**  `d_solv` takes the value 3 and then 4, and the step comes
one size earlier for three S-boxes than for two: at ℓ = 1 it is 3 at n = 7, 9 and 4 at n = 11, 13,
15; at ℓ = 2 it is 3 at n = 7 and 4 from n = 9 on.  And — this is the point — **it is the same for
`S_quad` and for a random Boolean quadratic system of the same shape at every size reached**, with
`random` occasionally *higher* (4 against 3 at n = 7, ℓ = 2).  The structural degree falls of T1
(2n per S-box) therefore do **not** buy a lower solving degree at these sizes; whatever the ideal
structure does to the Hilbert function, F4 does not see it as an easier computation.

Five values of n satisfy the ≥ 4 requirement of `COMPLEXITY_ACCOUNTING §4` for stating that
observation, but **not** for extrapolating a growth law: the data is a single step from 3 to 4 and
then a plateau, which is exactly the shape that must not be fitted.  No comparison with B1's
`d_reg = 19/25/44` is drawn here beyond noting that a plateau at 4 over n = 11..15 gives no
evidence either way at n = 128 — the toy systems are over-determined enough (9n equations in 3n
variables) that F4 reaches a single solution almost immediately, and that regime may or may not be
the one the real parameters are in.

The Gröbner bases themselves are essentially linear (`gb_max_degree` 1 or 2, `gb_size ≈ nv`): at
these toy sizes the systems have a single solution and F4 finds it at degree 4.

### Block orders (E10)

Two facts about Magma 2.24 shape what this experiment can and cannot say, and both were found the
hard way:

1. `BooleanPolynomialRing(n, order)` accepts an order **name** only; an `"elim"` order raises
   *"Bad argument types"*.  Block orders therefore run on `PolynomialRing(GF(2), n, "elim", S)` with
   the field equations added by hand — a different, slower engine.  Every row records which ring was
   used, and the only fair baseline for a block-order row is `genref` (grevlex on the *same* ring),
   not E04's Boolean-ring `grevlex`.
2. The declared elimination order **is** honoured in the returned basis (a direct check: the same
   toy system gives a 30-element basis under grevlex and a 15-element one under `elim`), but
   `SetVerbose("Faugere", 2)` reports `Order: Graded Reverse Lexicographical (bit vector word)` in
   both cases — **Magma's GF(2) F4 runs in grevlex and converts afterwards.**  A step-degree log
   under a block order is therefore a grevlex log, and `d_solv` is not comparable across orders with
   this tool.  The `d_solv` column of `results/E10_server_block_orders.csv` is grevlex's value *by
   construction*, not by measurement, and must not be read as an order comparison.

A third fact limits the experiment independently of the tool: these toy instances have a **single**
solution over F₂, and the reduced Gröbner basis of such an ideal is the same in every monomial
order.  Consistently, `gb_size` and `gb_max_degree` are identical across `genref`, `blockPT` and
`blockT` in every completed row, and so are the times to within noise.  **T2.2 gets no traction on
instances of this kind**; a real block-order comparison needs either a positive-dimensional
instance or a tool that computes directly in the block order (Sage/PolyBoRi does — §8).

---

## 6. E05 — resultant of the two-valued-inverse reduced system

The pipeline works end to end: `reduced_system_toy.py --export` writes the two sparse equations in
the format the E05 README specifies, `magma/reduced_system_resultant.m` rebuilds them over
`GF(2,n)`, **verifies that the true key is a root of both before doing anything expensive**, and
then times `Resultant(Eq1, Eq2, t1)`.  At n = 8 (exps 3,5,7, branch (ε₀,ε₁,ε₂) = (1,0,1)) the
system is

    Eq1:  64 terms, deg_{t0} = deg_{t1} = 254 = 2ⁿ−2
    Eq2: 678 terms, deg_{t0} = 254, deg_{t1} = 128
    Bézout bound  deg_{t0}·deg_{t1}:  64516 for Eq1 (2^{2n} = 65536), 32512 for Eq2

and the Sylvester matrix Magma has to take the determinant of is therefore about
`deg_{t1}(Eq1) + deg_{t1}(Eq2) = 382` square with entries of degree up to 254 in t₀.

**The measurement of this batch is that n = 8 does not finish.**  Magma ran `Resultant` for the
full **2-hour cap** without returning; the run is recorded with `censored = 1`, its numbers are
lower bounds, and the escalation stopped there (batch rule 0.2) — n = 9 and above were not
attempted.  `n = 7` was then launched to obtain at least one *completed* size, so that the jump has
a lower bound rather than nothing at all; see §A.7 for the rows as collected.

That is a weaker statement than the batch asked for — the T3 termination criterion needs the slope
of `log2(time)` over ≥ 4 sizes and there is at most one completed size — but it is a real and
useful one: **the direct resultant of the mixed-form reduced system is already out of reach at
n = 8**, which is 2⁴ = 16 times smaller than the smallest real parameter set in the exponent of the
Bézout bound.  The ROADMAP §3.3 expectation (T3.2, a negative result) is consistent with it; it is
not yet demonstrated by it.

---

## 7. Findings that contradict, or correct, what was expected

1. **A large `cross_excess` at ℓ = 3, n = 11 that vanishes at n = 13** (§2.4).  For most of this
   batch it looked like a shape effect — the excess localises onto the three first-layer S-boxes,
   exactly the structure AIM3-V has and AIM3-I/III does not — and the report was written that way
   before the n = 13 run came back with exactly 0.  It is now read as a small-parameter effect and
   is **not** in the verdict.  Recorded here because the intermediate reading was wrong and because
   the 0.9 saturation threshold evidently does not protect against everything at the smallest size
   a shape can be run at.
2. **`explain_n16.py --deg 4` over-reports cross relations by construction** (§4.1).  Its
   "supported on t₀ only / t₁ only" test is correct at degree 3 and wrong at degree 4.  Read
   naively it says 16 genuine cross-S-box quartics exist at n = 16; the correct count is 0.  The
   script was not changed beyond adding `--deg`/`--backend`; the corrected computation lives in
   `server/jobs/e12_criteria_subtract.py`.
3. **`cross_relations.py n 4` is meaningless at n = 14** — 16382 points against 24158 monomials
   forces ≥ 7776 relations for any instance whatsoever.  n = 16 is the first size at which the
   degree-4 question exists.
4. **A block-order F4 log cannot be obtained from Magma 2.24** (§5).  This invalidates the E10
   plan as written in the ROADMAP, not the E10 question.
5. **The per-block fall count in the big ring is not always E03's number.**  For generic and
   `{−3}`-class exponents the per-block value in (ℓ+1)n variables equals the 2n-variable value of
   E03 exactly; for the `{−2}` class it does not (e = 5 at n = 11 gives 6n in 44 variables against
   4n in 22; e = 7 at n = 15 gives 5n against 4n).  `cross_excess` is unaffected — it compares
   like with like — but any attempt to read the per-block column against `results/E03_summary.md`
   directly will trip over this.
6. **A timed-out Magma or msolve survives the driver.**  `subprocess.run(timeout=…)` kills
   `/usr/bin/time`, not its grandchild; one such orphan held 8 GB for half an hour before it was
   found. Both drivers now start the tool in its own process group and kill the group.
7. **msolve's default mode is not one F4 run** — it re-runs the whole computation about twenty
   times on the way to a rational parametrisation.  Any timing taken without `-g 2` is a timing of
   twenty runs.

---

## 8. Not done, and why

* **E06 at D = 4 for n = 13 and 15** — the kernel model's degree-4 Macaulay matrix is
  240396 × 1505583 (≈ 45 GB, ~6 h of M4RI) at n = 13 and ≈ 123 GB and about a day at n = 15, the
  latter above the one-Tier-C-job cap of `server/README.md`.  D = 4 is reported at n = 11 only;
  n = 13 and 15 were run to D = 3, which is where the `2·#eqs` claim lives.
* **A decoupled control for the ℓ = 3 anomaly** — four S-boxes on disjoint variable blocks, which
  is the only control that can prove the measurement itself is sound.  The random control cannot
  play that role because it has no falls at all.
* **A true block-order F4 profile** — needs Sage/PolyBoRi (which computes directly in a block
  order) or a positive-dimensional instance; Magma 2.24 cannot provide it (§5).
* **A uniform-GL_n matrix distribution** as a second control — every instance in this batch still
  uses the spec's unit-LU, as in Tier A.
* **ℓ = 3 at n ≥ 17, and ℓ = 2 at n ≥ 21** — the next non-saturated sizes; the D=4 matrices are
  4·10⁵ × 8·10⁵ and up.  A fourth value of n is what the ℓ = 3 finding needs before any growth law
  may be written down.
* **The E05 slope** — a single size (n = 8) that does not finish inside two hours is not a slope.
  The next batch should either (a) reduce the exponents modulo `2ⁿ−1` into a Kronecker
  substitution before calling `Resultant`, or (b) use `Resultant` on the *bilinear* form of Eq1
  (which separates as `W(t₀) + V(t₁)`) to eliminate t₁ by table lookup instead of by determinant —
  `count_roots` in `reduced_system_toy.py` already does exactly that in O(2ⁿ) and is the reason the
  Tier-A root counts exist at all.
* **The E10 block-order runs at n = 13 and 15 were stopped on purpose.**  Once it was established
  that Magma computes the GF(2) F4 phase in grevlex whatever order is declared (§5), those runs
  could not answer the question they were launched for, and they were consuming ~9 GB each while
  the Gate-1 jobs needed the memory.  The n = 7, 9, 11 rows are kept because they document the
  negative result.
* **No attack** is claimed anywhere in this batch.  Everything here is structure, modelling and
  cost measurement.

---

## 9. Suggested rows for `paper/claims_ledger.md` (not applied)

The ledger is not edited by this batch.  Proposed additions/changes, for a human to apply:

| claim | suggested status | support |
|---|---|---|
| `syz₃(S_e) = 8n` for a single S-box | [VERIFIED n = 11..20, every admissible e] (was n ≤ 16) | `results/E03_server_n17_20.csv` |
| first fall degree 4 with `2n + Δ(e)` new cubics, Δ given by the class rule | [CONJ], verified range widened to **n ≤ 20** | `results/server/20260831_E03/E03_summary_with_n17_20.md` |
| `outside_ideal = 2n` iff e = 2 (not e = 2⁻¹) | [OPEN], reproduced at n = 17, 19 | same |
| T1.3 direct sum of the D=4 falls | [VERIFIED ℓ = 1 at n = 17, 19; ℓ = 2 at n = 13, 15, 17; ℓ = 3 at n = 13; seeds 1–3; all subsets of the blocks] (was n = 11, 14) | `results/E04_server_gate1.csv`, `results/E04_server_gate1_pairs.csv` |
| the ℓ = 3, n = 11 excess | **[OPEN]** — record it, do not build on it; it is 0 at n = 13 and not stable under the exponent choice | `results/E04_server_gate1_pairs.csv` |
| E06 §2.2-L: kernel model has `dim(span_{D=3} ∩ deg≤2) = 2·#eqs`, the random bilinear control has none | [VERIFIED n = 7, 9, 11, 13, 15] (was n = 7, 9) | `results/E06_server_kernel.csv` |
| claim K, no cross-S-box relations | [VERIFIED n = 16 at degree ≤ **4**] (was degree 3) | `results/E12_server_criteria_subtract.csv` |
| measured `d_solv` of `S_quad` equals that of a same-shape random system | [VERIFIED n = 7, 9, 11, 13, 15] — new row | `results/E04_f4_degrees.csv` |
| block orders give no advantage on these instances | [VERIFIED, with the tool caveat of §5] — new row, phrased as a limitation of the measurement | `results/E10_server_block_orders.csv` |

---

<!--TABLES-->

## Appendix — generated tables

### A.1 Backend alignment

`server/jobs/verify_backend.py`: **76/76 comparisons agreed** → **BACKEND ALIGNED**.

| check | quantity | python | sage | agree |
|---|---|---:|---:|:--:|
| E03 n=11 e=3 | `D=3 rows` | 759 | 759 | yes |
| E03 n=11 e=3 | `D=3 cols` | 1794 | 1794 | yes |
| E03 n=11 e=3 | `D=3 rank` | 671 | 671 | yes |
| E03 n=11 e=3 | `D=3 syz` | 88 | 88 | yes |
| E03 n=11 e=3 | `D=3 dim_low` | 33 | 33 | yes |
| E03 n=11 e=3 | `D=4 rows` | 8382 | 8382 | yes |
| E03 n=11 e=3 | `D=4 cols` | 9109 | 9109 | yes |
| E03 n=11 e=3 | `D=4 rank` | 6006 | 6006 | yes |
| E03 n=11 e=3 | `D=4 syz` | 2376 | 2376 | yes |
| E03 n=11 e=3 | `D=4 dim_low` | 693 | 693 | yes |
| E03 n=11 e=3 | `falls_D4` | 22 | 22 | yes |
| E03 n=11 e=3 | `dim I(V)_<=3` | 693 | 693 | yes |
| E03 n=11 e=3 | `point-matrix rank` | 1101 | 1101 | yes |
| E03 n=13 e=7 | `D=3 rows` | 1053 | 1053 | yes |
| E03 n=13 e=7 | `D=3 cols` | 2952 | 2952 | yes |
| E03 n=13 e=7 | `D=3 rank` | 949 | 949 | yes |
| E03 n=13 e=7 | `D=3 syz` | 104 | 104 | yes |
| E03 n=13 e=7 | `D=3 dim_low` | 39 | 39 | yes |
| E03 n=13 e=7 | `D=4 rows` | 13728 | 13728 | yes |
| E03 n=13 e=7 | `D=4 cols` | 17902 | 17902 | yes |
| E03 n=13 e=7 | `D=4 rank` | 10348 | 10348 | yes |
| E03 n=13 e=7 | `D=4 syz` | 3380 | 3380 | yes |
| E03 n=13 e=7 | `D=4 dim_low` | 1040 | 1040 | yes |
| E03 n=13 e=7 | `falls_D4` | 91 | 91 | yes |
| E03 n=13 e=7 | `dim I(V)_<=3` | 1040 | 1040 | yes |
| E03 n=13 e=7 | `point-matrix rank` | 1912 | 1912 | yes |
| E12 n=11 | `D=3 rows` | 2244 | 2244 | yes |
| E12 n=11 | `D=3 cols` | 6018 | 6018 | yes |
| E12 n=11 | `D=3 rank` | 2068 | 2068 | yes |
| E12 n=11 | `D=3 syz` | 176 | 176 | yes |
| E12 n=11 | `D=3 dim_low` | 66 | 66 | yes |
| E12 n=11 | `D=4 rows` | 37092 | 37092 | yes |
| E12 n=11 | `D=4 cols` | 46938 | 46938 | yes |
| E12 n=11 | `D=4 rank` | 29315 | 29315 | yes |
| E12 n=11 | `D=4 syz` | 7777 | 7777 | yes |
| E12 n=11 | `D=4 dim_low` | 2112 | 2112 | yes |
| E12 n=11 | `falls_D4` | 44 | 44 | yes |
| E12 n=11 | `falls_D4 S-box0 alone` | 22 | 22 | yes |
| E12 n=11 | `falls_D4 S-box1 alone` | 22 | 22 | yes |
| E12 n=11 | `cross_excess` | 0 | 0 | yes |
| E06 n=7 | `D=2 rows` | 42 | 42 | yes |
| E06 n=7 | `D=2 cols` | 904 | 904 | yes |
| E06 n=7 | `D=2 rank` | 42 | 42 | yes |
| E06 n=7 | `D=2 syz` | 0 | 0 | yes |
| E06 n=7 | `D=2 dim_low` | 0 | 0 | yes |
| E06 n=7 | `D=3 rows` | 1806 | 1806 | yes |
| E06 n=7 | `D=3 cols` | 12384 | 12384 | yes |
| E06 n=7 | `D=3 rank` | 1806 | 1806 | yes |
| E06 n=7 | `D=3 syz` | 0 | 0 | yes |
| E06 n=7 | `D=3 dim_low` | 84 | 84 | yes |
| E06 n=7 | `D=4 rows` | 37968 | 37968 | yes |
| E06 n=7 | `D=4 cols` | 124314 | 124314 | yes |
| E06 n=7 | `D=4 rank` | 36994 | 36994 | yes |
| E06 n=7 | `D=4 syz` | 974 | 974 | yes |
| E06 n=7 | `D=4 dim_low` | 3730 | 3730 | yes |
| E03 csv n=11 e=3 | `syz3` | 88 | 88 | yes |
| E03 csv n=11 e=3 | `falls_D4` | 22 | 22 | yes |
| E03 csv n=11 e=3 | `mac3_rank` | 671 | 671 | yes |
| E03 csv n=11 e=3 | `IV_deg3` | 693 | 693 | yes |
| E03 csv n=13 e=7 | `syz3` | 104 | 104 | yes |
| E03 csv n=13 e=7 | `falls_D4` | 91 | 91 | yes |
| E03 csv n=13 e=7 | `mac3_rank` | 949 | 949 | yes |
| E03 csv n=13 e=7 | `IV_deg3` | 1040 | 1040 | yes |
| E12 csv n=11 | `mac3_rank` | 2068 | 2068 | yes |
| E12 csv n=11 | `mac3_syz` | 176 | 176 | yes |
| E12 csv n=11 | `falls_D4` | 44 | 44 | yes |
| E12 csv n=11 | `cross_excess` | 0 | 0 | yes |
| E06 csv n=7 D=2 | `rank` | 42 | 42 | yes |
| E06 csv n=7 D=2 | `syz` | 0 | 0 | yes |
| E06 csv n=7 D=2 | `dim_low` | 0 | 0 | yes |
| E06 csv n=7 D=3 | `rank` | 1806 | 1806 | yes |
| E06 csv n=7 D=3 | `syz` | 0 | 0 | yes |
| E06 csv n=7 D=3 | `dim_low` | 84 | 84 | yes |
| E06 csv n=7 D=4 | `rank` | 36994 | 36994 | yes |
| E06 csv n=7 D=4 | `syz` | 974 | 974 | yes |
| E06 csv n=7 D=4 | `dim_low` | 3730 | 3730 | yes |

### A.2 Gate 1 — is the D=4 degree fall a direct sum over the S-boxes?

#### A.2.1 `full_system_falls.py --per-block` (S_quad)

`saturation = dim(span_{D=4} ∩ deg≤3) / C(nv,≤3)`.  The batch rule admits a run into the verdict only when saturation < 0.9; the other rows are listed but marked.

| ℓ | n | seed | exps | nv | falls_D4 | per-block falls | Σ | cross_excess | saturation | in verdict |
|---:|---:|---:|---|---:|---:|---|---:|---:|---:|:--:|
| 1 | 13 | 1 | 2;3 | 26 | 1049 | 91 + 26 | 117 | **932** | 0.9983 | no (saturated) |
| 1 | 13 | 2 | 2;3 | 26 | 1051 | 91 + 26 | 117 | **934** | 0.9990 | no (saturated) |
| 1 | 13 | 3 | 2;3 | 26 | 1053 | 91 + 26 | 117 | **936** | 0.9997 | no (saturated) |
| 1 | 15 | 1 | 2;4 | 30 | 1875 | 105 + 30 | 135 | **1740** | 0.9777 | no (saturated) |
| 1 | 15 | 2 | 2;4 | 30 | 1875 | 105 + 30 | 135 | **1740** | 0.9777 | no (saturated) |
| 1 | 15 | 3 | 2;4 | 30 | 1875 | 105 + 30 | 135 | **1740** | 0.9777 | no (saturated) |
| 1 | 17 | 1 | 2;3 | 34 | 153 | 119 + 34 | 153 | **0** | 0.5245 | yes |
| 1 | 17 | 2 | 2;3 | 34 | 153 | 119 + 34 | 153 | **0** | 0.5245 | yes |
| 1 | 17 | 3 | 2;3 | 34 | 153 | 119 + 34 | 153 | **0** | 0.5245 | yes |
| 1 | 19 | 1 | 2;3 | 38 | 171 | 133 + 38 | 171 | **0** | 0.4699 | yes |
| 1 | 19 | 2 | 2;3 | 38 | 171 | 133 + 38 | 171 | **0** | 0.4699 | yes |
| 1 | 19 | 3 | 2;3 | 38 | 171 | 133 + 38 | 171 | **0** | 0.4699 | yes |
| 2 | 11 | 1 | 2;3;4 | 33 | 2914 | 77 + 22 + 22 | 121 | **2793** | 0.9997 | no (saturated) |
| 2 | 11 | 2 | 2;3;4 | 33 | 2913 | 77 + 22 + 22 | 121 | **2792** | 0.9995 | no (saturated) |
| 2 | 11 | 3 | 2;3;4 | 33 | 2914 | 77 + 22 + 22 | 121 | **2793** | 0.9997 | no (saturated) |
| 2 | 13 | 1 | 2;3;4 | 39 | 144 | 91 + 26 + 27 | 144 | **0** | 0.4548 | yes |
| 2 | 13 | 2 | 2;3;4 | 39 | 144 | 91 + 26 + 27 | 144 | **0** | 0.4548 | yes |
| 2 | 13 | 3 | 2;3;4 | 39 | 144 | 91 + 26 + 27 | 144 | **0** | 0.4548 | yes |
| 2 | 15 | 1 | 2;4;7 | 45 | 210 | 105 + 30 + 75 | 210 | **0** | 0.3980 | yes |
| 2 | 15 | 2 | 2;4;7 | 45 | 210 | 105 + 30 + 75 | 210 | **0** | 0.3980 | yes |
| 2 | 15 | 3 | 2;4;7 | 45 | 210 | 105 + 30 + 75 | 210 | **0** | 0.3980 | yes |
| 2 | 17 | 1 | 2;3;4 | 51 | 187 | 119 + 34 + 34 | 187 | **0** | 0.3492 | yes |
| 3 | 11 | 1 | 2;3;4;5 | 44 | 1826 | 77 + 22 + 22 + 66 | 187 | **1639** | 0.5208 | yes |
| 3 | 11 | 2 | 2;3;4;5 | 44 | 1826 | 77 + 22 + 22 + 66 | 187 | **1639** | 0.5208 | yes |
| 3 | 11 | 3 | 2;3;4;5 | 44 | 1826 | 77 + 22 + 22 + 66 | 187 | **1639** | 0.5208 | yes |
| 3 | 13 | 1 | 2;3;4;5 | 52 | 170 | 91 + 26 + 27 + 26 | 170 | **0** | 0.3417 | yes |

Controls on the same instances (no per-block decomposition is defined for them):

| ℓ | n | seed | system | falls_D4 | saturation |
|---:|---:|---:|---|---:|---:|
| 1 | 13 | 1 | P1only | 1196 | 0.643 |
| 1 | 13 | 1 | random | 845 | 0.9997 |
| 1 | 13 | 2 | P1only | 1196 | 0.643 |
| 1 | 13 | 2 | random | 845 | 0.9997 |
| 1 | 13 | 3 | P1only | 1196 | 0.643 |
| 1 | 13 | 3 | random | 845 | 0.9997 |
| 1 | 15 | 1 | P1only | 1620 | 0.5634 |
| 1 | 15 | 1 | random | 1735 | 0.9998 |
| 1 | 15 | 2 | P1only | 1620 | 0.5634 |
| 1 | 15 | 2 | random | 1735 | 0.9998 |
| 1 | 15 | 3 | P1only | 1620 | 0.5634 |
| 1 | 15 | 3 | random | 1735 | 0.9998 |
| 1 | 17 | 1 | P1only | 2108 | 0.5012 |
| 1 | 17 | 1 | random | 3009 | 0.9998 |
| 1 | 17 | 2 | P1only | 2108 | 0.5012 |
| 1 | 17 | 2 | random | 3009 | 0.9998 |
| 1 | 17 | 3 | P1only | 2108 | 0.5012 |
| 1 | 17 | 3 | random | 3009 | 0.9998 |
| 1 | 19 | 1 | P1only | 2660 | 0.4513 |
| 1 | 19 | 1 | random | 0 | 0.4844 |
| 1 | 19 | 2 | P1only | 2660 | 0.4513 |
| 1 | 19 | 2 | random | 0 | 0.4844 |
| 1 | 19 | 3 | P1only | 2660 | 0.4513 |
| 1 | 19 | 3 | random | 0 | 0.4844 |
| 2 | 11 | 1 | P1only | 1980 | 0.5155 |
| 2 | 11 | 1 | random | 2651 | 0.9998 |
| 2 | 11 | 2 | P1only | 1980 | 0.5155 |
| 2 | 11 | 2 | random | 2651 | 0.9998 |
| 2 | 11 | 3 | P1only | 1980 | 0.5155 |
| 2 | 11 | 3 | random | 2651 | 0.9998 |
| 2 | 13 | 1 | P1only | 2808 | 0.4403 |
| 2 | 13 | 1 | random | 0 | 0.4718 |
| 2 | 13 | 2 | P1only | 2808 | 0.4403 |
| 2 | 13 | 2 | random | 0 | 0.4718 |
| 2 | 13 | 3 | P1only | 2808 | 0.4403 |
| 2 | 13 | 3 | random | 0 | 0.4718 |
| 2 | 15 | 1 | P1only | 3780 | 0.3842 |
| 2 | 15 | 1 | random | 0 | 0.4079 |
| 2 | 15 | 2 | P1only | 3780 | 0.3842 |
| 2 | 15 | 2 | random | 0 | 0.4079 |
| 2 | 15 | 3 | P1only | 3780 | 0.3842 |
| 2 | 17 | 1 | P1only | 4896 | 0.3407 |
| 3 | 11 | 1 | P1only | 3608 | 0.3926 |
| 3 | 11 | 1 | random | 0 | 0.4173 |
| 3 | 11 | 2 | P1only | 3608 | 0.3926 |
| 3 | 11 | 2 | random | 0 | 0.4173 |
| 3 | 11 | 3 | P1only | 3608 | 0.3926 |
| 3 | 11 | 3 | random | 0 | 0.4173 |
| 3 | 13 | 1 | P1only | 5096 | 0.3344 |

#### A.2.2 The subset lattice, with the random control (`server/jobs/gate1_pairs.py`)

`cross_excess(B) = falls(B) − Σ_{i∈B} falls({i})` for every subset B of the S-box blocks, in the same nv variables.  The `random` rows are the *same shape* cut into blocks at the same positions: a non-zero `cross_excess` says something about AIM3 only if the control's is zero.

| ℓ | n | exps | system | blocks | eqs | D3 syz | falls_D4 | Σ singles | cross_excess | sat |
|---:|---:|---|---|---|---:|---:|---:|---:|---:|---:|
| 2 | 13 | 2;3;4 | S_quad | {0} | 39 | 104 | 91 | 91 | 0 | 0.1559 |
| 2 | 13 | 2;3;4 | S_quad | {1} | 39 | 104 | 26 | 26 | 0 | 0.1494 |
| 2 | 13 | 2;3;4 | S_quad | {2} | 39 | 104 | 27 | 27 | 0 | 0.1495 |
| 2 | 13 | 2;3;4 | S_quad | {0+1} | 78 | 208 | 117 | 117 | 0 | 0.3053 |
| 2 | 13 | 2;3;4 | S_quad | {0+2} | 78 | 208 | 118 | 118 | 0 | 0.3054 |
| 2 | 13 | 2;3;4 | S_quad | {1+2} | 78 | 208 | 53 | 53 | 0 | 0.2989 |
| 2 | 13 | 2;3;4 | S_quad | {0+1+2} | 117 | 312 | 144 | 144 | 0 | 0.4548 |
| 2 | 13 | 2;3;4 | random | {0} | 39 | 0 | 0 | 0 | 0 | 0.1573 |
| 2 | 13 | 2;3;4 | random | {1} | 39 | 0 | 0 | 0 | 0 | 0.1573 |
| 2 | 13 | 2;3;4 | random | {2} | 39 | 0 | 0 | 0 | 0 | 0.1573 |
| 2 | 13 | 2;3;4 | random | {0+1} | 78 | 0 | 0 | 0 | 0 | 0.3145 |
| 2 | 13 | 2;3;4 | random | {0+2} | 78 | 0 | 0 | 0 | 0 | 0.3145 |
| 2 | 13 | 2;3;4 | random | {1+2} | 78 | 0 | 0 | 0 | 0 | 0.3145 |
| 2 | 13 | 2;3;4 | random | {0+1+2} | 117 | 0 | 0 | 0 | 0 | 0.4718 |
| 2 | 15 | 2;4;7 | S_quad | {0} | 45 | 120 | 105 | 105 | 0 | 0.135 |
| 2 | 15 | 2;4;7 | S_quad | {1} | 45 | 120 | 30 | 30 | 0 | 0.13 |
| 2 | 15 | 2;4;7 | S_quad | {2} | 45 | 120 | 75 | 75 | 0 | 0.133 |
| 2 | 15 | 2;4;7 | S_quad | {0+1} | 90 | 240 | 135 | 135 | 0 | 0.265 |
| 2 | 15 | 2;4;7 | S_quad | {0+2} | 90 | 240 | 180 | 180 | 0 | 0.268 |
| 2 | 15 | 2;4;7 | S_quad | {1+2} | 90 | 240 | 105 | 105 | 0 | 0.263 |
| 2 | 15 | 2;4;7 | S_quad | {0+1+2} | 135 | 360 | 210 | 210 | 0 | 0.398 |
| 2 | 15 | 2;4;7 | random | {0} | 45 | 0 | 0 | 0 | 0 | 0.136 |
| 2 | 15 | 2;4;7 | random | {1} | 45 | 0 | 0 | 0 | 0 | 0.136 |
| 2 | 15 | 2;4;7 | random | {2} | 45 | 0 | 0 | 0 | 0 | 0.136 |
| 2 | 15 | 2;4;7 | random | {0+1} | 90 | 0 | 0 | 0 | 0 | 0.2719 |
| 2 | 15 | 2;4;7 | random | {0+2} | 90 | 0 | 0 | 0 | 0 | 0.2719 |
| 2 | 15 | 2;4;7 | random | {1+2} | 90 | 0 | 0 | 0 | 0 | 0.2719 |
| 3 | 11 | 2;3;4;5 | S_quad | {0} | 33 | 88 | 77 | 77 | 0 | 0.1035 |
| 3 | 11 | 2;3;4;5 | S_quad | {1} | 33 | 88 | 22 | 22 | 0 | 0.0997 |
| 3 | 11 | 2;3;4;5 | S_quad | {2} | 33 | 88 | 22 | 22 | 0 | 0.0997 |
| 3 | 11 | 2;3;4;5 | S_quad | {3} | 33 | 88 | 66 | 66 | 0 | 0.1028 |
| 3 | 11 | 2;3;4;5 | S_quad | {0+1} | 66 | 176 | 110 | 99 | **11** | 0.204 |
| 3 | 11 | 2;3;4;5 | S_quad | {0+2} | 66 | 176 | 99 | 99 | 0 | 0.2032 |
| 3 | 11 | 2;3;4;5 | S_quad | {0+3} | 66 | 176 | 143 | 143 | 0 | 0.2063 |
| 3 | 11 | 2;3;4;5 | S_quad | {1+2} | 66 | 176 | 44 | 44 | 0 | 0.1994 |
| 3 | 11 | 2;3;4;5 | S_quad | {1+3} | 66 | 176 | 88 | 88 | 0 | 0.2025 |
| 3 | 11 | 2;3;4;5 | S_quad | {2+3} | 66 | 176 | 88 | 88 | 0 | 0.2025 |
| 3 | 11 | 2;3;4;5 | S_quad | {0+1+2} | 99 | 264 | 265 | 121 | **144** | 0.313 |
| 3 | 11 | 2;3;4;5 | S_quad | {0+1+3} | 99 | 264 | 176 | 165 | **11** | 0.3068 |
| 3 | 11 | 2;3;4;5 | S_quad | {0+2+3} | 99 | 264 | 165 | 165 | 0 | 0.306 |
| 3 | 11 | 2;3;4;5 | S_quad | {1+2+3} | 99 | 264 | 110 | 110 | 0 | 0.3021 |
| 3 | 11 | 2;3;4;5 | S_quad | {0+1+2+3} | 132 | 352 | 1826 | 187 | **1639** | 0.5208 |
| 3 | 11 | 2;3;4;5 | random | {0} | 33 | 0 | 0 | 0 | 0 | 0.1043 |
| 3 | 11 | 2;3;4;5 | random | {1} | 33 | 0 | 0 | 0 | 0 | 0.1043 |
| 3 | 11 | 2;3;4;5 | random | {2} | 33 | 0 | 0 | 0 | 0 | 0.1043 |
| 3 | 11 | 2;3;4;5 | random | {3} | 33 | 0 | 0 | 0 | 0 | 0.1043 |
| 3 | 11 | 2;3;4;5 | random | {0+1} | 66 | 0 | 0 | 0 | 0 | 0.2086 |
| 3 | 11 | 2;3;4;5 | random | {0+2} | 66 | 0 | 0 | 0 | 0 | 0.2086 |
| 3 | 11 | 2;3;4;5 | random | {0+3} | 66 | 0 | 0 | 0 | 0 | 0.2086 |
| 3 | 11 | 2;3;4;5 | random | {1+2} | 66 | 0 | 0 | 0 | 0 | 0.2086 |
| 3 | 11 | 2;3;4;5 | random | {1+3} | 66 | 0 | 0 | 0 | 0 | 0.2086 |
| 3 | 11 | 2;3;4;5 | random | {2+3} | 66 | 0 | 0 | 0 | 0 | 0.2086 |
| 3 | 11 | 3;4;7;8 | S_quad | {0} | 33 | 88 | 22 | 22 | 0 | 0.0997 |
| 3 | 11 | 3;4;7;8 | S_quad | {0} | 33 | 88 | 22 | 22 | 0 | 0.0997 |
| 3 | 11 | 3;4;7;8 | S_quad | {1} | 33 | 88 | 22 | 22 | 0 | 0.0997 |
| 3 | 11 | 3;4;7;8 | S_quad | {1} | 33 | 88 | 22 | 22 | 0 | 0.0997 |
| 3 | 11 | 3;4;7;8 | S_quad | {2} | 33 | 88 | 23 | 23 | 0 | 0.0998 |
| 3 | 11 | 3;4;7;8 | S_quad | {2} | 33 | 88 | 23 | 23 | 0 | 0.0998 |
| 3 | 11 | 3;4;7;8 | S_quad | {3} | 33 | 88 | 23 | 23 | 0 | 0.0998 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+1} | 66 | 176 | 44 | 44 | 0 | 0.1994 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+2} | 66 | 176 | 45 | 45 | 0 | 0.1994 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+3} | 66 | 176 | 45 | 45 | 0 | 0.1994 |
| 3 | 11 | 3;4;7;8 | S_quad | {1+2} | 66 | 176 | 45 | 45 | 0 | 0.1994 |
| 3 | 11 | 3;4;7;8 | S_quad | {1+3} | 66 | 176 | 45 | 45 | 0 | 0.1994 |
| 3 | 11 | 3;4;7;8 | S_quad | {2+3} | 66 | 176 | 46 | 46 | 0 | 0.1995 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+1+2} | 99 | 264 | 144 | 67 | **77** | 0.3045 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+1+2} | 99 | 264 | 144 | 67 | **77** | 0.3045 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+1+3} | 99 | 264 | 67 | 67 | 0 | 0.2991 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+2+3} | 99 | 264 | 68 | 68 | 0 | 0.2992 |
| 3 | 11 | 3;4;7;8 | S_quad | {1+2+3} | 99 | 264 | 68 | 68 | 0 | 0.2992 |
| 3 | 11 | 3;4;7;8 | S_quad | {0+1+2+3} | 132 | 352 | 1859 | 90 | **1769** | 0.5231 |
| 3 | 13 | 2;3;4;5 | S_quad | {0} | 39 | 104 | 91 | 91 | 0 | 0.0875 |
| 3 | 13 | 2;3;4;5 | S_quad | {0} | 39 | 104 | 91 | 91 | 0 | 0.0875 |
| 3 | 13 | 2;3;4;5 | S_quad | {1} | 39 | 104 | 26 | 26 | 0 | 0.0847 |
| 3 | 13 | 2;3;4;5 | S_quad | {1} | 39 | 104 | 26 | 26 | 0 | 0.0847 |
| 3 | 13 | 2;3;4;5 | S_quad | {2} | 39 | 104 | 27 | 27 | 0 | 0.0848 |
| 3 | 13 | 2;3;4;5 | S_quad | {2} | 39 | 104 | 27 | 27 | 0 | 0.0848 |
| 3 | 13 | 2;3;4;5 | S_quad | {3} | 39 | 104 | 26 | 26 | 0 | 0.0847 |
| 3 | 13 | 2;3;4;5 | S_quad | {0+1} | 78 | 208 | 117 | 117 | 0 | 0.1722 |
| 3 | 13 | 2;3;4;5 | S_quad | {0+2} | 78 | 208 | 118 | 118 | 0 | 0.1722 |
| 3 | 13 | 2;3;4;5 | S_quad | {0+1+2} | 117 | 312 | 144 | 144 | 0 | 0.257 |
| 3 | 15 | 2;4;7;8 | S_quad | {0} | 45 | 120 | 105 | 105 | 0 | 0.0757 |
| 3 | 15 | 2;4;7;8 | S_quad | {1} | 45 | 120 | 30 | 30 | 0 | 0.0736 |
| 3 | 15 | 2;4;7;8 | S_quad | {2} | 45 | 120 | 90 | 90 | 0 | 0.0753 |

#### A.2.3 E12 first-layer system (two S-boxes in 3n variables)

| n | exps | nv | mac3 rank | syz (16n) | falls_D4 | S-box 0 | S-box 1 | cross_excess |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 13 | 3;4 | 39 | 2912 | 208 (208) | 53 | 26 | 27 | **0** |
| 14 | 3;5 | 42 | 3388 | 224 (224) | 56 | 28 | 28 | **0** |

### A.3 The calibrated fall-count rule at n = 17..20 (outside the range it was fitted on)

`predicted_falls_D4` (`src/aim3/structures.py`) was calibrated on 36 runs with n ≤ 16.  Every admissible e (gcd(e,n)=1, 2 ≤ e ≤ n−2) is run here; `censored=1` rows carry lower bounds and are excluded from the mismatch count.

| n | e | m = e⁻¹ | syz₃ | = 8n | falls_D4 | falls/n | predicted | class | match | outside ideal | dim I(V)≤3 | s |
|---:|---:|---:|---:|:--:|---:|---:|---:|---|:--:|---:|---:|---:|
| 17 | 2 | 9 | 136 | yes | 119 | 7.00 | 119 | `{2,2^-1}` | yes | 34 | 1802 | 110.961 |
| 17 | 3 | 6 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 122.674 |
| 17 | 4 | 13 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 118.767 |
| 17 | 5 | 7 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 120.472 |
| 17 | 6 | 3 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 117.392 |
| 17 | 7 | 5 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 112.236 |
| 17 | 8 | 15 | 136 | yes | 68 | 4.00 | 68 | `{-2,-2^-1}` | yes | 0 | 1717 | 125.046 |
| 17 | 9 | 2 | 136 | yes | 119 | 7.00 | 119 | `{2,2^-1}` | yes | 0 | 1768 | 119.397 |
| 17 | 10 | 12 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 123.92 |
| 17 | 11 | 14 | 136 | yes | 35 | 2.06 | 35 | `{-3,-3^-1}` | yes | 0 | 1684 | 145.813 |
| 17 | 12 | 10 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 165.805 |
| 17 | 13 | 4 | 136 | yes | 34 | 2.00 | 34 | `generic` | yes | 0 | 1683 | 196.911 |
| 17 | 14 | 11 | 136 | yes | 35 | 2.06 | 35 | `{-3,-3^-1}` | yes | 0 | 1684 | 195.333 |
| 17 | 15 | 8 | 136 | yes | 68 | 4.00 | 68 | `{-2,-2^-1}` | yes | 0 | 1717 | 174.263 |
| 18 | 5 | 11 | 144 | yes | 36 | 2.00 | 36 | `generic` | yes | 0 | 1890 | 256.673 |
| 18 | 7 | 13 | 144 | yes | 36 | 2.00 | 36 | `generic` | yes | 0 | 1890 | 265.297 |
| 18 | 11 | 5 | 144 | yes | 36 | 2.00 | 36 | `generic` | yes | 0 | 1890 | 274.37 |
| 18 | 13 | 7 | 144 | yes | 36 | 2.00 | 36 | `generic` | yes | 0 | 1890 | 322.98 |
| 19 | 2 | 10 | 152 | yes | 133 | 7.00 | 133 | `{2,2^-1}` | yes | 38 | 2242 | 589.884 |
| 19 | 3 | 13 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes | 0 | 2109 | 684.517 |
| 19 | 4 | 5 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes | 0 | 2109 | 895.212 |
| 19 | 5 | 4 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes | 0 | 2109 | 784.771 |
| 19 | 6 | 16 | 152 | yes | 39 | 2.05 | 39 | `{-3,-3^-1}` | yes | 0 | 2110 | 880.317 |
| 19 | 7 | 11 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes |  |  | 81.661 |
| 19 | 8 | 12 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes |  |  | 68.501 |
| 19 | 9 | 17 | 152 | yes | 76 | 4.00 | 76 | `{-2,-2^-1}` | yes |  |  | 68.483 |
| 19 | 10 | 2 | 152 | yes | 133 | 7.00 | 133 | `{2,2^-1}` | yes |  |  | 83.14 |
| 19 | 11 | 7 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes |  |  | 119.632 |
| 19 | 12 | 8 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes |  |  | 84.079 |
| 19 | 13 | 3 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes |  |  | 73.072 |
| 19 | 14 | 15 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes |  |  | 76.975 |
| 19 | 15 | 14 | 152 | yes | 38 | 2.00 | 38 | `generic` | yes |  |  | 77.973 |
| 19 | 16 | 6 | 152 | yes | 39 | 2.05 | 39 | `{-3,-3^-1}` | yes |  |  | 72.539 |
| 19 | 17 | 9 | 152 | yes | 76 | 4.00 | 76 | `{-2,-2^-1}` | yes |  |  | 103.57 |
| 20 | 3 | 7 | 160 | yes | 40 | 2.00 | 40 | `generic` | yes | 0 | 2340 | 1581.416 |
| 20 | 7 | 3 | 160 | yes | 40 | 2.00 | 40 | `generic` | yes | 0 | 2340 | 1798.025 |
| 20 | 9 | 9 | 160 | yes | 41 | 2.05 | 41 | `{n/2-1}` | yes | 0 | 2341 | 1897.927 |
| 20 | 11 | 11 | 160 | yes | 170 | 8.50 | 170 | `{1+n/2}` | yes |  |  | 133.793 |
| 20 | 13 | 17 | 160 | yes | 41 | 2.05 | 41 | `{-3,-3^-1}` | yes |  |  | 112.719 |
| 20 | 17 | 13 | 160 | yes | 41 | 2.05 | 41 | `{-3,-3^-1}` | yes |  |  | 133.595 |

**0 mismatch(es) in 40 completed runs** (n = [17, 18, 19, 20], 0 censored).

syz₃ = 8n in every run.

### A.4 E12 and E06 extensions

#### A.4.1 `cross_relations.py n 4` — quartic relations on (t₀,t₁)

| n | exps | deg | points | dim on (t₀,t₁) | monomials |
|---:|---|---:|---:|---:|---:|
| 16 | 3;5 | 4 | 65534 | 33 | 41449 |

#### A.4.2 `explain_n16.py --deg 4` — the crude (and at degree 4 incorrect) test

| n | exps | deg | dim total | dim(t₀ only) | dim(t₁ only) | wt(−r₀) | wt(−r₁) | genuine cross |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 16 | 3;5 | 4 | 33 | 0 | 17 | 5 | 3 | **16** |

#### A.4.3 Subtracting the *ideal* generated by the image criteria (`server/jobs/e12_criteria_subtract.py`)

At degree 4 it is not enough to ask whether a relation is supported on t₀ only or on t₁ only: if g = Tr(t₁^−r₁) is a cubic relation then every product `g·t₀ᵢ` is a quartic relation touching both blocks that is *not* a cross relation.  `dim C` is the degree-≤deg part of the ideal generated by the criteria; C ⊆ R, so `dim R − dim C` is the number of genuine cross-S-box relations.

| n | deg | exps | points | monomials | deg g₀ (wt) | deg g₁ (wt) | dim R | dim C | genuine cross |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 14 | 4 | 3;5 | 16382 | 24158 | 9 (9) | 11 (11) | 7776 | 0 | **7776** |
| 16 | 3 | 3;5 | 65534 | 5489 | 5 (5) | 3 (3) | 1 | 1 | **0** |
| 16 | 4 | 3;5 | 65534 | 41449 | 5 (5) | 3 (3) | 33 | 33 | **0** |

#### A.4.4 E06 kernel bilinear model and its random-bilinear control

| n | system | D | nv | eqs | rank | syz | dim_low | falls from D−1 | 2·#eqs | dim_low/n |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 11 | P1only | 2 | 33 | 33 | 33 | 0 | 0 |  | 66 | 0.0 |
| 11 | P1only | 3 | 33 | 33 | 1122 | 0 | 99 | 66 | 66 | 9.0 |
| 11 | P1only | 4 | 33 | 33 | 17919 | 627 | 3102 | 1980 | 66 | 282.0 |
| 11 | S_quad | 2 | 33 | 99 | 99 | 0 | 0 |  | 198 | 0.0 |
| 11 | S_quad | 3 | 33 | 99 | 3102 | 264 | 99 | 0 | 198 | 9.0 |
| 11 | S_quad | 4 | 33 | 99 | 41469 | 14169 | 6016 | 2914 | 198 | 546.9091 |
| 11 | kernel | 2 | 66 | 66 | 66 | 0 | 0 |  | 132 | 0.0 |
| 11 | kernel | 3 | 66 | 66 | 4422 | 0 | 132 | 66 | 132 | 12.0 |
| 11 | kernel | 4 | 66 | 66 | 143715 | 2277 | 9075 | 4653 | 132 | 825.0 |
| 11 | random_bilinear | 2 | 66 | 66 | 66 | 0 | 0 |  | 132 | 0.0 |
| 11 | random_bilinear | 3 | 66 | 66 | 4422 | 0 | 66 | 0 | 132 | 6.0 |
| 11 | random_bilinear | 4 | 66 | 66 | 143781 | 2211 | 4422 | 0 | 132 | 402.0 |
| 11 | random_quadratic | 2 | 33 | 99 | 99 | 0 | 0 |  | 198 | 0.0 |
| 11 | random_quadratic | 3 | 33 | 99 | 3366 | 0 | 99 | 0 | 198 | 9.0 |
| 11 | random_quadratic | 4 | 33 | 99 | 46937 | 8701 | 6017 | 2651 | 198 | 547.0 |
| 13 | P1only | 2 | 39 | 39 | 39 | 0 | 0 |  | 78 | 0.0 |
| 13 | P1only | 3 | 39 | 39 | 1560 | 0 | 117 | 78 | 78 | 9.0 |
| 13 | S_quad | 2 | 39 | 117 | 117 | 0 | 0 |  | 234 | 0.0 |
| 13 | S_quad | 3 | 39 | 117 | 4368 | 312 | 117 | 0 | 234 | 9.0 |
| 13 | kernel | 2 | 78 | 78 | 78 | 0 | 0 |  | 156 | 0.0 |
| 13 | kernel | 3 | 78 | 78 | 6162 | 0 | 156 | 78 | 156 | 12.0 |
| 13 | random_bilinear | 2 | 78 | 78 | 78 | 0 | 0 |  | 156 | 0.0 |
| 13 | random_bilinear | 3 | 78 | 78 | 6162 | 0 | 78 | 0 | 156 | 6.0 |
| 13 | random_quadratic | 2 | 39 | 117 | 117 | 0 | 0 |  | 234 | 0.0 |
| 13 | random_quadratic | 3 | 39 | 117 | 4680 | 0 | 117 | 0 | 234 | 9.0 |
| 15 | P1only | 2 | 45 | 45 | 45 | 0 | 0 |  | 90 | 0.0 |
| 15 | P1only | 3 | 45 | 45 | 2070 | 0 | 135 | 90 | 90 | 9.0 |
| 15 | S_quad | 2 | 45 | 135 | 135 | 0 | 0 |  | 270 | 0.0 |
| 15 | S_quad | 3 | 45 | 135 | 5850 | 360 | 136 | 1 | 270 | 9.0667 |
| 15 | kernel | 2 | 90 | 90 | 90 | 0 | 0 |  | 180 | 0.0 |
| 15 | kernel | 3 | 90 | 90 | 8190 | 0 | 180 | 90 | 180 | 12.0 |
| 15 | random_bilinear | 2 | 90 | 90 | 90 | 0 | 0 |  | 180 | 0.0 |
| 15 | random_bilinear | 3 | 90 | 90 | 8190 | 0 | 90 | 0 | 180 | 6.0 |
| 15 | random_quadratic | 2 | 45 | 135 | 135 | 0 | 0 |  | 270 | 0.0 |
| 15 | random_quadratic | 3 | 45 | 135 | 6210 | 0 | 135 | 0 | 270 | 9.0 |

### A.5 Real solving degree `d_solv` (F4 step degrees)

#### A.5.1 `d_solv` vs n, grevlex, both tools

| ℓ | system | seed | n=7 | n=9 | n=11 | n=13 | n=15 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | P1only | 1 | 3 | 3 | 4 | 4 | 4 |
| 1 | P1only | 2 | 3 | 3 | 4 | 4 | 4 |
| 1 | S_quad | 1 | 3 | 3 | 4 | 4 | 4 |
| 1 | S_quad | 2 | 3 | 3 | 4 | 4 | 4 |
| 1 | random | 1 | 3 | 3 | 4 | 4 | 4 |
| 1 | random | 2 | 3 | 3 | 4 | 4 | 4 |
| 2 | P1only | 1 | 3 | 4 | 4 | 4 | – |
| 2 | P1only | 2 | 3 | 4 | 4 | – | – |
| 2 | S_quad | 1 | 3 | 4 | 4 | 4 | – |
| 2 | S_quad | 2 | 3 | 4 | 4 | – | – |
| 2 | random | 1 | 4 | 4 | 4 | – | – |
| 2 | random | 2 | 4 | 4 | 4 | – | – |

A cell `a/b` in bold would mean msolve and Magma disagreed; a plain number means they agreed (or only one tool ran).

#### A.5.2 Largest F4 matrix and cost, grevlex

| n | ℓ | system | seed | tool | d_solv | steps | max rows | max cols | s | peak RSS (MB) |
|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 7 | 1 | P1only | 1 | magma | 3 | 5 | 602 | 470 | 0.1 | 16 |
| 7 | 1 | P1only | 1 | msolve | 3 | 5 | 784 | 652 | 0.1 | 6 |
| 7 | 1 | P1only | 2 | magma | 3 | 5 | 602 | 470 | 0.1 | 16 |
| 7 | 1 | P1only | 2 | msolve | 3 | 5 | 785 | 653 | 0.1 | 6 |
| 7 | 1 | S_quad | 1 | magma | 3 | 4 | 631 | 470 | 0.1 | 16 |
| 7 | 1 | S_quad | 1 | msolve | 3 | 4 | 826 | 666 | 0.0 | 6 |
| 7 | 1 | S_quad | 2 | magma | 3 | 4 | 631 | 470 | 0.1 | 16 |
| 7 | 1 | S_quad | 2 | msolve | 3 | 4 | 826 | 666 | 0.1 | 6 |
| 7 | 1 | random | 1 | magma | 3 | 4 | 633 | 470 | 0.1 | 16 |
| 7 | 1 | random | 1 | msolve | 3 | 4 | 826 | 666 | 0.1 | 6 |
| 7 | 1 | random | 2 | magma | 3 | 4 | 632 | 470 | 0.1 | 16 |
| 7 | 1 | random | 2 | msolve | 3 | 4 | 826 | 666 | 0.1 | 6 |
| 7 | 2 | P1only | 1 | magma | 3 | 6 | 3496 | 1562 | 0.1 | 18 |
| 7 | 2 | P1only | 1 | msolve | 3 | 6 | 3811 | 1979 | 0.2 | 12 |
| 7 | 2 | P1only | 2 | magma | 3 | 6 | 3480 | 1562 | 0.1 | 18 |
| 7 | 2 | P1only | 2 | msolve | 3 | 6 | 3812 | 1979 | 0.2 | 12 |
| 7 | 2 | S_quad | 1 | magma | 3 | 5 | 3471 | 1562 | 0.1 | 18 |
| 7 | 2 | S_quad | 1 | msolve | 3 | 5 | 3804 | 2003 | 0.2 | 11 |
| 7 | 2 | S_quad | 2 | magma | 3 | 5 | 3462 | 1562 | 0.1 | 18 |
| 7 | 2 | S_quad | 2 | msolve | 3 | 5 | 3803 | 2003 | 0.2 | 11 |
| 7 | 2 | random | 1 | magma | 4 | 5 | 14960 | 5122 | 0.2 | 28 |
| 7 | 2 | random | 1 | msolve | 4 | 5 | 17021 | 7183 | 0.4 | 28 |
| 7 | 2 | random | 2 | magma | 4 | 5 | 15240 | 5402 | 0.2 | 28 |
| 7 | 2 | random | 2 | msolve | 4 | 5 | 17381 | 7543 | 0.5 | 28 |
| 9 | 1 | P1only | 1 | magma | 3 | 6 | 2096 | 988 | 0.1 | 17 |
| 9 | 1 | P1only | 1 | msolve | 3 | 6 | 2326 | 1296 | 0.1 | 9 |
| 9 | 1 | P1only | 2 | magma | 3 | 6 | 2119 | 988 | 0.1 | 17 |
| 9 | 1 | P1only | 2 | msolve | 3 | 6 | 2350 | 1296 | 0.1 | 9 |
| 9 | 1 | S_quad | 1 | magma | 3 | 5 | 2160 | 988 | 0.1 | 17 |
| 9 | 1 | S_quad | 1 | msolve | 3 | 5 | 2437 | 1312 | 0.1 | 9 |
| 9 | 1 | S_quad | 2 | magma | 3 | 5 | 2162 | 988 | 0.1 | 18 |
| 9 | 1 | S_quad | 2 | msolve | 3 | 5 | 2442 | 1312 | 0.1 | 9 |
| 9 | 1 | random | 1 | magma | 3 | 4 | 1029 | 988 | 0.1 | 18 |
| 9 | 1 | random | 1 | msolve | 3 | 4 | 1350 | 1312 | 0.1 | 10 |
| 9 | 1 | random | 2 | magma | 3 | 4 | 1029 | 988 | 0.1 | 18 |
| 9 | 1 | random | 2 | msolve | 3 | 4 | 1350 | 1312 | 0.1 | 10 |
| 9 | 2 | P1only | 1 | magma | 4 | 7 | 29050 | 20011 | 0.6 | 63 |
| 9 | 2 | P1only | 1 | msolve | 4 | 8 | 33898 | 25180 | 3.9 | 150 |
| 9 | 2 | P1only | 2 | magma | 4 | 7 | 29056 | 19995 | 0.6 | 63 |
| 9 | 2 | P1only | 2 | msolve | 4 | 8 | 34321 | 25546 | 5.3 | 150 |
| 9 | 2 | S_quad | 1 | magma | 4 | 6 | 28664 | 20030 | 0.5 | 62 |
| 9 | 2 | S_quad | 1 | msolve | 4 | 7 | 34138 | 25501 | 3.9 | 143 |
| 9 | 2 | S_quad | 2 | magma | 4 | 6 | 27606 | 18929 | 0.6 | 61 |
| 9 | 2 | S_quad | 2 | msolve | 4 | 7 | 33168 | 24483 | 4.7 | 141 |
| 9 | 2 | random | 1 | magma | 4 | 5 | 30062 | 15807 | 1.0 | 86 |
| 9 | 2 | random | 1 | msolve | 4 | 6 | 35325 | 21070 | 8.3 | 200 |
| 9 | 2 | random | 2 | magma | 4 | 5 | 30062 | 15807 | 0.8 | 88 |
| 9 | 2 | random | 2 | msolve | 4 | 6 | 35325 | 21070 | 7.8 | 201 |
| 11 | 1 | P1only | 1 | magma | 4 | 6 | 15758 | 8467 | 0.2 | 31 |
| 11 | 1 | P1only | 1 | msolve | 4 | 7 | 18313 | 11320 | 0.8 | 43 |
| 11 | 1 | P1only | 2 | magma | 4 | 6 | 15659 | 8229 | 0.2 | 32 |
| 11 | 1 | P1only | 2 | msolve | 4 | 7 | 18150 | 11039 | 0.9 | 43 |
| 11 | 1 | S_quad | 1 | magma | 4 | 5 | 15580 | 8625 | 0.2 | 31 |
| 11 | 1 | S_quad | 1 | msolve | 4 | 6 | 20180 | 13225 | 0.7 | 43 |
| 11 | 1 | S_quad | 2 | magma | 4 | 5 | 15503 | 8436 | 0.2 | 31 |
| 11 | 1 | S_quad | 2 | msolve | 4 | 6 | 18687 | 11620 | 0.8 | 41 |
| 11 | 1 | random | 1 | magma | 4 | 5 | 17077 | 6393 | 0.3 | 32 |
| 11 | 1 | random | 1 | msolve | 4 | 5 | 19596 | 8912 | 0.6 | 38 |
| 11 | 1 | random | 2 | magma | 4 | 5 | 16852 | 6194 | 0.3 | 32 |
| 11 | 1 | random | 2 | msolve | 4 | 5 | 19229 | 8571 | 0.7 | 37 |
| 11 | 2 | P1only | 1 | magma | 4 | 6 | 52505 | 45716 | 5.2 | 232 |
| 11 | 2 | P1only | 1 | magma | 4 | 6 | 52505 | 45716 | 4.1 | 232 |
| 11 | 2 | P1only | 1 | msolve | 4 | 7 | 62367 | 56041 | 54.2 | 961 |
| 11 | 2 | P1only | 1 | msolve | 4 | 7 | 62367 | 56041 | 48.2 | 966 |
| 11 | 2 | P1only | 2 | magma | 4 | 6 | 52229 | 45456 | 4.5 | 231 |
| 11 | 2 | P1only | 2 | msolve | 4 | 7 | 62624 | 56305 | 56.4 | 961 |
| 11 | 2 | S_quad | 1 | magma | 4 | 5 | 52152 | 45930 | 4.4 | 232 |
| 11 | 2 | S_quad | 1 | magma | 4 | 5 | 52152 | 45930 | 4.3 | 232 |
| 11 | 2 | S_quad | 1 | msolve | 4 | 6 | 65618 | 59396 | 52.7 | 923 |
| 11 | 2 | S_quad | 1 | msolve | 4 | 6 | 65618 | 59396 | 56.2 | 924 |
| 11 | 2 | S_quad | 2 | magma | 4 | 5 | 49156 | 42950 | 5.4 | 228 |
| 11 | 2 | S_quad | 2 | msolve | 4 | 6 | 60278 | 54072 | 69.4 | 914 |
| 11 | 2 | random | 1 | magma | 4 | 5 | 53245 | 38125 | 5.4 | 295 |
| 11 | 2 | random | 1 | msolve | 4 | 6 | 64108 | 48988 | 115.6 | 1593 |
| 11 | 2 | random | 2 | magma | 4 | 5 | 53245 | 38125 | 5.2 | 295 |
| 11 | 2 | random | 2 | msolve | 4 | 6 | 64108 | 48988 | 160.1 | 1590 |
| 13 | 1 | P1only | 1 | magma | 4 | 6 | 26151 | 16904 | 0.6 | 65 |
| 13 | 1 | P1only | 1 | msolve | 4 | 7 | 30561 | 21732 | 5.2 | 178 |
| 13 | 1 | P1only | 2 | magma | 4 | 6 | 26218 | 16971 | 0.6 | 65 |
| 13 | 1 | P1only | 2 | msolve | 4 | 7 | 31456 | 22624 | 6.9 | 180 |
| 13 | 1 | S_quad | 1 | magma | 4 | 5 | 25302 | 16529 | 0.5 | 63 |
| 13 | 1 | S_quad | 1 | msolve | 4 | 6 | 30722 | 21949 | 5.1 | 170 |
| 13 | 1 | S_quad | 2 | magma | 4 | 5 | 26217 | 17444 | 0.5 | 64 |
| 13 | 1 | S_quad | 2 | msolve | 4 | 6 | 34195 | 25422 | 5.2 | 176 |
| 13 | 1 | random | 1 | magma | 4 | 5 | 26832 | 13183 | 0.7 | 70 |
| 13 | 1 | random | 1 | msolve | 4 | 6 | 31330 | 17681 | 4.8 | 136 |
| 13 | 1 | random | 2 | magma | 4 | 5 | 27000 | 13337 | 0.6 | 71 |
| 13 | 1 | random | 2 | msolve | 4 | 6 | 31681 | 18018 | 4.0 | 137 |
| 13 | 2 | P1only | 1 | magma | 4 | 8 | 93423 | 90207 | 117.6 | 1647 |
| 13 | 2 | P1only | 1 | msolve | 4 | 9 | 111293 | 110272 | 1091.8 | 6496 |
| 13 | 2 | S_quad | 1 | magma | 4 | 7 | 92156 | 86658 | 131.7 | 1637 |
| 13 | 2 | S_quad | 1 | msolve | 4 | 8 | 110386 | 104993 | 1209.7 | 6342 |
| 15 | 1 | P1only | 1 | magma | 4 | 7 | 40388 | 30576 | 5.4 | 156 |
| 15 | 1 | P1only | 1 | msolve | 4 | 8 | 48871 | 39616 | 35.9 | 738 |
| 15 | 1 | P1only | 2 | magma | 4 | 7 | 40469 | 30657 | 6.8 | 156 |
| 15 | 1 | P1only | 2 | msolve | 4 | 8 | 48924 | 39667 | 36.5 | 734 |
| 15 | 1 | S_quad | 1 | magma | 4 | 6 | 40419 | 31243 | 6.3 | 156 |
| 15 | 1 | S_quad | 1 | msolve | 4 | 7 | 52910 | 43734 | 32.8 | 641 |
| 15 | 1 | S_quad | 2 | magma | 4 | 6 | 40428 | 31252 | 5.5 | 157 |
| 15 | 1 | S_quad | 2 | msolve | 4 | 7 | 51841 | 42665 | 47.4 | 639 |
| 15 | 1 | random | 1 | magma | 4 | 5 | 40512 | 25133 | 1.7 | 164 |
| 15 | 1 | random | 1 | msolve | 4 | 6 | 48345 | 32966 | 30.9 | 563 |
| 15 | 1 | random | 2 | magma | 4 | 5 | 40310 | 24946 | 1.6 | 164 |
| 15 | 1 | random | 2 | msolve | 4 | 6 | 48126 | 32762 | 29.3 | 562 |

### A.6 Block orders (E10)

Magma 2.24's `BooleanPolynomialRing` takes an order *name* only, so the block orders run on `PolynomialRing(GF(2), nv, "elim", S)` with the field equations added by hand.  That ring is slower than the Boolean one, so the baseline for these rows is `genref` (= grevlex on the *same* general ring), not the `grevlex` of §5.

| n | ℓ | system | seed | order | ring | d_solv | steps | max rows | s |
|---:|---:|---|---:|---|---|---:|---:|---:|---:|
| 7 | 1 | P1only | 1 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 85 | 0.1 |
| 7 | 1 | P1only | 1 | `blockT` (elimlast:7) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 85 | 0.1 |
| 7 | 1 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 85 | 0.1 |
| 7 | 1 | P1only | 2 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 85 | 0.1 |
| 7 | 1 | P1only | 2 | `blockT` (elimlast:7) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 85 | 0.1 |
| 7 | 1 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 85 | 0.1 |
| 7 | 1 | S_quad | 1 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 144 | 0.1 |
| 7 | 1 | S_quad | 1 | `blockT` (elimlast:7) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 144 | 0.1 |
| 7 | 1 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 144 | 0.1 |
| 7 | 1 | S_quad | 2 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 631 | 0.1 |
| 7 | 1 | S_quad | 2 | `blockT` (elimlast:7) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 631 | 0.1 |
| 7 | 1 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 4 | 631 | 0.1 |
| 7 | 1 | random | 1 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 147 | 0.1 |
| 7 | 1 | random | 1 | `blockT` (elimlast:7) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 147 | 0.1 |
| 7 | 1 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 147 | 0.1 |
| 7 | 1 | random | 2 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 147 | 0.1 |
| 7 | 1 | random | 2 | `blockT` (elimlast:7) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 147 | 0.1 |
| 7 | 1 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 147 | 0.1 |
| 7 | 2 | P1only | 1 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 3496 | 0.1 |
| 7 | 2 | P1only | 1 | `blockT` (elimlast:14) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 3496 | 0.1 |
| 7 | 2 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 3496 | 0.1 |
| 7 | 2 | P1only | 2 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 3480 | 0.1 |
| 7 | 2 | P1only | 2 | `blockT` (elimlast:14) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 3480 | 0.1 |
| 7 | 2 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 3480 | 0.1 |
| 7 | 2 | S_quad | 1 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 3471 | 0.1 |
| 7 | 2 | S_quad | 1 | `blockT` (elimlast:14) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 3471 | 0.1 |
| 7 | 2 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 3471 | 0.1 |
| 7 | 2 | S_quad | 2 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 3462 | 0.1 |
| 7 | 2 | S_quad | 2 | `blockT` (elimlast:14) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 3462 | 0.1 |
| 7 | 2 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 3462 | 0.1 |
| 7 | 2 | random | 1 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 294 | 0.3 |
| 7 | 2 | random | 1 | `blockT` (elimlast:14) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 294 | 0.3 |
| 7 | 2 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 294 | 0.2 |
| 7 | 2 | random | 2 | `blockPT` (elim:7) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 294 | 0.3 |
| 7 | 2 | random | 2 | `blockT` (elimlast:14) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 294 | 0.4 |
| 7 | 2 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 294 | 0.3 |
| 9 | 1 | P1only | 1 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 2096 | 0.1 |
| 9 | 1 | P1only | 1 | `blockT` (elimlast:9) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 2096 | 0.1 |
| 9 | 1 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 6 | 2096 | 0.1 |
| 9 | 1 | P1only | 2 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 129 | 0.2 |
| 9 | 1 | P1only | 2 | `blockT` (elimlast:9) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 129 | 0.2 |
| 9 | 1 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 129 | 0.2 |
| 9 | 1 | S_quad | 1 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 2160 | 0.1 |
| 9 | 1 | S_quad | 1 | `blockT` (elimlast:9) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 2160 | 0.1 |
| 9 | 1 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 2160 | 0.1 |
| 9 | 1 | S_quad | 2 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 2162 | 0.1 |
| 9 | 1 | S_quad | 2 | `blockT` (elimlast:9) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 2162 | 0.1 |
| 9 | 1 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 5 | 2162 | 0.1 |
| 9 | 1 | random | 1 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 225 | 0.1 |
| 9 | 1 | random | 1 | `blockT` (elimlast:9) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 225 | 0.1 |
| 9 | 1 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 225 | 0.1 |
| 9 | 1 | random | 2 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 225 | 0.1 |
| 9 | 1 | random | 2 | `blockT` (elimlast:9) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 225 | 0.1 |
| 9 | 1 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 3 | 3 | 225 | 0.1 |
| 9 | 2 | P1only | 1 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 29050 | 0.6 |
| 9 | 2 | P1only | 1 | `blockT` (elimlast:18) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 29050 | 0.7 |
| 9 | 2 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 29050 | 0.6 |
| 9 | 2 | P1only | 2 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 29056 | 0.6 |
| 9 | 2 | P1only | 2 | `blockT` (elimlast:18) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 29056 | 0.6 |
| 9 | 2 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 29056 | 0.6 |
| 9 | 2 | S_quad | 1 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 28664 | 0.6 |
| 9 | 2 | S_quad | 1 | `blockT` (elimlast:18) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 28664 | 0.6 |
| 9 | 2 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 28664 | 0.6 |
| 9 | 2 | S_quad | 2 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 27606 | 0.6 |
| 9 | 2 | S_quad | 2 | `blockT` (elimlast:18) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 27606 | 0.5 |
| 9 | 2 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 27606 | 0.6 |
| 9 | 2 | random | 1 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 459 | 0.8 |
| 9 | 2 | random | 1 | `blockT` (elimlast:18) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 459 | 0.8 |
| 9 | 2 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 459 | 0.8 |
| 9 | 2 | random | 2 | `blockPT` (elim:9) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 459 | 0.8 |
| 9 | 2 | random | 2 | `blockT` (elimlast:18) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 459 | 0.9 |
| 9 | 2 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 459 | 0.8 |
| 11 | 1 | P1only | 1 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 15758 | 0.3 |
| 11 | 1 | P1only | 1 | `blockT` (elimlast:11) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 15758 | 0.3 |
| 11 | 1 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 15758 | 0.3 |
| 11 | 1 | P1only | 2 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 186 | 0.7 |
| 11 | 1 | P1only | 2 | `blockT` (elimlast:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 186 | 0.7 |
| 11 | 1 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 186 | 0.7 |
| 11 | 1 | S_quad | 1 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 15580 | 0.3 |
| 11 | 1 | S_quad | 1 | `blockT` (elimlast:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 15580 | 0.3 |
| 11 | 1 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 15580 | 0.3 |
| 11 | 1 | S_quad | 2 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 15503 | 0.3 |
| 11 | 1 | S_quad | 2 | `blockT` (elimlast:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 15503 | 0.3 |
| 11 | 1 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 15503 | 0.2 |
| 11 | 1 | random | 1 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 319 | 0.3 |
| 11 | 1 | random | 1 | `blockT` (elimlast:11) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 319 | 0.3 |
| 11 | 1 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 319 | 0.4 |
| 11 | 1 | random | 2 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 319 | 0.3 |
| 11 | 1 | random | 2 | `blockT` (elimlast:11) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 319 | 0.3 |
| 11 | 1 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 319 | 0.3 |
| 11 | 2 | P1only | 1 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 52505 | 4.5 |
| 11 | 2 | P1only | 1 | `blockT` (elimlast:22) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 52505 | 5.9 |
| 11 | 2 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 52505 | 4.3 |
| 11 | 2 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 52505 | 4.9 |
| 11 | 2 | P1only | 2 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 52229 | 4.4 |
| 11 | 2 | P1only | 2 | `blockT` (elimlast:22) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 52229 | 4.3 |
| 11 | 2 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 52229 | 4.3 |
| 11 | 2 | S_quad | 1 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 52152 | 5.9 |
| 11 | 2 | S_quad | 1 | `blockT` (elimlast:22) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 52152 | 5.2 |
| 11 | 2 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 52152 | 6.5 |
| 11 | 2 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 52152 | 4.8 |
| 11 | 2 | S_quad | 2 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 49156 | 4.5 |
| 11 | 2 | S_quad | 2 | `blockT` (elimlast:22) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 49156 | 4.4 |
| 11 | 2 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 49156 | 4.6 |
| 11 | 2 | random | 1 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 660 | 2.4 |
| 11 | 2 | random | 1 | `blockT` (elimlast:22) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 660 | 2.4 |
| 11 | 2 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 660 | 2.4 |
| 11 | 2 | random | 2 | `blockPT` (elim:11) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 660 | 2.6 |
| 11 | 2 | random | 2 | `blockT` (elimlast:22) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 660 | 2.7 |
| 11 | 2 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 660 | 2.7 |
| 13 | 1 | P1only | 1 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 26151 | 0.6 |
| 13 | 1 | P1only | 1 | `blockT` (elimlast:13) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 26151 | 0.6 |
| 13 | 1 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 26151 | 0.6 |
| 13 | 1 | P1only | 2 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 26218 | 0.6 |
| 13 | 1 | P1only | 2 | `blockT` (elimlast:13) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 26218 | 0.6 |
| 13 | 1 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 26218 | 0.6 |
| 13 | 1 | S_quad | 1 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 25302 | 0.6 |
| 13 | 1 | S_quad | 1 | `blockT` (elimlast:13) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 25302 | 0.6 |
| 13 | 1 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 25302 | 0.6 |
| 13 | 1 | S_quad | 2 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 26217 | 0.6 |
| 13 | 1 | S_quad | 2 | `blockT` (elimlast:13) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 26217 | 0.6 |
| 13 | 1 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 5 | 26217 | 0.6 |
| 13 | 1 | random | 1 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 429 | 0.6 |
| 13 | 1 | random | 1 | `blockT` (elimlast:13) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 429 | 0.6 |
| 13 | 1 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 429 | 0.7 |
| 13 | 1 | random | 2 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 429 | 0.6 |
| 13 | 1 | random | 2 | `blockT` (elimlast:13) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 429 | 0.6 |
| 13 | 1 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 429 | 0.6 |
| 13 | 2 | P1only | 1 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 8 | 93423 | 245.6 |
| 13 | 2 | P1only | 1 | `blockT` (elimlast:26) | PolynomialRing(GF(2)) + field equations | 4 | 8 | 93423 | 183.8 |
| 13 | 2 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 8 | 93423 | 155.3 |
| 13 | 2 | S_quad | 1 | `blockPT` (elim:13) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 92156 | 169.0 |
| 13 | 2 | S_quad | 1 | `blockT` (elimlast:26) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 92156 | 119.4 |
| 13 | 2 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 92156 | 144.2 |
| 15 | 1 | P1only | 1 | `blockPT` (elim:15) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 40388 | 6.2 |
| 15 | 1 | P1only | 1 | `blockT` (elimlast:15) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 40388 | 6.6 |
| 15 | 1 | P1only | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 40388 | 6.5 |
| 15 | 1 | P1only | 2 | `blockPT` (elim:15) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 40469 | 6.3 |
| 15 | 1 | P1only | 2 | `blockT` (elimlast:15) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 40469 | 6.6 |
| 15 | 1 | P1only | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 7 | 40469 | 6.5 |
| 15 | 1 | S_quad | 1 | `blockPT` (elim:15) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 40419 | 6.6 |
| 15 | 1 | S_quad | 1 | `blockT` (elimlast:15) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 40419 | 6.6 |
| 15 | 1 | S_quad | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 40419 | 6.5 |
| 15 | 1 | S_quad | 2 | `blockPT` (elim:15) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 40428 | 6.5 |
| 15 | 1 | S_quad | 2 | `blockT` (elimlast:15) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 40428 | 6.5 |
| 15 | 1 | S_quad | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 6 | 40428 | 6.1 |
| 15 | 1 | random | 1 | `blockPT` (elim:15) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 555 | 1.4 |
| 15 | 1 | random | 1 | `blockT` (elimlast:15) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 555 | 1.4 |
| 15 | 1 | random | 1 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 555 | 1.5 |
| 15 | 1 | random | 2 | `blockPT` (elim:15) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 555 | 1.4 |
| 15 | 1 | random | 2 | `blockT` (elimlast:15) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 555 | 1.5 |
| 15 | 1 | random | 2 | `genref` (gen:grevlex) | PolynomialRing(GF(2)) + field equations | 4 | 4 | 555 | 1.3 |

### A.7 E05 — resultant of the two-valued-inverse reduced system

| n | exps | ε | Eq1 terms | Eq2 terms | deg_t1(Eq1) | deg_t1(Eq2) | resultant s | deg | terms | roots | key is a root | factors | peak RSS (GB) | censored |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|:--:|---:|---:|:--:|
| 8 |  |  |  |  |  |  |  |  |  |  |  |  | 0.0 | yes |

