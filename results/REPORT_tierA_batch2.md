# Tier-A batch 2 — report

**Date** 2026-08-31 · **Branch** `claude/tier-a-batch-2` (from `main` after batch 1 merged, `daa531e`)
· **Tier** A only, pure Python standard library, one 4-core VM.

Write scope respected: only `docs/`, `paper/`, `src/aim3/`, the scripts and READMEs of
E01/E03/E05/E06/E09/E11, `results/E0*_theory_*.csv` and this file. `legacy/`, `server/`,
`results/server/` and the E04/E12 scripts were not touched; no `--backend` option was added.
`tests/` gained two regression tests for the new `src/aim3` modules (8/8 pass).

---

## 1. Wording corrections (task 1)

| where | was | now |
|---|---|---|
| `docs/ROADMAP.md` Status update | "核模型 D=3 降次 = 2·#eqs = 12n" | dim(span∩deg≤2) = 2·#eqs = 12n, **of which 6n is the 6n input equations themselves**; the *new* quadratics are #eqs = 6n |
| `experiments/E06_kernel_bilinear_model/README.md` | "84 = 2·#eqs degree falls"; "D=3 fall count is 2·#eqs = 12n" | same correction; the results table gained an explicit **new** column (= `falls_from_prev_D` in the CSV) |
| `paper/claims_ledger.md` rows L, L3 | as above | as above; L3 restated as "the control's degree-≤2 part is exactly its 6n input equations, i.e. 0 new" |

No datum changed: the CSV column `falls_from_prev_D` always held the correct quantity (42 at n=7,
54 at n=9 — i.e. 6n). `experiments/E12_cross_sbox_relations/README.md` already states the n=9
`cross_excess = 1` residual and its explanation (both exponents degenerate, per-S-box counts 5n and
7n, a single element that does not grow with n) — checked, no change needed, and E12 is outside this
batch's write scope in any case.

**One stale statement left in place:** `results/REPORT_tierA_batch1.md` §3.4 and §5.8 still use the
old phrasing. That file is outside this batch's write scope and is the historical record of batch 1;
the correction is recorded here and in the three files above.

## 2. The degree-fall generators (task 2) — the intended description is refuted

Full write-up: **`docs/THEORY_NOTES.md`**. Script:
`experiments/E03_syzygy_and_degree_fall/fall_generators.py` (`--sbox`, `--kernel`, `--orbits`,
`--all`), new field-polynomial layer `src/aim3/fieldpoly.py`, new `fall_space` / `span_rank` /
`reduce_modulo` / `labelled_rows` in `src/aim3/boolmac.py`. Data:
`results/E03_theory_fall_generators.csv`, `results/E03_theory_orbit_lemma.csv`,
`results/E06_theory_fall_generators.csv`.

### 2.1 What is proved

Writing $\operatorname{Tr}_\gamma$ for the coordinate functionals, every row of $\mathrm{Mac}_4$ is
$\operatorname{Tr}(\gamma\,M\,P_k)$ for a **field** monomial $M=x^ay^b$ with wt(a)+wt(b) ≤ 2 (a
Boolean bit is a trace, so a degree-≤2 bit-monomial is a degree-≤2 *function*, and those are exactly
the $\mathbb F_{2^n}$-span of such $M$). The degree-4 part of $\operatorname{Tr}(\gamma
x^{2^i}y^{2^j}P_1)$ is $\operatorname{Tr}(\gamma x^{2^i+1}y^{2^j+1})$, homogeneous of degree 4, so a
fall inside that family is a pure Frobenius coincidence.

**Orbit lemma [PROVED, checked n=5..20].** For 1 ≤ i,j ≤ n−1 the pairs $(2^i+1,2^j+1)$ and
$(2^{i'}+1,2^{j'}+1)$ are in one Frobenius orbit iff $(i',j')=(i,j)$, or $i=j$ and
$(i',j')=(n-i,n-i)$. Hence exactly $c=\lfloor (n-1)/2\rfloor$ coincidences, and the ansatz module
($n^2+2n$ generators × $n$ coefficients) has

> top-degree rank $= n((n-1)^2-c)$,  kernel $= n(n^2+2n-(n-1)^2+c)$.

At n=11 that predicts rank 1045 and kernel 528 — and the measurement gives **exactly** 1045 and 528.
The accounting is complete: 473 = 11·43 kernel dimensions from generators with no degree-4 part,
55 = 11·5 from the diagonal coincidences.

### 2.2 What is refuted

**Every one of those kernel elements already lies in $\mathrm{Mac}_3$.** The module proposed for
Theorem T1.2 — $\{x^{2^i}y^{2^j}P_1\}\cup\{x^{2^k}P_2\}\cup\{y^{2^k}P_3\}$ — explains **0 of the 2n**
new cubics, at (n,e) = (11,3) and (13,3). Nor does any single one of the 18 multiplier families
$\{1,x^{2^i},y^{2^j},x^{2^i}y^{2^j},x^{2^i+2^j},y^{2^i+2^j}\}\times\{P_1,P_2,P_3\}$ produce a new
cubic on its own — all 18 measured, all 0.

**What does work — leave-one-out.** Dropping one family at a time from the full row space: 16 of the
18 leave the count at 2n; dropping `x^{2^i}y^{2^j}·P₂` leaves n, and dropping `y^{2^i+2^j}·P₃` leaves
n. So **falls₄ = 2n splits as n + n**, and both blocks need **weight-2** multipliers on P₂/P₃ — the
refuted ansatz used weight-1 ones. [VERIFIED (11,3), all 18 families measured.]

**And why those two.** Written as multiples of P₁, `x^{2^k}P₂` and `y^{2^k}P₃` have *no* degree-4 part
(Boolean degree ≤3 — they are inside Mac₃, which is the first half of why the ansatz fails);
`x^{2^i}y^{2^j}P₁` has exactly *one* degree-4 monomial, so its only cancellations are the ⌊(n−1)/2⌋
diagonal Frobenius coincidences, and those were measured to land in Mac₃ too (the second half); while

    x^{2^i}y^{2^j}P₂ = x^{2^i+2}y^{2^j+1} + **x^{2^i+1+2^e}y^{2^j}** + x^{2^i+1}y^{2^j}
    y^{2^i+2^j}P₃   = x y^{2^i+2^j+2}    + **x^{2^e}y^{2^i+2^j+1}**  + y^{2^i+2^j+1}

each carry **two** degree-4 monomials, the second involving 2^e. Cancellation can happen inside one
generator, and the condition sees the exponent — **which is how falls₄ can depend on the class of e
at all**; a single-monomial family could not. This is the concrete starting point Theorem T1.2 needs,
and it is in `docs/THEORY_NOTES.md` §2.6–2.7 with the symbolic verification.

**A structural obstruction, worth stating in the paper.** falls₄ = 2n invites "two field cubics,
n coordinates each". That cannot be the general shape: the classes {−3,(−3)⁻¹} and {n/2−1} give
**2n+1**, and 2n+1 is not a multiple of n. Any correct description must produce a lone extra element
in those classes, not a third n-block.

### 2.3 The kernel model's 6n new quadratics

[VERIFIED n=7,9] The 6n new quadratics the D=3 pass adds to the kernel model are exactly the
coordinates of

    a_i·K_i = t_i a_i² + a_i^{2^{e_i}+1}   and   t_i·K_i = t_i² a_i + t_i a_i^{2^{e_i}},  i = 0,1,2

(with t₂ = pt + ct) — six field polynomials, n coordinates each. Unexplained remainder: **0**.

They do not raise the degree because a_i² and t_i² are F₂-linear and wt(2^e+1) = 2: **no cancellation
is involved**, unlike the D=4 falls. Two consequences: (i) it is the **K equations alone** that
produce all 6n — the B equations contribute none, so the brief's guess "t_i·(B) against a_i·(K)" is
half right; (ii) it explains E06's control — a random bilinear system has no Frobenius structure, so
multiplying by a variable really does give a cubic, which is why the control measures 0.

The S-box relations P₁⁽ⁱ⁾ are *implied* by the model (B_i + x_i·K_i = a_i·P₁⁽ⁱ⁾, one line) but
P₁⁽ⁱ⁾ itself is **not** in the D=3 row space and explains **none** of the 6n — only the cubic
multiple a_i·P₁⁽ⁱ⁾ is there.

## 3. ROADMAP §3.3 rewritten (task 3)

§3.3 now describes what `reduced_system_toy.py` actually builds, with a note at the top saying why
the previous formulation cannot be expanded: its second equation composes a non-linearised sparse
polynomial with another one (≈ 2ⁿ terms). The new text gives the mixed form (two-valued inverse on
the first two S-boxes, the last S-box forward through its own P₁), **2 branch bits rather than 3**,
the measured size table (Eq1 ≈ n² terms, Eq2 ≈ n³, deg_{t₀} = deg_{t₁} = 2ⁿ−2, Bézout ≈ 2^{2n}),
and the pseudo-root counts. T3.1 is upgraded to [VERIFIED n=8..14]; T3.2 records that the mixed form
hands the attacker one bit (4 branches, not 8) without changing the 2^{2n} order. **Only §3.3 was
touched**; the Status-update line changed from "needs rewriting" to "rewritten".

## 4. The bit cost of one AIM3 evaluation (task 4)

The spec PDF was reachable and Table 5 (p.42) was read: per forward evaluation,

| scheme | F_{2ⁿ} mul (count, depth) | F_{2ⁿ} square (count, depth) | Grover total (spec) |
|---|---|---|---|
| AIM3-I | 47, 24 | 399, 254 | 2^165.6 |
| AIM3-III | 59, 26 | 651, 382 | 2^232.5 |
| AIM3-V | 73, 28 | 1056, 510 | 2^298.0 |

**Independent cross-check of the reading:** the same addition-chain rule (Itoh–Tsujii needs
⌊log₂(n−1)⌋+wt(n−1)−1 multiplications, the Mersenne part ⌊log₂ e⌋+wt(e)−1, plus one per S-box) gives
50 / 63 / 77 with no sharing, against the spec's 47 / 59 / 73 — the 3–4 gap is the sharing between
the two x^{2^k−1} towers. (Note the *reference implementation* uses the naive n−1 multiplications,
not Itoh–Tsujii; Table 5 is about the optimised circuit, and it is the right input here.)

With an explicit bit model (squaring and reduction 4n; multiplication n²+(n−1)²+4n schoolbook or
6·n^{log₂3}+4n Karatsuba; n×n matrix–vector n²+n(n−1); addition n) and the observation that in a
single-key search the iv is fixed so the SHAKE tape is **not** a per-candidate cost:

| scheme | bits/eval (schoolbook) | (Karatsuba) | log₂ | B3 = 2ⁿ · that |
|---|---|---|---|---|
| AIM3-I | 1 887 791 | 976 414 | 20.85 / 19.90 | 2^148.85 / 2^147.90 |
| AIM3-III | 5 167 931 | 2 312 744 | 22.30 / 21.14 | 2^214.30 / 2^213.14 |
| AIM3-V | 11 473 993 | 4 816 758 | 23.45 / 22.20 | 2^279.45 / 2^278.20 |

So one evaluation costs 2^19.9–2^23.5 bit operations — **6 to 9 bits more than the retired 2^14
placeholder**. Written into `docs/COMPLEXITY_ACCOUNTING.md` §1 (with the model spelled out) and §2
(B3 in bit operations); `experiments/E11_hybrid_estimates/hybrid_estimates.py` now prints the real
B3 instead of 2^{n+14}. New: `aim3_eval_cost.py`, `results/E11_theory_eval_cost.csv`.

## 5. Controls and scenarios (task 5)

### 5.1 The linear-layer distribution does not matter

`experiments/E01_sbox_structures/matrix_distribution_control.py` rebuilds the E12 first-layer system
(6n quadratics in 3n variables) with unit-LU, uniform GL_n(F₂), and random invertible of the same
density. `--check` verifies the builder reproduces `E12/first_layer_ideal.py` bit for bit, so a
difference would have to be the distribution.

| n | distribution | D=3 rank | syz₃ (2·8n) | D=4 falls (4n) | saturation |
|---:|---|---:|---:|---:|---:|
| 9 | unit-LU / GL_n / density | 1368 / 1368 / 1368 | 144 / 144 / 144 | 109 / 109 / 108 | 0.447 |
| 11 | unit-LU / GL_n / density | 2068 / 2068 / 2068 | 176 / 176 / 176 | **44 / 44 / 44** | 0.351 |

At n=11 — non-saturated, generic exponents — all three are **identical**, including falls₄ = 4n.
So unit-LU is not a weakness and not a convenience assumption: everything measured in batches 1–2 is
a property of the S-box. (ROADMAP §4.3's control requirement, discharged for the first layer.)

### 5.2 Reused key: the coupling is exactly C(q,2)·n

`experiments/E09_multitarget_weak_iv/reused_key_toy.py`: q instances sharing one pt, each with its
own iv/matrices/ct. Against the direct sum of q independent instances:

| n | ℓ | q | syz₃ | direct sum | **excess** | C(q,2)·n | saturation |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 9 | 2 | 2 | 441 | 432 | **9** | 9 | 0.167 |
| 9 | 2 | 3 | 675 | 648 | **27** | 27 | 0.147 |
| 9 | 2 | 4 | 918 | 864 | **54** | 54 | 0.133 |
| 7 | 2 | 2 | 343 | 336 | **7** | 7 | 0.910 |
| 7 | 2 | 3 | 525 | 504 | **21** | 21 | 0.948 |
| 9 | 1 | 2 | 297 | 288 | **9** | 9 | 0.918 |
| 7 | 1 | 2 | 497 | 224 | 273 | 7 | **0.996 saturated** |
| 9 | 1 | 4 | 717 | 576 | 141 | 54 | **0.947** |

**syz₃ excess = C(q,2)·n — one n-block of syzygies per *pair* of instances** [CONJ, verified in the
six non-saturated runs]. The two most saturated runs deviate and are reported as such rather than
averaged in. Reading: key reuse does couple the instances, but only linearly in n per pair; it is a
quantified gain, not a collapse. This is the first concrete number for track T6.

### 5.3 E11

The kernel model's 6n/6n shape was added to the generic curves (guess k bits of U, so n_U = 3n−k,
n_W = 3n), capped at k ≤ n because pt determines every other variable. Its generic FSS bound stays
far above B3 (2^1363 vs 2^148.9 at n=128) — consistent with E11's existing conclusion that
structure-agnostic bilinear bounds cannot beat brute force. E06 shows the model is *not* generic, so
this curve is an upper bound on a structure-blind solver, not a prediction.

## 6. Paper Base part (task 6)

`paper/main.tex` (stock `article` + `amsthm` so it compiles anywhere; one comment line switches it to
`iacrtrans`) and `paper/sections/{02_prelim,03_sbox,04_ideal}.tex`.

* **02_prelim** — AIM3 by the fact-sheet references, the encoding conventions and why toy experiments
  are convention-independent, Boolean ideals / Macaulay matrices / degree falls, an explicit
  *saturation* remark, and the two semi-regular baselines (E02).
* **03_sbox** — structures A–E with **full proofs**: closed form; image and trace criterion (including
  why Tr(y^{−(r+1)}) = Tr(y^{−r}), which is where the criterion's usual statement comes from); kernel
  variable and the alternating form; the polarised identity and the difference-quotient form of the
  S-box relation; the weight theorem r = Σ 2^{ie mod n}, wt(r) = e⁻¹ mod n, with the telescoping proof.
* **04_ideal** — ⟨P₁,P₂,P₃⟩ = ⟨P₁⟩ with proof; the 3n structural syzygies with proof and the 5n as
  [VERIFIED]+[OPEN]; first-fall degree 4 and the class table; no cubic generators outside the ideal
  and the e=2 exception; §2's refutation as a numbered observation; the cross-S-box theorem with the
  n−m explanation of the n=16 dimension-1; the direct sum; the saturation remark for E04; and the
  matrix-distribution observation.

Every numbered statement carries a `%% ledger: <ID>` comment; a cross-check confirms no dangling IDs.
No attack complexity numbers appear. Server-dependent results are `\todo{}`.

## 7. Not done, and why

* **The 5n non-trivial syzygies at D=3 are still [OPEN].** §2's refutation is a reason to expect this
  to be hard too, and the leave-one-out tooling is now in place to test any proposed family cheaply.
* **The class rule and the C(q,2)·n rule remain [CONJ]** — both are fits within n ≤ 16 / q ≤ 4.
* **`refs.bib` is unchanged**; `main.tex` cites `AIMerV3` which must be added before it compiles with
  bibliography. The skeleton compiles as a document; the bibliography is a one-line fix left to
  whoever adds the third-party citations (roadmap §6 requires checking each against the original PDF).
* Sections 5, 6, 8 of the paper are `\todo{}` — they depend on the server batch.
* **`results/README.md` was not updated** with rows for the five new `E0*_theory_*.csv` files. That
  file is not in this batch's write scope (the scope lists `results/E0x_theory_*.csv` and this report
  only), and it is shared with the server branch. The new files are documented here and in
  `docs/THEORY_NOTES.md`; adding the rows is a one-commit follow-up for whoever merges the two
  branches.

## 8. Open decisions for a human

1. **Theorem T1.2 as stated in ROADMAP §3.1 should be re-scoped.** It asks for "the explicit
   generators of W" with the candidate mechanism "(x^{2^i}y^{2^j})·P₁ and x^{2^k}·P₂, y^{2^k}·P₃
   cancelling pairwise along Frobenius orbits". That candidate is now measured and refuted (§2.2).
   The theorem should either drop the candidate mechanism or be restated as a structural claim about
   the quotient (with the 2n+1 obstruction of §2.2 built in). This batch did **not** edit §3.1.
2. **Is C(q,2)·n worth a server run?** It is the first quantitative T6 result and would want q up to
   6 and n = 13 to be extrapolatable under COMPLEXITY_ACCOUNTING §4.
3. **Bibliography and paper template**: `iacrtrans` vs LNCS, and who supplies `refs.bib` entries.
