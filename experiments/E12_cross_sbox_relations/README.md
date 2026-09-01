# E12 — Cross-S-box relations of degree ≥ 3

The spec proves no *quadratic* relation exists between two S-box outputs t_i, t_j. `cross_relations.py 14 3` shows: **no cubic
relation on (t₀,t₁) either** (dim 0 of 3683 monomials on 16382 points). On (pt,t₀,t₁) the cubic relation space has dimension
3444 at n=14; the next step is to compare this with the degree-3 Macaulay rank of the union of the two S-box ideals in 3n
variables (extend `legacy/aim3_macaulay.py` to a first-layer-only system) — any excess = cross degree falls, feeding E04.
Also run deg=4 on (t₀,t₁) at n=16 (needs 2^n > C(32,≤4)=41449 ✓; use M4RI on the server for the 65k×41k rank).

---

## Log — 2026-08-31 (Tier A, batch 1)

**Status.** No cross-S-box cubic relation [VERIFIED n=14 **and n=16**]. The first-layer system's D=4
degree falls are an **exact direct sum** of the two S-box contributions [VERIFIED n=11, and
ideal-theoretically at n=14]. This is the Tier-A answer to **Gate 1**.

**What was run.**
* `first_layer_ideal.py` (new): the two first-layer S-boxes only — 6n quadratics {P1,P2,P3} in the
  3n variables (pt, t₀, t₁), built on exactly the instance that `cross_relations.py n deg seed`
  uses (same random stream). n = 9, 11 with `--D4`, n = 14 at D=3. →
  `results/E12_first_layer_ideal.csv`.
* `cross_relations.py 16 3` → `results/E12_cross_relations.csv`.
* `explain_n16.py` (new) → `results/E12_cross_relations_explained.csv`.

**Result 1 — n=14, the ideal-theoretic comparison the batch asked for.**

| quantity | value |
|---|---|
| D=3 Macaulay rank of the first-layer system (42 vars) | 3388 |
| syz₃ | 224 = 2·8n |
| dim of all cubics vanishing on {(pt,t₀,t₁)} (`cross_relations.py 14 3`) | 3444 |
| **excess = 3444 − 3388** | **56** |
| direct-sum prediction 2·(2n) = 4n | **56** |
| difference | **0** |

So the cubics of the joint vanishing ideal that D=3 does not yet produce are exactly the 2n + 2n
first-fall cubics of the two S-boxes — no more. [VERIFIED n=14]

**Result 2 — the per-S-box direct-sum control at D=4.** Comparing the joint fall count against 4n is
*not* by itself a direct-sum test, because the joint Macaulay matrix may multiply by cross-block
monomials that do not exist in E03's 2n-variable ring. Running each S-box's own 3n quadratics alone
in the same 3n variables fixes that:

| n | exps | falls(both) | falls(S₀) | falls(S₁) | **cross_excess** | verdict |
|---:|---|---:|---:|---:|---:|---|
| 11 | 3, 4 | 44 | 22 | 22 | **0** | **direct sum** |
| 9 | 4, 5 | 109 | 45 | 63 | 1 | see below |

At n = 11 the two exponents are generic (E03: 2n falls each), the run is far from saturated
(dim_low/C(33,≤3) = 2112/6018 = 0.35), and the D=4 falls of the joint system are **exactly** the sum
of the per-S-box falls — and also exactly 4n. At n = 9 both exponents are degenerate in the
E03-calibrated sense (e = 4 has 2e ≡ −1 mod 9, e = 5 has 2e ≡ +1), which is why the per-S-box counts
are 5n and 7n rather than 2n; the residual cross_excess is 1, i.e. a single element, not a growing
term. **No run shows cross-S-box falls growing with n.**

**Result 3 — n=16 reports dim = 1, and it is not a cross relation.** `cross_relations.py 16 3` gives
dimension 1 (not 0) for the cubics on (t₀,t₁) with e = (3,5). `explain_n16.py` decides it: the
relation space is supported on the **t₁ variables alone**, and Tr(t₁^{−r₁}) — the image criterion
of S_{e₁} — vanishes on every graph point and has Boolean degree wt(−r₁) = n − (e⁻¹ mod n) =
16 − 13 = **3**. So the single cubic is the *univariate* image criterion, not a cross-S-box relation.
At n = 14 the same criteria have degree 9 and 11, which is why the dimension there is 0.

This sharpens claim K: **the cubic relation space on (t₀,t₁) has dimension #{i : wt(−rᵢ) ≤ 3}, i.e.
it contains only univariate image criteria; there is no cross-S-box cubic** [VERIFIED n=14, 16].
It also links E12 to the E08 weight criterion: the low-degree relations that do exist are governed
by wt(−r) = n − (e⁻¹ mod n), exactly the quantity the exponent audit scores.

**Next steps.** Tier B: deg = 4 on (t₀,t₁) at n = 16 (needs 2^16 > C(32,≤4) = 41449 ✓; M4RI for the
65k×41k rank) — the same wt(−r) argument predicts that univariate criteria with wt(−rᵢ) ≤ 4 appear
there, so the control must subtract them before any cross-relation claim; the first-layer D=4 run at
n = 13, 14 (the Python run at n=14, D=4 is ~10⁵ rows × 1.2·10⁵ columns — server work).

## Server log — 2026-08-31 (server batch 1, branch `claude/server-batch-1`)

**First-layer system.** `first_layer_ideal.py 13 --D4` and `14 --D4` with `--backend sage`
(`results/E12_server_first_layer.csv`).

**Degree 4 on (t₀,t₁).** `cross_relations.py 16 4 --backend sage` (65534 × 41449 over F₂) and
`explain_n16.py 16 --deg 4 --backend sage`.  Both scripts gained the flags they needed; the
committed Tier-A CSVs are untouched (server rows go to `results/E12_server_*.csv`).

**A correction that matters.** `explain_n16.py` decides "univariate" by asking whether a relation is
supported on t₀-only or t₁-only monomials.  At degree 3 that is right.  At degree **4** it is not:
if g = Tr(t₁^{−r₁}) is a cubic relation, then every product `g·t₀ᵢ` is a quartic relation that
touches both blocks and is *not* a cross relation.  Read naively, the n = 16 deg-4 run reports 16
"genuine cross" relations; all 16 are exactly those products.  `server/jobs/e12_criteria_subtract.py`
does what this README prescribes instead — it builds the criteria from their ANF (so their Boolean
degree is measured, not assumed), spans the whole degree-≤4 part of the ideal they generate, and
subtracts it.  Rows: `results/E12_server_criteria_subtract.csv`.

**A size caveat, stated so it is not misread later.** At n = 14, deg = 4 the graph has 16382 points
but C(28,≤4) = 24158 monomials, so *any* instance has ≥ 7776 relations purely from the point
deficit; the measurement is meaningless there and is reported as such.  n = 16 (65534 points,
41449 monomials) is the smallest size at which the deg-4 question can be asked at all.

Result: `results/REPORT_server_batch1.md` §4.
