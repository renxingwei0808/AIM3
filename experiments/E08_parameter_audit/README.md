# E08 — Exponent and linear-layer audit (track T5)

**Two audit criteria (both new, both closed-form).**
1. *Inverse-exponent weight theorem*: wt(r) = e⁻¹ mod n, wt(−r) = n − wt(r), supp(r) = {ie mod n}. Controls the Boolean degree of
   the kernel map y ↦ y^r (two-valued inverse) and of the image test Tr(y^{−r}). Design rule: e⁻¹ mod n ≈ n/2.
2. *Frobenius-orbit degeneracy* (from E03): e = 2; e = 1 + n/2; k·e ≡ ±1, ±2 (mod n) for small k. These inflate syzygies/degree
   falls (4n–9n instead of 2n).

**Recommended parameters under the audit** (`exponent_audit.py`):
AIM3-I: e=3 (43/85), 7 (55/73), 11 (**35**/93); AIM3-III: 11 (35/157), 23 (167/**25**), 47 (143/**49**); AIM3-V: 19 (**27**/229).
None is degenerate in sense 2; several are unbalanced in sense 1 (AIM3-V final S-box e=19 has kernel-map degree 27 = 10.5 % of n).

**Next steps.** (a) Run the full audit for n = 128/192/256 (done, `results/E08_audit_*.csv`) and for toy n to correlate the
degeneracy tags with measured fall counts (E03) — target: a *theorem* "fall count = 2n + Δ(e)" with Δ explicit.
(b) Linear-layer distributions: unit-LU vs uniform GL_n vs same-density random — compare Macaulay ranks and falls (E04 inputs).
(c) Output: a "weak (n,e) classification" section for the paper and, if a class contains a *cheaper attack*, an AIM3-parameter
warning in the style of ePrint 2025/2272.

---

## Log — 2026-08-31 (Tier A, batch 1)

**Status.** Audit tables complete for n = 128, 192, 256 [VERIFIED]; the degeneracy criterion is now
**calibrated against measured fall counts** rather than guessed.

**What was run.** `python3 exponent_audit.py` (all three n) →
`results/E08_audit_{128,192,256}.csv`, now with two extra columns `falls_D4_pred` and `falls_over_n`
from `aim3.structures.predicted_falls_D4`.

**The degeneracy criterion changed.** E03's full sweep (n = 11..16, every admissible e, 36 runs, see
`results/E03_summary.md`) shows that the previous rule — tag k·e ≡ ±1, ±2 (mod n) for k ≤ 4 — is
wrong in both directions: 11 tagged exponents are perfectly generic and 3 untagged ones have 4n
falls. The measured fall count depends only on the pair {e, e⁻¹ mod n}:

| class | falls_D4 | # admissible e at n=128 / 192 / 256 |
|---|---|---|
| {2, 2⁻¹} | 7n | 0 / 0 / 0 (needs n odd) |
| {−2, (−2)⁻¹} | 4n | 0 / 0 / 0 (needs n odd) |
| {1+n/2} | 17n/2 | 1 (e=65) / 1 (e=97) / 1 (e=129) |
| {−3, (−3)⁻¹} | 2n+1 | 2 (e=85, 125) / 0 / 2 (e=85, 253) |
| {n/2−1} | 2n+1 | 1 (e=63) / 1 (e=95) / 1 (e=127) |
| generic | 2n | all the rest |

**Recommended parameters under the audit.** All ten recommended exponents (AIM3-I 3,7,11;
AIM3-III 11,23,47; AIM3-V 3,7,11,19) are in the **generic** class, i.e. exactly 2n falls per S-box —
now asserted by `exponent_audit.py` itself and by `tests/test_core.py::test_weight_theorem_recommended`.
Criterion 1 (weight balance) is unchanged and still flags AIM3-V's final S-box e=19
(wt(r) = 27 = 10.5 % of n) and AIM3-III's e=23 (wt(−r) = 25).

**Cross-reference with E03 (the "degenerate tag ↔ fall count" link the batch asked for).** The
`degenerate` column of `results/E08_audit_*.csv` and the `class` column of `results/E03_summary.md`
are produced by the *same* function, and `predicted_falls_D4` reproduces all 36 measured counts
exactly. Note that at even n the two most damaging classes ({2,2⁻¹} → 7n and {−2,(−2)⁻¹} → 4n) are
inaccessible because gcd(2,n) ≠ 1, so the recommended even n are structurally protected from them —
worth stating in the paper's design-rule section.

**Next steps.** (b) linear-layer distributions (unit-LU vs uniform GL_n vs same-density random) as an
E04 control — still open, server. (c) a *cheaper attack* inside a degenerate class (needed before the
classification can be called a weakness rather than an observation) — still open.
