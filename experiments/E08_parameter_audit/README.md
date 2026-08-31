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
