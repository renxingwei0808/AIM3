# E05 — Two-valued-inverse reduced-variable systems and sparse elimination (track T3)

**Claim being tested.** The spec says reduced-variable systems "no longer exist" because the S-box is not invertible. The
two-valued inverse G_{e,ε}(y) = y^r·H_e(y^{-r}) + ε·y^r + 1 (ε ∈ F₂) gives, for AIM3-I/III, a system of **2 field equations in
2 field unknowns (t₀,t₁) plus 3 branch bits**. The question is only its *cost*.

**Facts established (`G_stats.py`).** For every recommended (n,e): G has n+2 monomials in F_{2^n}[y]; its n H-terms lie in **n
distinct cyclotomic cosets**, with Boolean degrees from 0 up to n−1, and ordinary degree ≈ 2^n. So the reduced system is
field-sparse but Boolean-dense; the AIM2 low-degree phenomenon (x = y^{2^e−1}) is gone.

**Plan.**
1. Tier A: build the reduced system symbolically for toy n (n = 8..16) using `SboxStructures.H` for the section; count terms after
   composing with A₁∘A₀⁻¹ (≤ n(n+2) terms), verify the true key is a root for the right branch (ε₀,ε₁,ε₂), count spurious roots.
2. Tier C (Magma, `magma/reduced_system_resultant.m`): Resultant / Kronecker / sparse-resultant timings for n = 8..24, fit
   log₂(time) vs n; test polynomial-decomposition (Yang–Zheng–Yang) on G.
3. Compare with the AIM2-style reduced system of the same n to state "non-permutation patch = X bits".
4. Negative result to prove: any univariate polynomial obtained by eliminating t₁ has degree ≥ 2^{2n−O(1)} unless coset
   collapse occurs — root finding alone then costs ≥ 2^{2n}, i.e. this route cannot beat B1 for AIM3-I (2^{256} vs 2^{262.6} is
   within the model's error). State this carefully as a *bound for this route*, not a security proof.
