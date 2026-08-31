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

---

## Log — 2026-08-31 (Tier A, batch 1)

**Status.** The reduced-variable system is now **explicit and machine-checked** for n = 8..14
[VERIFIED]. Its ordinary degrees confirm the premise of the T3.2 negative result.

**What was run.** `reduced_system_toy.py` (new), n = 8..14, seeds 1 and 2 →
`results/E05_reduced_system_stats.csv`. Steps 1–4 of the plan above are done; the resultant (step 2)
is left to Tier C.

**Result 1 — H_e as a linearised polynomial.** `SboxStructures.H` is F₂-linear on the trace-zero
hyperplane; extending it by H(γ) := 0 on one trace-one basis vector γ pins the single degree of
freedom (linearised polynomials vanishing on the hyperplane are exactly c·Tr(u)) and makes the n×n
Moore system solvable. The resulting Σⱼ hⱼ u^{2^j} agrees with `H(u)` on 100 random trace-zero u for
every n and both seeds.

**Result 2 — G_{e,ε} as a sparse dictionary.** {G_{e,0}(t), G_{e,1}(t)} = {x, x+a} for 100 random
x ≠ 1, every n, both seeds.

**Result 3 — the shape of the reduced system.** With x₂ = b₂+B₀t₀+B₁t₁ and t₂ = A₀⁻¹(G_{e₀,ε₀}(t₀)+b₀)+ct:

| | Eq1 | Eq2 |
|---|---|---|
| form | (A₁A₀⁻¹)(G_{e₀,ε₀}(t₀)+b₀) + b₁ + G_{e₁,ε₁}(t₁) | x₂·t₂ + x₂^{2^{e₂}} + 1 |
| terms | ≈ n² (measured 64, 80, 100, 121, 145, 169, 197 for n = 8..14) | ≈ 1.3–1.6 n³ (678 … 4446) |
| deg in t₀ | 2ⁿ−2 | 2ⁿ−2 |
| deg in t₁ | 2ⁿ−2 | 2^{n−1} |
| Boolean degrees | ≈ n+1 terms at each degree 0..n−1 | spread over 0..n |

Bézout bound deg_{t₀}·deg_{t₁} = (2ⁿ−2)² ≈ **2^{2n}** for Eq1 (measured 64516 at n=8 vs 2^16 = 65536,
… 268 369 924 at n=14 vs 2^28). This is exactly the input the T3.2 lower-bound argument needs.

**Result 4 — a correction to the ROADMAP form of the system, and one fewer branch bit.**
The ROADMAP §3.3 second equation, `b₂+B₀t₀+B₁t₁ = G_{e₂,ε₂}(A₀⁻¹(G_{e₀,ε₀}(t₀)+b₀)+ct)`, **cannot be
expanded into a sparse polynomial**: it composes a non-linearised sparse polynomial (G_{e₂,ε₂}, whose
exponents are ≈ 2ⁿ) with another non-linearised sparse polynomial, and the expansion has ≈ 2ⁿ terms.
Writing the last S-box forward instead — through its own P1 relation x₂·t₂ = x₂^{2^{e₂}}+1 — keeps
everything sparse (only "linearised ∘ sparse", sparse products and additions are needed) and removes
ε₂ entirely, because the forward relation has no branch ambiguity. So the reduced system carries
**2 branch bits, not 3**: 4 branches to try, not 8. The ROADMAP form is still verified *pointwise* at
the true key for the correct (ε₀,ε₁,ε₂) by `verify()`.

**Result 5 — pseudo-solutions.** Exact root counts of the sparse system in (F*)² for the true branch,
computed in O(2ⁿ) by tabulating G_{e₁,ε₁}: 1 to 3 roots (true key + 0 to 2 spurious) across n = 8..14
and both seeds. So filtering is cheap; the cost of this route is the elimination, not the filtering.

**Input format for the server (Tier C).** `reduced_system_toy.py` builds `SPoly` objects whose `.t`
is a dict `{(i, j): c}` — exponent of t₀, exponent of t₁ (both mod 2ⁿ−1), coefficient in F_{2ⁿ} as an
integer in the §4.1.1 bit convention. For Magma, emit each equation as
`&+[ F!c * t0^i * t1^j : ... ]` over `F<w> := GF(2, n)` with the toy modulus from
`src/aim3/gf2n.py::TOY_MODULI`, and run `Resultant(eq1, eq2, t1)` / Kronecker substitution. The
matrices, constants and the true key come from the same `(n, seed)` — regenerate them with
`build_instance(n, exps, seed)` and print them, or dump them alongside. Both branches (ε₀,ε₁) ∈ F₂²
must be run; the true one is printed by the script.

**Next steps.** Tier C: resultant/Kronecker timings for n = 8..24 and the log₂(time)-vs-n slope
(the T3 termination criterion is slope ≥ 2); polynomial decomposition attempts on G; comparison with
the AIM2 reduced system of the same n to state "non-permutation patch = X bits".
