# E07 — Ore / skew-polynomial elimination (Magic-Pot generalisation, track T4)

**Idea.** Over F_{2^n}, each AIM3 relation is *Frobenius-linear* in one block and *scalar-affine* in the other:
first-layer  Eq_i: (φ^{e_i}∘A_i + t_i·A_i)(pt) = c_i(t_i)   — an element Λ_i(t_i) of the Ore ring R = F_{2^n}[t_i]⟨φ⟩ (φa = a²φ)
acting on pt;  last layer: Λ₂(pt)(t₀,t₁) = c₂(pt).  Magic Pot's F_{2^n}[pt]-linear algebra becomes module computations over R.

**Core questions.**
1. Common-root elimination of pt between Λ₀ and Λ₁ via LCLM / GCRD (Ore's Euclidean algorithm) and the *skew resultant*: what
   is the degree, in the conjugates t₀^{2^k}, t₁^{2^k}, of the resulting scalar condition? (Conjecture: ≈ n in each — the Dickson
   determinant of x^{2^e}+t·x is a norm-like polynomial; then Frobenius consistency must be imposed.)
2. Row-degree growth of weak-Popov / approximant bases over R with coefficients in F_{2^n}[t]: polynomial or exponential in n?
3. Is there a scalar s ∈ F_{2^n} (candidates: t₂ = pt+ct, a₂, or a normal-basis coordinate) in which *all* equations are
   polynomial of moderate degree? (Obstruction: every choice meets Frobenius closure somewhere; a negative theorem is publishable.)

**Tools.** Sage `OrePolynomialRing(K, K.frobenius_endomorphism(), 'phi')` (skeleton in `sage/skew_probe.sage`), Magma
`TwistedPolynomialRing` / `SkewPolynomialRing` if available on the server version. Tier B/C.
**Kill criterion.** If for n = 8..20 the LCLM/approximant row degrees grow like 2^{cn} with c ≥ 1/2, close the attack line and
write the impossibility argument instead (still a T4 deliverable).
