# Paper outline — working title
**Twisted-Bilinear Structure of AIM3: Kernel Coordinates, Degree Falls, and the Cost of Non-Invertibility**
(alternative if T3/T4 succeed: *New Elimination Techniques for AIM3: Two-Valued Inversion and Skew-Polynomial Elimination*)

Target venues: IACR ToSC/FSE (structural + attack), ASIACRYPT/EUROCRYPT (if a concrete complexity below B1 with a new framework), DCC (if theory-heavy).

## Abstract (draft skeleton)
AIM3 replaces AIM2's invertible Mersenne S-boxes by 2-to-1 "enhanced" S-boxes and adds an iv-dependent input affine layer;
its algebraic security argument rests on (i) the non-existence of reduced-variable systems and (ii) a semi-regular model of
3(ℓ+1)n quadratics. We give the enhanced S-box a closed form and a kernel-coordinate description: y = S_e(x) iff a = y^r
satisfies the alternating bilinear identity B_e(a,x) = a, so that S_e = M_e ∘ κ_e and every AIM3 relation is Frobenius-linear
in one block of variables and scalar-affine in the other ("twisted bilinear"). From this we derive: exact syzygy and
degree-fall counts (8n at degree 3, 2n at degree 4 per S-box), a closed-form weight theorem for the inverse exponent, a
constant-branch two-valued inverse giving reduced-variable systems, and [new elimination technique / attack with complexity …].
[Numbers, comparison with the designers' 2^{262.6}/2^{381.0}/2^{651.5}, and what remains above 2^λ.]

## Sections and supporting evidence
1. Introduction — positioning vs AIM/AIM2 cryptanalysis; contributions list (Base / Target / Stretch, see below).
2. Preliminaries — AIM3 (fact sheet §1), Boolean ideals, bilinear systems, Ore rings.
3. The enhanced Mersenne S-box: closed form, image, kernel coordinates, B_e(u,v)=uv(M_e(u)+M_e(v)), weight theorem. [E01, E08]
4. Ideal structure of S_quad: ⟨P1,P2,P3⟩=⟨P1⟩, S_quad = XL₁(P1), syz₃=8n, first-fall degree 4 with 2n falls, exceptional exponent classes. [E03, E12]
5. Twisted-bilinear models of AIM3: P1-system as (X|T) twisted bilinear; kernel model (6n); Dickson-matrix / MinRank view; structured Hilbert function and measured solving degrees vs the semi-regular baselines B1/B2. [E04, E06, E10, E11, E02]
6. Reduced-variable systems via the two-valued inverse; sparse elimination; the cost of non-invertibility (= X bits vs AIM2). [E05]
7. Skew-polynomial elimination (Magic Pot over F_{2^n}⟨φ⟩): framework, degree-growth theorem (positive or negative). [E07]
8. Parameter audit and multi-target considerations. [E08, E09]
9. Complexity summary table (B1/B2/B3 vs each technique, ω=2 and 2.81) and conclusions.

## Contribution tiers
* **Base (already supported):** §3, §4, E02 calibration, E12 cubic non-existence, E06 model + falls, audit criteria.
* **Target:** measured/extrapolated d_solv of the twisted-bilinear models beating B1 on ≥1 parameter set; quantified cost of the two-valued-inverse route; degree-fall closed form as a theorem.
* **Stretch:** skew elimination with polynomial row-degree growth ⇒ attack below B2; or a weak-parameter family with a concrete attack.
