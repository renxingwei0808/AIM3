# E09 — Multi-target, weak-iv and reused-key scenarios (track T6)

AIM3 explicitly targets multi-target one-wayness (spec §3.1, §6.1). Questions:
1. **Weak iv classes.** Is there a property Π of the generated (A_j, b_j) (e.g. A₁A₀⁻¹ Frobenius-sparse, common invariant
   subspace, low-rank part of the twisted-bilinear operator) that makes any track-T2/T3 attack cheaper by a factor 2^s, with
   density ρ(Π) among iv's? Amortised cost over N targets: min over Π of 2^{-s}·cost / (N·ρ(Π)) once N·ρ ≥ 1.
2. **Reused key with different iv** (misuse model, as in the Magic Pot paper): q pairs (iv_j, ct_j) with the same pt give
   3nq equations in n + 2nq variables (P1). Compute d_solv for q = 2,3 on toy n (E04 machinery) — does the extra pt-coupling help?
3. **Zero-input rejection leakage**: valid keys avoid ℓ+1 affine hyperplanes' points; quantify (negligible, but state it).
Tier A/B. Deliverable: a scenario table with amortised formulas; results in the single-key and multi-target models are reported separately.
