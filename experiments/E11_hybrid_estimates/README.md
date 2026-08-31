# E11 — Generic hybrid estimates for the bilinear view (track T2, accounting)

`hybrid_estimates.py` prints, for k guessed pt-bits, (a) the FSS bound d_reg ≤ min(n_X,n_T)+2 for affine bilinear systems and
(b) the y-XL linearisation degree. Conclusion: with *generic* bilinear bounds no hybrid beats exhaustive search (the curves
only fall below B1 = 262.6 when k ≥ ~118, i.e. next to brute force). Any T2 attack must therefore exploit **non-generic**
structure (falls, syzygies, twisted-bilinear rank defects) — this script quantifies the bar.
