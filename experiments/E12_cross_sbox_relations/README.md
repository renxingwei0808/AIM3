# E12 — Cross-S-box relations of degree ≥ 3

The spec proves no *quadratic* relation exists between two S-box outputs t_i, t_j. `cross_relations.py 14 3` shows: **no cubic
relation on (t₀,t₁) either** (dim 0 of 3683 monomials on 16382 points). On (pt,t₀,t₁) the cubic relation space has dimension
3444 at n=14; the next step is to compare this with the degree-3 Macaulay rank of the union of the two S-box ideals in 3n
variables (extend `legacy/aim3_macaulay.py` to a first-layer-only system) — any excess = cross degree falls, feeding E04.
Also run deg=4 on (t₀,t₁) at n=16 (needs 2^n > C(32,≤4)=41449 ✓; use M4RI on the server for the 65k×41k rank).
