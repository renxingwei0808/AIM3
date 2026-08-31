# E04 — Structured Hilbert function and true solving degree of the full AIM3 system

**Goal (track T1 core).** Replace "9n random quadratics" by the actual bigraded structure: (ℓ+1) S-box ideals, each with
syz₃ = 8n and 2n falls at D=4 (E03), glued by the linear layers. Determine, for toy instances, the solving degree d_solv of
`S_quad`, `P1-only`, and the kernel model (E06), and whether cross-S-box falls appear at D ≥ 4.

**Questions.**
* Q1 Is the D=4 fall count of the full system exactly (ℓ+1)·2n (direct sum, like syz₃ = (ℓ+1)·8n at D=3) or larger (cross falls)?
* Q2 d_solv(S_quad) vs d_solv(P1-only) vs random controls, for n = 7..16, ℓ = 1,2,3, several seeds and both unit-LU and uniform GL_n.
* Q3 Does a block order (pt-first or t-first) change the degree profile (link to E10)?

**Method.** Tier A: `legacy/aim3_macaulay.py n D` and `legacy/aim3_solvdeg.py n D l` (dense Macaulay, n≤11, D≤4).
Tier C: Magma F4 with `SetVerbose("Faugere",2)` to log the step degrees; msolve `-v 2` similarly. Skeleton in `magma/`.
Record per run: n, ℓ, exps, seed, matrix distribution, system, order, max step degree, #new polynomials per degree, time, peak RAM.

**Deliverable.** A table "d_solv vs n" for each system and a fitted model; the honest extrapolation rule (docs/COMPLEXITY_ACCOUNTING.md §4)
decides whether a closed-form bigraded Hilbert series can be claimed. This table is the new algebraic baseline of the paper.
