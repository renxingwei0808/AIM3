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

---

## Log — 2026-08-31 (Tier A, batch 1)

**Status.** D=3 direct sum [VERIFIED n=7,9,11, ℓ=1,2]. **D=4 direct sum: NOT decidable at Tier A —
every full-system run saturates.** Gate 1 is answered instead by E12 (see below).

**What was run.** `full_system_falls.py` (new) via `bash run_tierA.sh`: n = 7, 9, 11; ℓ = 2 (three
S-boxes, 3n variables) and ℓ = 1 (two S-boxes, 2n variables); seeds 1 and 2; systems `S_quad`,
`P1only` and a same-shape random quadratic control with a planted root; `--per-block` control (each
S-box's own 3n quadratics alone in the same variables). → `results/E04_full_system_falls.csv`.

**Result 1 — D=3 is an exact direct sum.** syz₃(S_quad) = (ℓ+1)·8n exactly, on every non-saturated
run: 168 = 3·8·7, 216 = 3·8·9, 264 = 3·8·11 (ℓ=2) and 144 = 2·8·9, 176 = 2·8·11 (ℓ=1). The random
control has 0 syzygies at the same shape. (ℓ=1, n=7 is already saturated at D=3 and is excluded.)

**Result 2 — D=4 saturates, so `falls_D4` of the full system measures nothing here.** S_quad has
3n(ℓ+1) equations in only (ℓ+1)n variables, so at D=4 its Macaulay row space already fills the
degree-≤3 space: `dim(span ∩ deg≤3) / C(nv,≤3)` is 0.997–0.9997 in every run, and the **same-shape
random control saturates identically** (e.g. n=11, ℓ=2: S_quad 6016/6018, random 6017/6018). The
apparent "cross-falls" verdict of the naive (ℓ+1)·2n criterion is therefore an artefact of the
monomial basis, not evidence of cross-S-box structure. Runs in this state are recorded with
`verdict = indeterminate(saturated)` and must not be read either way. The first non-saturated full
system would be n ≥ 13 at D=4 (rows/cols ≈ 0.99 at n=13, 0.85 at n=15) — Tier B/C.

**Result 3 — the naive criterion is wrong even without saturation.** Comparing the joint fall count
with (ℓ+1)·2n is not a direct-sum test: the joint Macaulay matrix also multiplies by cross-block
monomials (pt·t₁, t₀·t₁, …) that do not exist in E03's 2n-variable single-S-box ring, so a single
S-box embedded in (ℓ+1)n variables already gives far more than 2n falls (measured: 45–77 at n=7,9).
The `--per-block` control replaces it with `cross_excess = falls(all) − Σᵢ falls(S-box i alone)`.

**Gate 1 — answered by E12, not by E04.** The first-layer system (6n quadratics in 3n variables,
`experiments/E12_cross_sbox_relations/first_layer_ideal.py`) is *not* saturated and gives
`cross_excess = 0` at n=11 and an ideal-theoretic excess of exactly 4n at n=14. See that README.

**Next steps (server).**
* Tier B/C: n = 13, 15 at D=4 with M4RI (not saturated there) and the `--per-block` control;
  ℓ = 3 for the AIM3-V shape; uniform GL_n as a second matrix distribution.
* Tier C: the actual d_solv via F4 step degrees (`magma/`, `SetVerbose("Faugere",2)`) — the
  saturation problem does not affect a step-degree log, which is the better observable for Q2.
