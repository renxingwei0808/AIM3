# E06 — Kernel-variable bilinear model (track T2)

**Model (Focus Report §4.3).** Variables U = (pt, t₀, t₁), W = (a₀, a₁, a₂) (kernel elements a_i = t_i^{r_i}); 6n equations,
all F₂-bilinear in (U, W): B_{e_i}(a_i, A_i pt + b_i) = a_i, t_i a_i = a_i^{2^{e_i}}, and the same for the last S-box with
x₂ = b₂ + B₀t₀ + B₁t₁, t₂ = pt + ct. Branches are absorbed (B_e alternating ⇒ both preimages x, x+a satisfy the same equation).

**Important framing (do not overclaim).** In *field* variables this is a reparametrisation (still ℓ+1 field unknowns), not a
variable reduction. Its value is the *shape*: every equation is Frobenius-linear in one block and scalar-affine in the other
("twisted bilinear"), which is the input format of block-F4/F5, y-XL, MinRank/Support-Minors and Dickson-matrix methods.

**First observation (n=7, `kernel_bilinear_toy.py 7 3`).** 42 equations / 42 variables: no syzygies at D=3 but **84 = 2·#eqs
degree falls to degree ≤2**, versus S_quad (168 syzygies, 175 falls incl. the P2/P3 duplicates) and P1-only (0, 63).

**Plan.** (1) Tier A: falls vs n for n = 5..9, D ≤ 4; identify the fall generators (likely t_i·B-equation vs a_i·K-equation
cancellations). (2) Tier C: F4 with block orders {U > W} and {W > U}; y-XL multiplying only by W-monomials; compare solving
degrees with S_quad on identical instances. (3) Write the U-linear form N(W)·(U,1)ᵀ = 0 and test Support-Minors / KS modelling.
(4) Feed the measured fall structure into E04. Success = a model whose d_solv grows slower than S_quad's with n.

---

## Log — 2026-08-31 (Tier A, batch 1)

**Status.** The D=3 fall count of the kernel model is 2·#eqs = 12n [VERIFIED n=7,9] (n=5 gives
12n−1). New falls do appear at D=4 [VERIFIED n=5,7]. The **random bilinear control of the same shape
has zero falls at both degrees**, so the behaviour is structural, not a property of bilinearity.

**What was run.** `kernel_bilinear_toy.py n Dmax` with the new `random_bilinear_system(n)` control
(6n equations, blocks U = 3n / W = 3n, random U×W bilinear + affine, planted root) and a random
quadratic control of the S_quad shape: n = 5 and 7 at D ≤ 4, n = 9 at D ≤ 3 →
`results/E06_kernel_model_falls.csv`. `dim_low` below is dim(span ∩ deg ≤ D−1).

| system | n | D=3 syz | D=3 dim_low | 2·#eqs | D=4 syz | D=4 dim_low | D=4 falls |
|---|---:|---:|---:|---:|---:|---:|---:|
| kernel (6n/6n) | 5 | 1 | 59 | 60 | 598 | 1866 | 937 |
| kernel (6n/6n) | 7 | 0 | **84** | 84 | 974 | 3730 | 1924 |
| kernel (6n/6n) | 9 | 0 | **108** | 108 | — | — | — |
| random bilinear (6n/6n) | 5 | 0 | 30 | 60 | 465 | 930 | **0** |
| random bilinear (6n/6n) | 7 | 0 | 42 | 84 | 903 | 1806 | **0** |
| random bilinear (6n/6n) | 9 | 0 | 54 | 108 | — | — | — |
| S_quad (9n/3n) | 7 | 168 | 175 | — | 7696 | 1558 | 340 |
| P1-only (3n/3n) | 7 | 0 | 63 | — | 273 | 1218 | 756 |

**Answers to the three questions.**
1. *Is the D=3 fall count always 2·#eqs = 12n?* Yes at **n = 7 and n = 9** (84 and 108, exactly).
   At n = 5 it is 59 = 12n−1, together with a single syzygy — so the identity is [VERIFIED n=7,9]
   and **fails by one at n=5**; it is not stated for n ≤ 5.
2. *Does D=4 produce new falls?* Yes, and a lot: dim_low jumps to 1866 (n=5) and 3730 (n=7), i.e.
   937 and 1924 new elements beyond rank_{D=3}. These runs are **not** saturated (dim_low/C(6n,≤3)
   = 0.41 and 0.30), unlike every full-system D=4 run in E04.
3. *What does the random bilinear control do?* Nothing: at D=3 its row space contains only the 6n
   equations themselves, and at D=4 `dim_low` equals rank_{D=3} exactly — **zero falls at both
   degrees, for n = 5, 7, 9**. The kernel model's falls are therefore not a generic consequence of
   the bipartite bilinear shape.

**Next steps.** Tier B/C: n = 11..16 at D = 3, 4 with M4RI; identify the fall generators (candidate:
tᵢ·(B-equation) vs aᵢ·(K-equation) cancellations); block orders {U > W} and {W > U} in Magma F4;
y-XL multiplying only by W-monomials. Feed the fall structure into E04.

## Server log — 2026-08-31 (server batch 1, branch `claude/server-batch-1`)

`kernel_bilinear_toy.py` gained two options: `--backend python|sage` (the M4RI rank engine, verified
identical to the pure-Python one by `server/jobs/verify_backend.py` at n = 7 for D = 2, 3, 4) and
`--export DIR`, which writes the kernel model, `S_quad`, `P1only` and both random controls in msolve
and Magma syntax for the F4 step-degree runs of E04 task 6.

Sizing, measured rather than guessed: the D = 4 Macaulay matrix of the kernel model is
145992 × 768692 at n = 11 (≈ 14 GB dense) and 240396 × 1505583 at n = 13 (≈ 45 GB) — the n = 15
matrix would be ≈ 123 GB *and* about a day of M4RI, so it is out of reach for this batch and is
listed as not done.

Rows: `results/E06_server_kernel.csv`.  D ≤ 4 was run at n = 11; at n = 13 and n = 15 only D ≤ 3,
because the D = 4 matrices are 45 GB / 123 GB and 6 h / ~1 day of M4RI — out of proportion to this
batch and, at n = 15, above the one-Tier-C-job cap.  The D = 3 result `dim_low = 2·#eqs` continues
to hold, and the same-shape random bilinear control continues to show none of it.
Numbers and the comparison with the random bilinear control:
`results/REPORT_server_batch1.md` §4.
