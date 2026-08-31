# E06 — Kernel-variable bilinear model (track T2)

**Model (Focus Report §4.3).** Variables U = (pt, t₀, t₁), W = (a₀, a₁, a₂) (kernel elements a_i = t_i^{r_i}); 6n equations,
all F₂-bilinear in (U, W): B_{e_i}(a_i, A_i pt + b_i) = a_i, t_i a_i = a_i^{2^{e_i}}, and the same for the last S-box with
x₂ = b₂ + B₀t₀ + B₁t₁, t₂ = pt + ct. Branches are absorbed (B_e alternating ⇒ both preimages x, x+a satisfy the same equation).

**Important framing (do not overclaim).** In *field* variables this is a reparametrisation (still ℓ+1 field unknowns), not a
variable reduction. Its value is the *shape*: every equation is Frobenius-linear in one block and scalar-affine in the other
("twisted bilinear"), which is the input format of block-F4/F5, y-XL, MinRank/Support-Minors and Dickson-matrix methods.

**First observation (n=7, `kernel_bilinear_toy.py 7 3`).** 42 equations / 42 variables: no syzygies at D=3 and
**dim(span ∩ deg≤2) = 84 = 2·#eqs**, versus S_quad (168 syzygies, 175) and P1-only (0, 63).
*Careful with the wording* (corrected in batch 2): 84 is the whole degree-≤2 part of the row space, and 42 of
that is the 42 original bilinear equations themselves. The number of **new** quadratics produced at D=3 is
84 − 42 = 42 = #eqs = 6n, not 12n.

**Plan.** (1) Tier A: falls vs n for n = 5..9, D ≤ 4; identify the fall generators (likely t_i·B-equation vs a_i·K-equation
cancellations). (2) Tier C: F4 with block orders {U > W} and {W > U}; y-XL multiplying only by W-monomials; compare solving
degrees with S_quad on identical instances. (3) Write the U-linear form N(W)·(U,1)ᵀ = 0 and test Support-Minors / KS modelling.
(4) Feed the measured fall structure into E04. Success = a model whose d_solv grows slower than S_quad's with n.

---

## Log — 2026-08-31 (Tier A, batch 1)

**Status.** At D=3 the kernel model has dim(span ∩ deg≤2) = 2·#eqs = 12n [VERIFIED n=7,9] (12n−1 at n=5).
Since the 6n equations are themselves of degree 2 and lie in the row space, the **new** quadratics produced
at D=3 number 12n − 6n = **#eqs = 6n** [VERIFIED n=7,9] (6n−1 at n=5, together with one syzygy).
New falls do appear at D=4 [VERIFIED n=5,7]. The **random bilinear control of the same shape produces no
new elements at either degree**, so the behaviour is structural, not a property of bilinearity.

**What was run.** `kernel_bilinear_toy.py n Dmax` with the new `random_bilinear_system(n)` control
(6n equations, blocks U = 3n / W = 3n, random U×W bilinear + affine, planted root) and a random
quadratic control of the S_quad shape: n = 5 and 7 at D ≤ 4, n = 9 at D ≤ 3 →
`results/E06_kernel_model_falls.csv`. `dim_low` below is dim(span ∩ deg ≤ D−1).

`dim_low` = dim(span ∩ deg≤D−1); `new` = dim_low(D) − rank(D−1), i.e. the elements the pass actually adds
(this is the `falls_from_prev_D` column of the CSV).

| system | n | D=3 syz | D=3 dim_low | 2·#eqs | D=3 **new** (=#eqs?) | D=4 syz | D=4 dim_low | D=4 new |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| kernel (6n/6n) | 5 | 1 | 59 | 60 | 29 (=6n−1) | 598 | 1866 | 937 |
| kernel (6n/6n) | 7 | 0 | **84** | 84 | **42 = 6n** | 974 | 3730 | 1924 |
| kernel (6n/6n) | 9 | 0 | **108** | 108 | **54 = 6n** | — | — | — |
| random bilinear (6n/6n) | 5 | 0 | 30 | 60 | **0** | 465 | 930 | **0** |
| random bilinear (6n/6n) | 7 | 0 | 42 | 84 | **0** | 903 | 1806 | **0** |
| random bilinear (6n/6n) | 9 | 0 | 54 | 108 | **0** | — | — | — |
| S_quad (9n/3n) | 7 | 168 | 175 | — | 112 | 7696 | 1558 | 340 |
| P1-only (3n/3n) | 7 | 0 | 63 | — | 42 | 273 | 1218 | 756 |

**Answers to the three questions.**
1. *Is dim(span ∩ deg≤2) at D=3 always 2·#eqs = 12n?* Yes at **n = 7 and n = 9** (84 and 108, exactly).
   At n = 5 it is 59 = 12n−1, together with a single syzygy — so the identity is [VERIFIED n=7,9]
   and **fails by one at n=5**; it is not stated for n ≤ 5. Equivalently, the D=3 pass produces
   **6n new quadratics** (= one per equation), not 12n: half of the degree-≤2 part of the row space is
   the 6n input equations themselves.
2. *Does D=4 produce new falls?* Yes, and a lot: dim_low jumps to 1866 (n=5) and 3730 (n=7), i.e.
   937 and 1924 new elements beyond rank_{D=3}. These runs are **not** saturated (dim_low/C(6n,≤3)
   = 0.41 and 0.30), unlike every full-system D=4 run in E04.
3. *What does the random bilinear control do?* Nothing: at D=3 the degree-≤2 part of its row space is
   exactly the 6n input equations (so **0 new** quadratics, against the kernel model's 6n), and at D=4
   `dim_low` equals rank_{D=3} exactly — **zero new elements at both degrees, for n = 5, 7, 9**. The
   kernel model's falls are therefore not a generic consequence of the bipartite bilinear shape.

**Next steps.** Tier B/C: n = 11..16 at D = 3, 4 with M4RI; identify the fall generators (candidate:
tᵢ·(B-equation) vs aᵢ·(K-equation) cancellations); block orders {U > W} and {W > U} in Magma F4;
y-XL multiplying only by W-monomials. Feed the fall structure into E04.
