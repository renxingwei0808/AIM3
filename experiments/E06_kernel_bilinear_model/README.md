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
