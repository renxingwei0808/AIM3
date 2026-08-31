# Claims ledger

Status tags: **proved** · **verified n≤…** (computed, for the parameters actually run) · **conj**
(experimental support, no proof) · **open**. A claim is [VERIFIED] only for the parameters in the
"Evidence" column. Every row links to a CSV or MD file under `results/` that carries
`tier, host, commit, seconds, censored`.

| ID | Claim | Status | Evidence | Verified by |
|---|---|---|---|---|
| A | S_e(x)=(x+1)^{2^e}/x; projective normal form s^{2^e}+ys+y=0 | proved | E01 | |
| B | Im S_e = {0} ∪ {Tr(y^{-r})=0}, size 2^{n-1} | proved | E01 | |
| C | two-valued inverse G_{e,ε} = a(H_e(1/a)+ε)+1 | proved | E01, `src/aim3/structures.py` | |
| D | kernel variable, B_e alternating, B_e(u,v)=uv(M_e(u)+M_e(v)), S_e = M_e∘κ_e | proved | E01 | |
| E | wt(r) = e^{-1} mod n, supp(r) arithmetic progression | proved | E01/E08 | |
| F | ⟨P1,P2,P3⟩ = ⟨P1⟩ (Boolean), S_quad = XL_1(P1) | proved | Focus Report §2.6 | |
| G | syz_3 = 8n per S-box | **verified n=11..16, every admissible e (36 runs)** (3n proved) | `results/E03_sbox_degree_fall.csv`, `results/E03_summary.md` | batch 1 |
| H | first-fall degree 4; falls_D4 = 2n per S-box for generic e | **verified n=11..16, every admissible e** | `results/E03_summary.md` | batch 1 |
| H2 | no cubic generators outside ⟨P1,P2,P3⟩ (`outside_ideal` = 0) | **verified n=11..16 for all e except e=2** | `results/E03_summary.md` | batch 1 |
| H3 | for e = 2 exactly 2n cubics of I(V_e) are *not* reached at D=4, while for e = 2^{-1} none are | verified n=11,13,15 | `results/E03_summary.md` | batch 1 |
| I | falls_D4 depends only on {e, e^{-1} mod n}; it is 7n on {2,2⁻¹}, 4n on {−2,(−2)⁻¹}, 17n/2 on {1+n/2}, 2n+1 on {−3,(−3)⁻¹} and {n/2−1}, and 2n otherwise | **conj** (reproduces all 36 measured counts, n=11..16; untested for n>16) | `results/E03_summary.md`, `aim3.structures.predicted_falls_D4` | batch 1 |
| I2 | the extra fall of the 2n+1 class is produced by 3e ≡ **−1** (mod n); 3e ≡ +1 gives the generic 2n | conj (n=11..16) | `results/E03_summary.md` | batch 1 |
| I3 | all ten recommended AIM3 exponents lie in the generic class (predicted falls_D4 = 2n) | verified (closed form) | `results/E08_audit_{128,192,256}.csv`, `tests/test_core.py` | batch 1 |
| J | Table 4 = generic series + C(N-k+d,d)^2; Boolean series gives 218.4/312.5/529.1 | verified | `results/E02_baselines.csv` | |
| K | no cubic cross relation on (t0,t1) | **verified n=14 and n=16** | `results/E12_cross_relations.csv`, `results/E12_cross_relations_explained.csv` | batch 1 |
| K2 | the cubic relation space on (t0,t1) consists exactly of the univariate image criteria Tr(t_i^{-r_i}); its dimension is #{i : wt(−r_i) ≤ 3} | verified n=14, 16 | `results/E12_cross_relations_explained.csv` | batch 1 |
| L | kernel bilinear model at D=3: 0 syzygies and dim(span ∩ deg≤2) = 2·#eqs = 12n, of which 6n is the input equations themselves, so the pass produces **6n new quadratics** (= #eqs, one per equation) | **verified n=7, 9** (n=5 gives 12n−1 / 6n−1 with 1 syzygy) | `results/E06_kernel_model_falls.csv` | batch 1, wording corrected batch 2 |
| L2 | the kernel model produces new falls at D=4 (937 at n=5, 1924 at n=7), in a non-saturated regime | verified n=5, 7 | `results/E06_kernel_model_falls.csv` | batch 1 |
| L3 | a random bilinear system of the same shape (6n eqs, blocks 3n/3n, planted root) produces **zero new** elements at D=3 and D=4 (its degree-≤2 part is exactly its 6n equations) — so L/L2 are structural, not generic to bilinearity | verified n=5, 7, 9 | `results/E06_kernel_model_falls.csv` | batch 1, wording corrected batch 2 |
| M | syz_3 of the full system is the direct sum: syz_3(S_quad) = (ℓ+1)·8n | **verified n=7,9,11 for ℓ=2 and n=9,11 for ℓ=1** (n=7, ℓ=1 is saturated) | `results/E04_full_system_falls.csv` | batch 1 |
| M2 | the D=4 falls of the first-layer two-S-box system are the exact direct sum of the per-S-box falls (`cross_excess` = 0) | **verified n=11** (generic exponents, non-saturated) | `results/E12_first_layer_ideal.csv` | batch 1 |
| M3 | ideal-theoretic form of M2: dim I(V)_{≤3}(pt,t0,t1) − rank Mac_3 = 4n exactly | **verified n=14** (3444 − 3388 = 56 = 4n) | `results/E12_first_layer_ideal.csv` | batch 1 |
| M4 | the full 9n-equation system's D=4 fall count is **not measurable at Tier A**: the Macaulay row space saturates the degree-≤3 space (0.997–0.9997 of C(nv,≤3)) for n ≤ 11, and the same-shape random control saturates identically | verified n=7,9,11, ℓ=1,2 | `results/E04_full_system_falls.csv` | batch 1 |
| N | the two-valued-inverse reduced system exists explicitly over F_{2^n}[t0,t1] with **2** branch bits (not 3) and the true key is a root | **verified n=8..14, 2 seeds** | `results/E05_reduced_system_stats.csv` | batch 1 |
| N2 | Eq1 has ≈ n² terms, Eq2 ≈ n³; deg_{t0} = deg_{t1} = 2^n−2 for Eq1, so the Bézout bound is ≈ 2^{2n} | verified n=8..14 | `results/E05_reduced_system_stats.csv` | batch 1 |
| N3 | the ROADMAP §3.3 form of the second equation cannot be expanded sparsely (non-linearised ∘ non-linearised, ≈ 2^n terms); only the mixed "outer-inverse / inner-forward" form is sparse | proved (by inspection of the exponents) + verified pointwise n=8..14 | `experiments/E05_two_valued_inverse/reduced_system_toy.py` | batch 1 |
| N4 | the sparse reduced system has 1–3 roots in (F*)² per branch (true key + 0–2 spurious), so pseudo-solution filtering is cheap | verified n=8..14, 2 seeds | `results/E05_reduced_system_stats.csv` | batch 1 |
| O | AIM3 encoding conventions are `lsb / A / col / LU`; the reference implementation's KAT is reproduced for AIM3-I/III/V | **verified** (records 0–2 of each parameter set) | `results/E00_kat_check.md`, `experiments/E00_reference_and_facts/kat_vectors.json` | batch 1 |
| P | the linear-layer distribution is irrelevant to the algebra: unit-LU, uniform GL_n(F_2) and same-density random invertible matrices give the same D=3 rank and syz_3 and the same D=4 fall count for the first-layer system | **verified n=9, 11, two seeds** (at n=11 all three give rank 2068, syz_3 = 176 = 2·8n, falls_D4 = 44 = 4n) | `results/E01_theory_matrix_distribution.csv` | batch 2 |
| Q | reused key (same pt under q different iv): the D=3 syzygy count exceeds the direct sum of q independent instances by exactly **C(q,2)·n** — one n-block per pair of instances | **conj**, verified in the non-saturated runs (n,ℓ,q) = (9,2,2), (9,2,3), (9,2,4), (7,2,2), (7,2,3), (9,1,2); fails in the two most saturated runs, which are flagged | `results/E09_theory_reused_key.csv` | batch 2 |
| R | the natural ansatz for the D=4 fall generators — the F_{2^n}-module spanned by x^{2^i}y^{2^j}P1, x^{2^k}P2, y^{2^k}P3 — explains **none** of the 2n new cubics; every one of its falls already lies in Mac_3, and no single one of the 18 multiplier families {1, x^{2^i}, y^{2^j}, x^{2^i}y^{2^j}, x^{2^i+2^j}, y^{2^i+2^j}} × {P1,P2,P3} yields a new cubic on its own | **verified (n,e) = (11,3), (13,3)** — a refutation, not a construction | `results/E03_theory_fall_generators.csv`, `docs/THEORY_NOTES.md` | batch 2 |
| S | one forward AIM3 evaluation costs 2^20.85 / 2^22.30 / 2^23.45 bit operations (schoolbook) or 2^19.90 / 2^21.14 / 2^22.20 (Karatsuba), from the addition chains of spec Table 5; the retired placeholder 2^14 understated this by 6–9 bits | verified (arithmetic from the spec table, cross-checked by an independent Itoh–Tsujii derivation: 50/63/77 vs the spec's 47/59/73 multiplications) | `results/E11_theory_eval_cost.csv`, `docs/COMPLEXITY_ACCOUNTING.md` §1 | batch 2 |

