# AIM3 fact sheet (every entry carries a spec page reference; AIMer v3.0 English specification, 59 pages)

Status legend: [SPEC] quoted from the specification · [VERIFIED] independently checked in this repo · [PENDING] to check.

## 0. Encoding conventions
| # | Fact | Ref | Status |
|---|---|---|---|
| 0.1 | Field element ↔ bitstring: a₀‖…‖a_{n−1} ↔ Σ aᵢXⁱ ↔ (a₀,…,a_{n−1})ᵀ; hex example A0 00…01 ↔ X¹²⁷+X¹²⁵+1 | p.16 §4.1.1 | [SPEC] |
| 0.2 | Tape bit order inside SHAKE bytes; A·v vs Aᵀ·v; row- vs column-major tape matrix; A=L·U vs U·L | p.19–20 Alg.4 | **[VERIFIED E00]** |
| 0.2a | **Pinned convention: `lsb / A / col / LU`** — bit k of the tape is bit (k mod 8) of byte k/8; each squeezed field element is one **column** of the tape matrix; L = strictlower+I, U = strictupper+I; A = L·U; applied as A·pt. (`lsb / AT / row / LU` is the same matrix.) | KAT + `Reference_Implementation/aim3.c` | [VERIFIED E00] |
| 0.3 | Key layout: `pk = iv‖ct`, `sk = pt‖iv‖ct`; a field element is little-endian bytes | `sign.c::crypto_sign_keypair`, `field_common.c::gf_from_bytes` | [VERIFIED E00] |
| 0.4 | Reproduces the official KAT for AIM3-I/III/V, records 0–2 of each (`results/E00_kat_check.md`, `experiments/E00_reference_and_facts/kat_vectors.json`, regression test `test_kat_conventions_pinned`) | KAT `PQCsignKAT_{48,72,96}.rsp` | [VERIFIED E00] |

## 1. Definition
| # | Fact | Ref | Status |
|---|---|---|---|
| 1.1 | AIM3(pt,iv) = S[e_ℓ] ∘ Lin[iv] ∘ S[e₀..e_{ℓ−1}] ∘ Lin_in[iv](pt) ⊕ pt | p.12 §3.1 | [SPEC][VERIFIED impl] |
| 1.2 | S[e](x) = x^{2^e−1} + x^{2^n−2}, gcd(e,n)=1, S(0)=0 | p.13 | [VERIFIED] |
| 1.3 | Lin_in[iv](x) = (A₀x+b₀)‖…‖(A_{ℓ−1}x+b_{ℓ−1}); Lin[iv](x₀..) = b_ℓ + Σ A_{ℓ+i}(xᵢ) | p.13 | [VERIFIED impl] |
| 1.4 | Moduli: X¹²⁸+X⁷+X²+X+1, X¹⁹²+X⁷+X²+X+1, X²⁵⁶+X¹⁰+X⁵+X²+1 | p.13 | [SPEC] |
| 1.5 | Table 2: I (128; 3,7,11), III (192; 11,23,47), V (256; 3,7,11,19) | p.14 | [SPEC] |
| 1.6 | GenerateLinear: tape = SHAKE_λ(iv, 2ℓn²+(ℓ+1)n); L_j, U_j unit-triangular, A_j = L_j·U_j; b_j from the tail of the tape | p.20 Alg.4 | [VERIFIED impl + KAT, E00] |
| 1.7 | NoZeroInpSB outputs ⊥ if any S-box input is 0; keygen resamples (pt,iv) | p.21 Alg.5, p.20 | [VERIFIED impl] |
| 1.8 | SHAKE128 for λ=128, SHAKE256 for 192/256 | p.17 §4.1.2 | [SPEC] |

## 2. Official security claims about AIM3 that this project addresses
| # | Claim | Ref | Our status |
|---|---|---|---|
| C1 | dim QE(V) = 3n for every recommended S-box (P1,P2,P3; lower bound by trace argument, upper bound by computer) | p.37–39 | [VERIFIED] agree |
| C2 | "the reduced-variable systems of higher degree available for AIM2 no longer exist, since building them requires inverting the S-box" | p.39 | reduced systems **exist** and are now explicit and machine-checked for n=8..14 (E05): 2 field unknowns and **2** branch bits (not 3 — the last S-box is written forward), Eq1 ≈ n² terms, Eq2 ≈ n³ terms, deg_{t₀}·deg_{t₁} ≈ 2^{2n}; cost still to be quantified on the server |
| C3 | S_quad = 3(ℓ+1)n quadratics in (ℓ+1)n variables; Table 4: k=57/102/150, d_reg=19/25/44, log₂T=262.6/381.0/651.5 (ω=2) | p.39–40 | reproduced exactly with the **generic** series (E02); Boolean series gives 218.4/312.5/529.1 |
| C4 | Semi-regularity assumed; complexity from eq.(3) with ω=2 | p.12, p.35 | semi-regularity **false**: 8n syzygies at D=3 [VERIFIED n=11..16, every admissible e] and 2n falls at D=4 per S-box for every recommended exponent (E03) |
| C5 | "S_quad contains no quadratic equations in two distinct S-box outputs" | p.39–40 | true at degree 2 (n=11..13) **and 3** (n=14 **and n=16**, E12). The single cubic seen at n=16 is the *univariate* image criterion Tr(t₁^{−r₁}), degree wt(−r₁)=3, not a cross relation |
| C6 | Magic Pot not available: pt enters through Lin_in of degree 2^{n−1} over F_{2^n} | p.40 | correct for the original algorithm; skew version is open (E07) |
| C7 | Fast exhaustive search / polynomial method / linearisation-by-guessing not applicable | p.40–41 | not contested |
| C8 | Resultant attack (Sun et al.) mitigated by the input affine layer | p.41 | to be re-examined with the two-valued inverse (E05) |
| C9 | Multi-user: random iv's give no advantage except collisions | p.35 | to be examined (E09) |
| C10 | Brute force: 2^n forward evaluations; Grover cost 2^{165.6}/2^{232.5}/2^{298.0} | p.41–42 | reference only |

## 3. Structures established in this project (proof status per Focus Report + this repo)
A closed form S_e(x) = (x+1)^{2^e}/x [PROVED]; B image = trace-zero criterion, |Im|=2^{n−1} [PROVED]; C two-valued inverse with linear
section [PROVED]; D kernel variable a=y^r, alternating B_e, B_e(u,v)=uv(M_e(u)+M_e(v)) [PROVED]; E wt(r)=e⁻¹ mod n [PROVED];
F ⟨P1,P2,P3⟩=⟨P1⟩ in the Boolean ring, S_quad = XL₁(P1) [PROVED]; G syz₃ = 8n [VERIFIED n≤11, 3n proved]; H first-fall degree 4
with 2n falls [VERIFIED n=11..16, every admissible e, generic class]; I degenerate exponent classes, now a calibrated rule on the pair
{e, e⁻¹ mod n}: 7n / 4n / 17n·½ / 2n+1 / 2n [CONJ, fitted to 36 runs at n=11..16]; J Table 4 formula [VERIFIED];
K no cross-S-box cubic [VERIFIED n=14, 16]; direct sum of the D=4 falls of the first layer [VERIFIED n=11, ideal-theoretically n=14].
