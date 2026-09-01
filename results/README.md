# results/
CSV files written by the experiment scripts through `experiments/_common.py::record()`, which appends the
provenance columns **`tier, host, commit, seconds, censored`** to every row (CLAUDE.md ground rule 2).
`censored = 1` means the run hit its time budget: the numbers in that row are **lower bounds, never values**
(ground rule 4), and the columns that were not reached are left empty.

| file | written by | content |
|---|---|---|
| `E00_kat_check.md` | `E00_reference_and_facts/kat_pin.py --md` | which encoding convention reproduces the official AIMer v3.0 KAT, per parameter set |
| `E02_baselines.csv` | `E02_official_baseline/table4_reproduce.py --csv` | official Table 4 reproduction and the baselines B1/B2; the (log2T, k, d_reg) triples are now flattened into scalar columns |
| `E03_sbox_degree_fall.csv` | `E03_syzygy_and_degree_fall/sbox_degree_fall.py` (driver `sweep.sh`) | per-(n,e): dim I(V)_{≤3}, Macaulay ranks at D=3,4, syz₃, falls_D4, cubics outside the ideal, degeneracy class. Append-only |
| `E03_summary.md` | `E03_syzygy_and_degree_fall/summarise.py --md` | the full 36-run table plus the calibrated fall-count rule and its checks |
| `E04_full_system_falls.csv` | `E04_full_system_hilbert/full_system_falls.py` (driver `run_tierA.sh`) | full-system Macaulay statistics at D=3,4 for S_quad / P1-only / random control, with `low_cols`, `saturation`, the per-block falls and `cross_excess` |
| `E05_reduced_system_stats.csv` | `E05_two_valued_inverse/reduced_system_toy.py` | term counts, ordinary degrees in t₀/t₁, Boolean-degree histograms, Bézout bounds and exact root counts of the two-valued-inverse reduced system |
| `E06_kernel_model_falls.csv` | `E06_kernel_bilinear_model/kernel_bilinear_toy.py` | kernel model vs S_quad vs P1-only vs random bilinear vs random quadratic, per degree D |
| `E08_audit_{128,192,256}.csv` | `E08_parameter_audit/exponent_audit.py` | full exponent audit: wt(r), wt(−r), balance, degeneracy class, predicted falls_D4, recommended flag |
| `E12_cross_relations.csv` | `E12_cross_sbox_relations/cross_relations.py` | dimension of the low-degree relation space on (t₀,t₁) and on (pt,t₀,t₁) |
| `E12_first_layer_ideal.csv` | `E12_cross_sbox_relations/first_layer_ideal.py` | first-layer Macaulay ranks, the point-set ideal dimension, and the per-S-box direct-sum control |
| `E12_cross_relations_explained.csv` | `E12_cross_sbox_relations/explain_n16.py` | splits the cubic relation space into univariate image criteria vs genuine cross relations |
| `REPORT_tierA_batch1.md` | written by hand | the batch report: anchors, results, Gate 1, anomalies, what is left for the server |
| `REPORT_tierA_batch2.md` | written by hand | batch-2 report: the refuted T1.2 ansatz, the two indispensable families, the kernel-model generators, the bit cost |
| `REPORT_server_batch1.md` | written by hand | server batch-1 report: Gate 1 at n = 13..19, real `d_solv`, E05 resultants, E12 at degree 4 |

**Tier-A batch 2 — theory tables** (all written through `record()`, same provenance columns):

| file | written by | content |
|---|---|---|
| `E01_theory_matrix_distribution.csv` | `E01_sbox_structures/matrix_distribution_control.py` | first-layer D=3 rank/syzygies and D=4 falls under unit-LU vs uniform GL_n vs same-density random invertible matrices |
| `E03_theory_fall_generators.csv` | `E03_syzygy_and_degree_fall/fall_generators.py --sbox/--all` | per (n,e): what the proposed T1.2 ansatz explains (0), the per-family sweep, and the leave-one-out result |
| `E03_theory_fall_basis.csv` | `fall_generators.py --basis` | one row per element of a basis of the 2n new cubics: kernel-vector size, multiplier weights, which of P1/P2/P3 it uses |
| `E03_theory_orbit_lemma.csv` | `fall_generators.py --orbits` | the Frobenius-orbit lemma and the closed form for the ansatz module's rank/kernel, n = 5..20 |
| `E06_theory_fall_generators.csv` | `fall_generators.py --kernel` | the kernel model's 6n new quadratics at D=3: the P1_i candidate (fails) vs a_i·K_i and t_i·K_i (explains everything) |
| `E09_theory_reused_key.csv` | `E09_multitarget_weak_iv/reused_key_toy.py` | q instances sharing one pt: D=3 rank, syzygies, saturation, and the excess over the direct sum |
| `E11_theory_eval_cost.csv` | `E11_hybrid_estimates/aim3_eval_cost.py --csv` | bit cost of one forward AIM3 evaluation from the spec Table 5 addition chains, schoolbook and Karatsuba |

**Server batch 1 — summary tables copied back** (large matrices and logs stay on the server under
`results/server/<date>_<exp>/`):

| file | content |
|---|---|
| `E03_server_n17_20.csv` | the E03 single-S-box sweep extended to n = 17..20 |
| `E03_server_D5.csv` | the single-S-box Macaulay statistics at D = 5 |
| `E04_server_gate1.csv` | Gate 1: full-system D=4 falls and `cross_excess` at ℓ = 1,2,3 and n = 11..19, with saturation |
| `E04_server_gate1_pairs.csv` | the same over *every subset* of the S-box blocks — the strong form of the direct-sum test |
| `E04_f4_degrees.csv` | real solving degrees from msolve and Magma F4 step logs (see the exponent caveat, COMPLEXITY_ACCOUNTING §5) |
| `E05_server_resultants.csv` | E05 resultant timings, degrees, root and factor counts (n = 7 completes, n = 8 censored) |
| `E06_server_kernel.csv` | the kernel model and its random bilinear control at n = 11, 13, 15 |
| `E10_server_block_orders.csv` | block-order runs — read together with the tool caveat: Magma's GF(2) F4 computes in grevlex whatever order is declared |
| `E12_server_first_layer.csv` | the first-layer system at n = 13, 14 |
| `E12_server_cross_relations_deg4.csv` | the degree-4 relation space on (t0,t1) at n = 16 |
| `E12_server_cross_relations_explained_deg4.csv` | its naive univariate/cross split (the reading that needed correcting) |
| `E12_server_criteria_subtract.csv` | the corrected test: subtract the whole ideal generated by the image criteria — 0 genuine cross relations at degree ≤ 4 |

Large matrices / logs from the server go to `results/server/<date>_<exp>/` (git-ignored if > 5 MB; keep the summary CSV).
