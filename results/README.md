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

Large matrices / logs from the server go to `results/server/<date>_<exp>/` (git-ignored if > 5 MB; keep the summary CSV).
