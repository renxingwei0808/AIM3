# Experiments — index and conventions

Each `Exx_*/` directory is one self-contained experiment with a `README.md` (goal, hypothesis, method, status, compute tier,
outputs, next steps) and scripts. Results go to `../results/` as CSV (schema in `../results/README.md`). Never overwrite a
result file produced on another tier; append with a `tier` column.

| ID | Name | Track | Tier | Status (2026-08-31) |
|---|---|---|---|---|
| E00 | Reference implementation & fact sheet | all | A | ref impl written; KAT pinning pending |
| E01 | S-box structures A–G (closed form, image, kernel, weights) | T5 | A | verified n≤11 (legacy scripts) |
| E02 | Official Table 4 reproduction & baselines B1/B2 | T1 | A | **done** (exact reproduction) |
| E03 | Syzygies (8n) and degree falls (2n at D=4) per S-box | T1 | A→C | verified n≤16; theorem pending |
| E04 | Full-system structured Hilbert function / solving degree | T1 | B/C | to do (server) |
| E05 | Two-valued inverse reduced system & resultants | T3 | A→C | stats done; resultants pending |
| E06 | Kernel-variable bilinear (6n) model | T2 | A→C | generator done; falls observed |
| E07 | Ore / skew-polynomial elimination | T4 | B/C | skeleton |
| E08 | Exponent & linear-layer audit | T5 | A | audit script done |
| E09 | Multi-target / weak-iv / reused-key | T6 | A/B | to do |
| E10 | Block-order Gröbner (FreeLunch-style) | T2 | C | skeleton |
| E11 | Generic hybrid estimates (FSS / y-XL) | T2 | A | done |
| E12 | Cross-S-box relations of degree ≥3 | T1 | A/B | cubic (t0,t1)=0 verified n=14 |

Compute tiers: **A** = pure Python, runs in Claude Cowork / laptop (n ≤ 16, minutes). **B** = server via Claude Code:
Sage / msolve / M4RI (Boolean n ≤ 32). **C** = server Magma (F4 with block orders, resultants, skew polynomials).

Run anything from the repo root, e.g. `python3 experiments/E03_syzygy_and_degree_fall/sbox_degree_fall.py 12 5`.
