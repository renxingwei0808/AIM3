# AIM3 — new algebraic analysis techniques for the AIMer v3.0 one-way function

Public research repository for an independent algebraic cryptanalysis of **AIM3** (the symmetric primitive of AIMer v3.0).
The goal is *not* to patch the designers' analysis but to develop **new algebraic techniques** from the new structures of AIM3
(kernel coordinates, twisted-bilinear form, two-valued inversion, skew-polynomial elimination) and to obtain independently
accounted complexity results against the official security claims. Technical roadmap: **`docs/ROADMAP.md`**.

## Layout
```
docs/         ROADMAP.md (technical route), FACT_SHEET.md (spec facts with page refs), COMPLEXITY_ACCOUNTING.md, EXPERIMENT_TEMPLATE.md
src/aim3/     pure-Python reference AIM3 (Alg. 4/5, any n) + S-box structure toolkit (kernel, section, two-valued inverse, audit)
legacy/       the six original verification scripts (unchanged; results in the Focus Report were produced with them)
experiments/  E00..E12, one directory per experiment (README = design + status; scripts; magma/ sage/ skeletons for the server)
results/      CSV outputs (schema in results/README.md)
server/       plan, sizing tables and job templates for the 2×Xeon 6230R / 256 GB server (Magma / Sage / msolve / M4RI)
paper/        outline, claims ledger, bibliography
tests/        fast regression tests (python3 -m pytest tests/)
```

## Quick start
```
python3 -m pip install -r requirements.txt          # (pytest only; core code has no dependencies)
python3 src/aim3/aim3_ref.py                        # full-size AIM3-I/III/V keygen + convention-independent self-test
python3 experiments/E02_official_baseline/table4_reproduce.py     # reproduces official Table 4 exactly
python3 experiments/E03_syzygy_and_degree_fall/sbox_degree_fall.py 12 5   # 8n syzygies + 2n degree falls
python3 experiments/E06_kernel_bilinear_model/kernel_bilinear_toy.py 7 3  # kernel bilinear model vs S_quad
python3 experiments/E08_parameter_audit/exponent_audit.py        # exponent audit for n=128/192/256
bash scripts/smoke_test.sh
```

## Status board (2026-08-31)
| Track | Name | State |
|---|---|---|
| T1 | Structured Hilbert function / true solving degree | syz₃=8n, first-fall 4 with 2n falls verified n≤16; theorem + full-system falls pending (E03/E04) |
| T2 | Twisted-bilinear & kernel models, block orders, hybrid | 6n model built; falls observed; server F4 runs pending (E06/E10/E11) |
| T3 | Two-valued inverse reduced systems | statistics done; resultant experiments pending (E05) |
| T4 | Skew-polynomial (Ore) elimination | skeleton (E07) |
| T5 | Parameter audit / weak exponent classes | criteria + audit script done (E08) |
| T6 | Multi-target / reused key | planned (E09) |

## Provenance and honesty rules
Every claim carries a status tag ([PROVED]/[VERIFIED n≤…]/[CONJ]/[OPEN]); every number in the paper traces to a results CSV
(commit, tier, host). Baselines B1/B2/B3 are defined in `experiments/E02_official_baseline/README.md`. No full-size number is
quoted before the KAT pinning of E00.

## Working with AI agents
`CLAUDE.md` describes the conventions for Claude Cowork (tier A work, writing) and Claude Code on the local machine driving
the remote server over SSH (tiers B/C, Magma). Human review is required before any experiment is marked *done*.
