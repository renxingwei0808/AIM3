# legacy/ — original verification scripts (unchanged)
These six scripts produced the numbers in `AIM3_Focus_Report.md` (2026-08-31) and the earlier report (2026-08-28).
They are kept verbatim for provenance; new code lives in `src/aim3/` and `experiments/`. Run from this directory or via
`experiments/_common.py` (which puts this directory on `sys.path`). Field implementation and toy moduli (n ≤ 16) are shared with
`src/aim3/gf2n.py`, so results are directly comparable.

| file | purpose |
|---|---|
| `aim3_sbox_probe.py` | exhaustive image / fiber / trace criterion / two-valued inverse checks (n=5,7,8) |
| `aim3_lab.py` | structures A/B/C/E; wt(r) closed form at recommended parameters; cross-S-box quadratic relations (n=11..13) |
| `aim3_kernel_model.py` | kernel variable identities and dim QE(x,a) |
| `aim3_macaulay.py` | symbolic S_quad / P1-only toy AIM3 systems and Boolean Macaulay ranks |
| `aim3_sbox_syz.py` | the 8n syzygy data (single S-box, D=3) |
| `aim3_solvdeg.py` | solving-degree comparison AIM3 vs random |
