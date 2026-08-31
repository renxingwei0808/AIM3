# E00 — Reference implementation, KAT pinning, and the AIM3 fact sheet

**Goal.** A trusted base: (1) `src/aim3/aim3_ref.py` reproduces the official AIMer v3.0 test vectors; (2) every factual
claim about AIM3 used in the paper is listed in `docs/FACT_SHEET.md` with a spec page reference.

**What is already done.** `aim3_ref.py` implements Algorithm 4 (GenerateLinear: unit-LU matrices from a SHAKE tape of
2ℓn²+(ℓ+1)n bits) and Algorithm 5 (NoZeroInpSB) for any n, and passes convention-independent self-tests at full size
(P1, kernel identities, second preimage, trace criterion).

**What is NOT done (required before quoting any full-size number).** Two encoding conventions are not fixed by the spec text
and are parameters of `Conventions`: `tape_bit_order` (lsb/msb within a tape byte) and `matrix_apply` (A·v vs Aᵀ·v; Alg. 4 says
matrices are produced "in their transposed form"). Pin them as follows (tier A or B, ~1 hour):
1. Download the official reference implementation / KAT files from https://aimer-signature.org (or the KpqC/NIST submission package).
2. For each parameter set take one KAT `(pt, iv, ct)`; loop over the 4 convention combinations; keep the one that reproduces `ct`.
3. Record the winning `Conventions(...)` in `docs/FACT_SHEET.md` §0 and set it as the default in `aim3_ref.py`.
4. Also confirm that the KAT's `iv` byte string is fed to SHAKE exactly as stored (no length prefix / domain byte).

**Outputs.** `results/E00_kat_check.md` (which convention matched, for each of the three parameter sets).

**Note.** All *toy-n* experiments in this repository are convention-independent (they use random unit-LU matrices from a PRNG),
so the KAT pinning only matters for experiments on the real parameter sets (E08 statistics, E09 weak-iv search).
