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

---

## Log — 2026-08-31 (Tier A, batch 1) — **KAT pinning DONE**

**Status.** The encoding conventions are **pinned** against the official AIMer v3.0 KAT for all three
parameter sets [VERIFIED]. Full-size iv-dependent numbers may now be quoted.

**Source.** `https://aimer-signature.org` links the reference implementation at
`github.com/samsungsds-opensource/AIMer`; the repository holds `AIMer.zip` (commit `734cc76`,
*"update AIMer V3.0 reference implementation"*) containing `Reference_Implementation/` and `KAT/`.
`Reference_Implementation/aim3.c` / `aim3.h` carry the AIM3 exponents (3,7,11) / (11,23,47) /
(3,7,11,19), i.e. it is the v3.0 primitive. (Only that host was reachable from this session; the
GitHub API is out of scope here, the repository was read through the session's anonymous git proxy.)

**What the reference code fixes** (read, not guessed):
* `sign.c::crypto_sign_keypair`: `pk = iv || ct`, `sk = pt || iv || ct`.
* `field_common.c::gf_from_bytes`: a field element is little-endian bytes, so the §4.1.1 bit i ↔ X^i
  convention already used in `src/aim3/gf2n.py` is correct, with `int.from_bytes(buf, "little")`.
* `hash.h`: SHAKE128 for λ=128, SHAKE256 for 192/256 — as assumed.
* `aim3.c::squeeze_invertible_matrix`: each squeezed field element is stored as one **column**, the
  bits above the diagonal go to L and the bits below to U (in that stored/transposed form), and
  `mat[i] = U[i]·L`, i.e. **A = L·U** with L unit-lower and U unit-upper in the true matrix, applied
  as **A·pt**.

That last point is a convention axis the previous `Conventions` did not have: it distinguishes
A = L·U from A = U·L, and row-major from column-major filling of the tape matrix. `Conventions` now
has four axes — `tape_bit_order`, `matrix_apply`, `tape_matrix_order`, `lu_order` — and
`ALL_CONVENTIONS` enumerates the 16 combinations.

**Result.** `kat_pin.py` searched all 16 for the first three KAT records of every parameter set:

| set | n | matching conventions |
|---|---:|---|
| AIM3-I | 128 | `lsb / A / col / LU` and `lsb / AT / row / LU` |
| AIM3-III | 192 | same |
| AIM3-V | 256 | same |

The two survivors are the *same* matrix (transposing both the tape reading and the application
cancels), so the convention is unique. It is now the **default** in `src/aim3/aim3_ref.py`.
The `f` and `s` variants of each level share the KAT seed and hence the same AIM3 records; all six
directories were checked.

**Reproducibility.** The first three records of each set are committed to `kat_vectors.json`, so
`python3 experiments/E00_reference_and_facts/kat_pin.py --vectors` reproduces the table with no
network, and `tests/test_core.py::test_kat_conventions_pinned` fails if the default conventions ever
stop reproducing the official vectors. Output: `results/E00_kat_check.md`.

**Also done.** `generate_linear` was rewritten to assemble A from byte-aligned tape rows with
bitmask row operations (O(n²) big-int XORs instead of an O(n³) Python loop), so a full-size
evaluation is ~0.2 s at n=256 instead of minutes; this is what makes the 16-way search practical.

**Next step.** `docs/FACT_SHEET.md` §0.2 is updated. The remaining E00 item is the cost of one AIM3
evaluation in bit operations (Table 5 addition chain) for COMPLEXITY_ACCOUNTING §1 — still using the
placeholder 2^14.
