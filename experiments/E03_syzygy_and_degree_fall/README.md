# E03 — Syzygies (8n at D=3) and degree falls (2n at D=4) of a single S-box

**Hypotheses.** For a single enhanced-Mersenne S-box with the official quadratics {P1,P2,P3} (3n equations, 2n variables):
* (H1) [Focus Report, verified n≤11] the degree-3 Macaulay matrix has exactly 8n syzygies (3n trivial ones xP1=P2, yP1=P3,
  yP2=xP3, plus 5n from Frobenius-orbit cancellations).
* (H2) [**new, verified n=11..16, generic e**] the degree-4 Macaulay matrix produces exactly **2n** new cubics (first-fall
  degree 4, 2n falls), and these exhaust the cubic part of the vanishing ideal I(V_e) — i.e. no cubic generators outside
  ⟨P1,P2,P3⟩. Verified: (n,e) = (11,3),(11,4),(13,3),(14,3),(14,5),(16,3) → 2n; some (n,e) give 2n+1.
* (H3) [**new**] degenerate exponent classes inflate the count: e=2 → 9n; 2e≡−1 (mod n) → 4n (n=13,e=6); e=1+n/2 → ≈8.5n
  (n=12,e=7). None of the recommended exponents is in these classes (checked in E08).

**Data so far** (`results/E03_sbox_degree_fall.csv` is appended by the script):

| n | e | dim I(V)_{≤3} | Mac3 rank | syz3 | falls D4→3 | outside ideal |
|---|---|---|---|---|---|---|
| 11 | 3 | 693 | 671 | 88 | 22 | 0 |
| 12 | 5 | 829 | 804 | 96 | 25 | 0 |
| 12 | 7 | 906 | 804 | 96 | 102 | 0 |
| 13 | 3 | 975 | 949 | 104 | 26 | 0 |
| 13 | 6 | 1001 | 949 | 104 | 52 | 0 |
| 14 | 3/5 | 1134 | 1106 | 112 | (2n=28 by excess) | – |
| 16 | 3 | 1488 | 1456 | 128 | (32) | – |

**Method.** `sbox_degree_fall.py n e` — (i) rank of the evaluation matrix of all monomials of degree ≤3 on the 2^n−1 graph
points (meaningful when 2^n−1 > C(2n,≤3), i.e. n ≥ 14); (ii) Boolean Macaulay ranks at D=3,4 with multipliers of degree ≤ D−2.

**Next steps.**
1. Tier A: fill the table for all admissible e at n=11..14 (D=4 takes seconds to minutes).
2. Tier C (server, Magma/M4RI): D=4 at n=16..24 and D=5 at n=11..14 (see `../../server/SIZING.md`); `magma/sbox_macaulay.m` is a skeleton.
3. Theory (T1): write the 2n cubics explicitly. Candidate mechanism: products (x^{2^i}y^{2^j})·P1 whose degree-4 parts cancel
   pairwise along Frobenius orbits (same mechanism as the 5n syzygies). Then prove: first-fall degree 4, fall count 2n, and a
   closed form for the exceptional classes. Target statement: *the Boolean ideal of one S-box has bigraded Hilbert function
   H(d) = ... for d ≤ 4*, to be fed into E04.

---

## Log — 2026-08-31 (Tier A, batch 1)

**Status.** syz₃ = 8n [VERIFIED n=11..16, every admissible e]. First-fall degree 4 [VERIFIED n=11..16].
Fall count = 2n for generic e, with four exceptional classes [CONJ, fitted to 36 runs].

**What was run.** `bash sweep.sh` (n = 11,12,13,14,15,16; every e with gcd(e,n)=1 and 2 ≤ e ≤ n−2;
36 runs, 0 censored, ≤ 81 s each) → `results/E03_sbox_degree_fall.csv`;
`python3 summarise.py --md` → **`results/E03_summary.md`** (the full table and the analysis).
The batch asked only for n = 11..14; n = 15 and 16 were added because the fast Macaulay routine
(`src/aim3/boolmac.py`) made them cost minutes, and they are what pins the 2n vs 2n+1 rule.

**Results.**
1. **syz₃ = 8n in all 36 runs**, with no exception.
2. **`outside_ideal` = 0 in 33 of 36 runs**: the D=4 Macaulay matrix exhausts the cubic part of
   I(V_e), so there are no cubic generators outside ⟨P1,P2,P3⟩. The three exceptions are exactly
   **e = 2** (n = 11, 13, 15), where 2n cubics of the ideal are *not* reached at D=4. Note the
   asymmetry: e = 2⁻¹ has the same fall count but `outside_ideal` = 0. [OPEN] why.
3. **falls_D4 depends only on the unordered pair {e, e⁻¹ mod n}** (checked on every inverse pair in
   the table) and equals 2n except for:

   | class of e | falls_D4 |
   |---|---|
   | {2, 2⁻¹} | 7n |
   | {−2, (−2)⁻¹} | 4n |
   | {1+n/2} (n even) | 17n/2 |
   | {−3, (−3)⁻¹} | 2n+1 |
   | {n/2−1} (n even) | 2n+1 |
   | generic | 2n |

4. **The 2n vs 2n+1 discriminator is found** (the batch allowed us to answer [OPEN]): the extra fall
   appears iff 3e ≡ −1 (mod n) or 3e⁻¹ ≡ −1 (mod n) or e = n/2 − 1. The **sign matters** —
   3e ≡ +1 gives the generic 2n ((11,4), (13,9), (14,5), (16,11) all check out.)

**Correction to H3 as previously stated.** e = 2 gives **7n**, not 9n; "2e ≡ ±1 → 4n" is wrong as
written — 2e ≡ −1 gives 4n but 2e ≡ +1 gives 7n (because 2e ≡ 1 means e = 2⁻¹, i.e. the class of
e = 2); e = 1+n/2 gives exactly 17n/2 = 8.5n.

**Why `aim3.structures.degenerate_class` was changed.** The previous rule tagged every e with
k·e ≡ ±1, ±2 (mod n) for k ≤ 4. Against the sweep it is wrong in both directions: 11 of the 36
exponents carried a tag and are perfectly generic (e.g. (11,3), (13,9), (14,5), (16,11)), and 3
exponents with 4n falls carried no tag at all ((11,9), (13,11), (15,13)). It has been replaced by
the class rule above, which reproduces all 36 measured fall counts exactly; the old behaviour is
still available as `degenerate_class(n, e, legacy_kmax=4)` so that older CSV rows can be
reproduced. `predicted_falls_D4(n, e)` returns the predicted count. All recommended AIM3 exponents
remain in the generic class (asserted by `exponent_audit.py` and `tests/test_core.py`).

**Next steps.**
* Theory (T1, unchanged): write the 2n cubics explicitly; the class rule above is the shape the
  theorem must reproduce, including the sign asymmetry of 3e ≡ ±1 and the e ↔ e⁻¹ symmetry.
* **Server (Tier B/C):** the class rule is a fit to n ≤ 16 and must be tested outside it.
  `bash experiments/E03_syzygy_and_degree_fall/sweep.sh 17 18 19 20` under Magma/M4RI, and
  D = 5 at n = 11..14 (`magma/sbox_macaulay.m`, still an untested skeleton — validate against
  `results/E03_sbox_degree_fall.csv` at n ≤ 13 first).
* [OPEN] the e = 2 asymmetry in `outside_ideal`.
