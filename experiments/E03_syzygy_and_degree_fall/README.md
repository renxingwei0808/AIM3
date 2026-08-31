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
