# Dense Macaulay matrix sizing (rows = #eqs · #multipliers of degree ≤ D−2, cols = #monomials of degree ≤ D; bits)

Full toy AIM3, S_quad (9n equations, N = 3n variables):

| n | N | D | cols C(N,≤D) | rows 9n·C(N,≤D−2) | bits | GB | tool |
|---|---|---|---|---|---|---|---|
| 11 | 33 | 4 | 46,376 | 55,638 | 2.6e9 | 0.3 | Python ok (slow) / M4RI |
| 11 | 33 | 5 | 283,000 | 592,000 | 1.7e11 | 21 | M4RI |
| 13 | 39 | 5 | 660,000 | 1,160,000 | 7.7e11 | 96 | M4RI, one job |
| 16 | 48 | 4 | 214,000 | 169,000 | 3.6e10 | 4.5 | M4RI |
| 11 | 33 | 6 | 1.4e6 | 4.6e6 | 6.4e12 | 800 | **F4/msolve only** |

Single S-box (3n equations, 2n variables): D=4 at n=20 → cols C(40,≤4)=102,091, rows 60·821=49,260 → 0.6 GB (M4RI, minutes);
D=5 at n=16 → cols 242,825, rows 48·5,489=263,000 → 8 GB.

Rule: above ~50 GB dense, switch to F4 (Magma / msolve) and record step degrees instead of full ranks.
