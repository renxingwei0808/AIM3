#!/usr/bin/env python3
"""Generic complexity curves for the *bilinear* view of the P1 system (blocks X = pt bits, T = S-box outputs):
  (a) FSS bound for affine bilinear systems: d_reg <= min(n_X, n_T) + 2 after guessing k bits of pt;
  (b) y-XL linearisation degree  d_T >= (n_X+1) n_T / (m - n_X - 1)  (multiply only by T-monomials);
  (c) exhaustive search reference 2^n * cost(AIM3) with cost(AIM3) ~ 2^14 bit-ops (placeholder, refine in E00).
These are *generic* formulas; they bound what structure-agnostic bilinear solvers can do and show that no
hybrid of this family beats brute force. Output: log2 costs vs k for the three parameter sets.
"""
from math import comb, log2, ceil
def curves(n, ell):
    nT = ell * n; m = (ell + 1) * n
    out = []
    for k in range(0, n + 1, 4):
        nX = n - k
        dF = nX + 2
        cF = k + 2 * log2(comb(nX + nT, dF)) if dF <= nX + nT else float("inf")
        dT = ceil((nX + 1) * nT / max(m - nX - 1, 1))
        cX = k + 2 * log2(m * comb(nT + dT - 1, dT)) if dT <= nT else float("inf")
        out.append((k, dF, round(cF, 1), dT, round(cX, 1)))
    return out
if __name__ == "__main__":
    for name, n, ell in [("AIM3-I", 128, 2), ("AIM3-III", 192, 2), ("AIM3-V", 256, 3)]:
        print(f"== {name}: brute force ~ 2^{n} AIM3 calls (~2^{n+14} bit ops)")
        print("   k   FSS d_reg  log2 cost  |  y-XL d_T  log2 cost")
        for k, dF, cF, dT, cX in curves(n, ell):
            print(f"  {k:4d}  {dF:5d}   {cF:9.1f}  |  {dT:5d}   {cX:9.1f}")
