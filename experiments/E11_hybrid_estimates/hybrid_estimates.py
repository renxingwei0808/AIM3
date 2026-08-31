#!/usr/bin/env python3
"""Generic complexity curves for the *bilinear* view of the P1 system (blocks X = pt bits, T = S-box outputs):
  (a) FSS bound for affine bilinear systems: d_reg <= min(n_X, n_T) + 2 after guessing k bits of pt;
  (b) y-XL linearisation degree  d_T >= (n_X+1) n_T / (m - n_X - 1)  (multiply only by T-monomials);
  (c) exhaustive search reference B3 = 2^n * cost(AIM3), with cost(AIM3) now MEASURED from the addition
      chains of spec Table 5 by `aim3_eval_cost.py` (2^19.9..2^23.5 bit ops depending on n and on whether
      field multiplication is costed schoolbook or Karatsuba) -- it used to be a 2^14 placeholder;
  (d) the E06 kernel-variable model's shape (6n equations, blocks U = (pt,t0,t1) and W = (a0,a1,a2), each
      3n bits): the same FSS bound after guessing k bits of the U block, so n_U = 3n - k and n_W = 3n.
      This is the *generic* bilinear bound for that shape; E06 shows the model is NOT generic (it has 6n
      degree falls at D=3 where a random bilinear system of the same shape has none), so this curve is an
      upper bound on what a structure-agnostic solver would do, not a prediction for the real model.
These are *generic* formulas; they bound what structure-agnostic bilinear solvers can do and show that no
hybrid of this family beats brute force. Output: log2 costs vs k for the three parameter sets.
"""
import sys, os
from math import comb, log2, ceil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aim3_eval_cost import eval_cost, SPEC_TABLE5
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


def kernel_curve(n):
    """E06 kernel model: 6n bilinear equations, U = (pt,t0,t1) and W = (a0,a1,a2), 3n bits each.
    Guess k bits of U (n_U = 3n - k, n_W = 3n) and apply the same affine-bilinear FSS bound.

    k is capped at n: pt alone determines t0, t1, a0, a1, a2, so guessing more than n bits of U buys
    nothing an attacker could not get by guessing pt and evaluating forward (that is B3)."""
    nW = 3 * n
    out = []
    for k in range(0, n + 1, 8):
        nU = 3 * n - k
        d = min(nU, nW) + 2
        c = k + 2 * log2(comb(nU + nW, d)) if d <= nU + nW else float("inf")
        out.append((k, nU, d, round(c, 1)))
    return out
if __name__ == "__main__":
    for name, n, ell in [("AIM3-I", 128, 2), ("AIM3-III", 192, 2), ("AIM3-V", 256, 3)]:
        sb = eval_cost(name, "schoolbook")
        ka = eval_cost(name, "karatsuba")
        print(f"== {name}: B3 = 2^{n} AIM3 evaluations = 2^{sb['log2_B3']} bit ops (schoolbook) / "
              f"2^{ka['log2_B3']} (Karatsuba)   [E11 aim3_eval_cost.py, spec Table 5]")
        print("   k   FSS d_reg  log2 cost  |  y-XL d_T  log2 cost")
        for k, dF, cF, dT, cX in curves(n, ell):
            print(f"  {k:4d}  {dF:5d}   {cF:9.1f}  |  {dT:5d}   {cX:9.1f}")
        kc = kernel_curve(n)
        best = min(kc, key=lambda r: r[3])
        print(f"   kernel model (6n/6n, blocks 3n|3n), generic FSS bound (k <= n): best is "
              f"k={best[0]} -> n_U={best[1]}, d_reg {best[2]}, log2 cost {best[3]}  "
              f"({'above' if best[3] > sb['log2_B3'] else 'BELOW'} B3 = 2^{sb['log2_B3']})")
