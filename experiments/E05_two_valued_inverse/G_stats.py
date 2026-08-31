#!/usr/bin/env python3
"""Exponent statistics of the two-valued inverse G_{e,eps}(y) = y^r H_e(y^{-r}) + eps y^r + 1.
As an element of F_{2^n}[y], the H-part has the n exponents r(1-2^j) mod (2^n-1); we report their Boolean
degrees (Hamming weights), the number of distinct cyclotomic cosets, and the ordinary degree.
Usage: python3 G_stats.py            (recommended parameters)   |  python3 G_stats.py n e
"""
import sys
def coset_rep(x, n):
    q1 = (1 << n) - 1; v = x % q1; s = {v}
    for _ in range(n):
        v = (v * 2) % q1; s.add(v)
    return min(s)
def stats(n, e):
    q1 = (1 << n) - 1
    r = pow((1 << e) - 1, -1, q1)
    exps = [(r * (1 - (1 << j))) % q1 for j in range(n)]
    wts = [bin(x).count("1") for x in exps]
    return dict(n=n, e=e, wt_r=bin(r).count("1"), terms=n + 2, min_wt=min(wts), max_wt=max(wts),
                cosets=len({coset_rep(x, n) for x in exps}), maxexp_over_2n=round(max(exps) / 2 ** n, 4))
if __name__ == "__main__":
    if len(sys.argv) == 3:
        print(stats(int(sys.argv[1]), int(sys.argv[2])))
    else:
        for name, n, es in [("AIM3-I", 128, [3, 7, 11]), ("AIM3-III", 192, [11, 23, 47]), ("AIM3-V", 256, [3, 7, 11, 19])]:
            for e in es:
                print(name, stats(n, e))
