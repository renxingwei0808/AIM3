#!/usr/bin/env python3
"""
Single enhanced-Mersenne-S-box: Boolean quadratic system and its degree-3
syzygy module, compared with a random quadratic system of the same shape.

Variables: x bits 0..n-1, y bits n..2n-1.
Equations: P1 : xy + x^{2^e} + 1
           P2 : x^2 y + x^{2^e+1} + x
           P3 : x y^2 + x^{2^e} y + y
"""
from __future__ import annotations
import itertools, random
from math import gcd
from aim3_lab import F
from aim3_macaulay import (AffExpr, prod_coords, lin_coords, add_polys,
                           monomial_index, poly_to_row, rank_bits, random_system)


def sbox_system(n, e):
    fld = F(n)
    X = AffExpr(fld, {k: 1 << k for k in range(n)}, 0)
    Y = AffExpr(fld, {n + k: 1 << k for k in range(n)}, 0)
    X2 = X.frob(1); Xe = X.frob(e); Y2 = Y.frob(1)
    one = AffExpr(fld, {}, 1)
    P1 = add_polys(prod_coords(X, Y), lin_coords(Xe), lin_coords(one))
    P2 = add_polys(prod_coords(X2, Y), prod_coords(Xe, X), lin_coords(X))
    P3 = add_polys(prod_coords(X, Y2), prod_coords(Xe, Y), lin_coords(Y))
    return fld, P1, P2, P3


def check_vanishing(fld, eqs, e):
    n = fld.n
    for _ in range(30):
        x = random.randrange(1, fld.q)
        y = fld.pw(x, (1 << e) - 1) ^ fld.inv(x)
        sol = [(x >> k) & 1 for k in range(n)] + [(y >> k) & 1 for k in range(n)]
        for p in eqs:
            s = 0
            for m in p:
                v = 1
                for j in m: v &= sol[j]
                s ^= v
            if s: return False
    return True


def macaulay_stats(eqs, nv, D):
    idx, ncol, by_deg = monomial_index(nv, D)
    topmask = 0
    for c in by_deg[D]: topmask |= 1 << c
    mults = [frozenset(m) for d in range(D - 1)
             for m in itertools.combinations(range(nv), d)]
    rows = []
    for f in eqs:
        for mu in mults:
            g = set()
            for m in f: g ^= {m | mu}
            if not g: continue
            if max(len(m) for m in g) > D: continue
            rows.append(poly_to_row(g, idx))
    rk = rank_bits(rows)
    rk_top = rank_bits([r & topmask for r in rows])
    return len(rows), rk, rk - rk_top, ncol


if __name__ == "__main__":
    random.seed(1)
    print("single S-box, variables (x,y) = 2n bits, equations P1|P2|P3 = 3n")
    print(f"{'n':>3} {'e':>3} | {'rows3':>6} {'rank3':>6} {'syz3':>6} "
          f"{'syz3/n':>7} | {'rnd rows':>8} {'rnd rank':>8} {'rnd syz':>7} | "
          f"{'fall(P1)':>8}")
    for n in (5, 6, 7, 8, 9, 10, 11):
        for e in range(2, n - 1):
            if gcd(e, n) != 1: continue
            fld, P1, P2, P3 = sbox_system(n, e)
            allq = P1 + P2 + P3
            assert check_vanishing(fld, allq, e), (n, e)
            nv = 2 * n
            rows, rk, low, ncol = macaulay_stats(allq, nv, 3)
            syz = rows - rk
            reqs, _ = random_system(nv, 3 * n, seed=n * 100 + e)
            rrows, rrk, rlow, _ = macaulay_stats(reqs, nv, 3)
            # P1-only : how many new degree<=2 elements appear at D=3 ?
            _, rk1, low1, _ = macaulay_stats(P1, nv, 3)
            print(f"{n:3d} {e:3d} | {rows:6d} {rk:6d} {syz:6d} {syz/n:7.2f} | "
                  f"{rrows:8d} {rrk:8d} {rrows-rrk:7d} | {low1:8d}  (=n+{low1-n})")
