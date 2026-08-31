#!/usr/bin/env python3
"""Solving-degree comparison for a two-S-box toy AIM3 (l=1)."""
import itertools, sys
from aim3_macaulay import (build_aim3, random_system, monomial_index,
                           poly_to_row, rank_bits)


def stats(eqs, nv, D):
    idx, ncol, by_deg = monomial_index(nv, D)
    masks = []
    acc = 0
    for d in range(D + 1):
        for c in by_deg[d]: acc |= 1 << c
        masks.append(acc)          # masks[d] = columns of degree <= d
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
    out = []
    for d in range(D):                      # dim(span ∩ deg<=d)
        top = ~masks[d]
        rk_top = rank_bits([r & top for r in rows])
        out.append(rk - rk_top)
    return len(rows), rk, ncol, out


def run(name, eqs, nv, Dmax):
    print(f"  {name}  ({len(eqs)} eqs, {nv} vars)")
    for D in range(2, Dmax + 1):
        nrow, rk, ncol, inter = stats(eqs, nv, D)
        solved = "  *** SOLVED (all linear forms recovered)" if inter[1] >= nv else ""
        print(f"    D={D}: rows={nrow:6d} cols={ncol:6d} rank={rk:6d} "
              f"| dim(∩deg<=1)={inter[1]:4d} dim(∩deg<=2)={inter[2] if D>2 else '-':>5} "
              f"{solved}")
    print()


if __name__ == "__main__":
    n = int(sys.argv[1]); Dmax = int(sys.argv[2]); nexp = int(sys.argv[3])
    from math import gcd
    cand = [e for e in range(2, n - 1) if gcd(e, n) == 1]
    exps = cand[:nexp]
    print(f"toy AIM3 n={n}, l={len(exps)-1}, exponents={exps}")
    nv, eqsP1, eqsAll, sol = build_aim3(n, exps)
    run("AIM3 S_quad (P1+P2+P3)", eqsAll, nv, Dmax)
    run("AIM3 P1-only          ", eqsP1, nv, Dmax)
    r1, _ = random_system(nv, len(eqsAll), seed=11)
    run("random, S_quad shape  ", r1, nv, Dmax)
    r2, _ = random_system(nv, len(eqsP1), seed=13)
    run("random, P1 shape      ", r2, nv, Dmax)
