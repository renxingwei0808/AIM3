#!/usr/bin/env python3
"""Cross-S-box relations: dimension of Boolean relations of degree <= deg vanishing on the joint graph
{(t0,t1)} and {(pt,t0,t1)} of a toy AIM3 first layer (unit-LU matrices). The spec only excludes *quadratic*
cross relations; this measures cubic (and higher) ones. Meaningful when 2^n >> C(#vars, <=deg).
Usage: python3 cross_relations.py n deg [seed]     e.g.  n=14 deg=3 (2n=28 vars: 3683 monomials vs 16383 points)
"""
import sys, os, random, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS
from math import gcd, comb
from aim3_lab import F, sbox, rand_unit_lu, apply_mat, rank_gf2
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "E03_syzygy_and_degree_fall"))
from sbox_degree_fall import relation_dim
if __name__ == "__main__":
    n, deg = int(sys.argv[1]), int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed); fld = F(n)
    es = [e for e in range(3, n - 1) if gcd(e, n) == 1][:2]
    A0, A1 = rand_unit_lu(n, rng), rand_unit_lu(n, rng); b0, b1 = rng.randrange(fld.q), rng.randrange(fld.q)
    pts_t, pts_all = [], []
    for pt in range(fld.q):
        x0 = apply_mat(A0, pt, n) ^ b0; x1 = apply_mat(A1, pt, n) ^ b1
        if x0 == 0 or x1 == 0: continue
        t0, t1 = sbox(fld, x0, es[0]), sbox(fld, x1, es[1])
        pts_t.append(t0 | (t1 << n)); pts_all.append(pt | (t0 << n) | (t1 << (2 * n)))
    d_t, tot_t, _ = relation_dim(pts_t, 2 * n, deg)
    print(f"n={n} e={es} deg<={deg}: (t0,t1): dim={d_t} of {tot_t} monomials on {len(pts_t)} points  "
          f"[expected 0 if no cross relations; note the trace-image constraints Tr(t_i^-r)=0 have degree wt(-r)]")
    if 3 * n <= 48 and comb(3 * n, deg) < 2 ** n:
        d_a, tot_a, _ = relation_dim(pts_all, 3 * n, deg)
        print(f"   (pt,t0,t1): dim={d_a} of {tot_a}   [degree-3 span of the two S-box ideals predicts 2*(6n^2-5n)+? ; excess = cross falls]")
