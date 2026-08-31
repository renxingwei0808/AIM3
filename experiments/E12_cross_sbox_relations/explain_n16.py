#!/usr/bin/env python3
"""E12 — why `cross_relations.py 16 3` reports dim = 1 while n=14 reports 0.

`cross_relations.py n 3` measures the space of Boolean cubics vanishing on the joint graph
{(t0, t1)} of the two first-layer S-boxes.  At n=14 (e = 3, 5) the dimension is 0; at n=16 (e = 3, 5)
it is 1.  The candidate explanation is NOT a cross-S-box relation but the UNIVARIATE image criterion

    t_i in Im S_{e_i}  <=>  Tr(t_i^{-r_i}) = 0,

whose Boolean degree is wt(-r_i) = n - (e_i^{-1} mod n).  For (n,e) = (16,5): e^{-1} mod 16 = 13, so
wt(-r) = 3 and the criterion IS a cubic; for every exponent at n=14 the degree is 9 or 11, too high.

This script decides it: it computes the relation space, and checks whether it is spanned by a
polynomial in the t1 variables alone that agrees with Tr(t1^{-r_1}).

Usage: python3 explain_n16.py [n] [seed]
"""
import sys, os, random, time
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3_lab import F, sbox, rand_unit_lu, apply_mat, rank_gf2
import itertools


def anf_of_trace_criterion(fld, e, n):
    """the Boolean function t -> Tr(t^{-r}) on t != 0, as a truth table over all 2^n values of t
    (value at t = 0 is set to 1, since 0 is in the image but 0^{-r} is undefined -- the graph never
    contains t = 0 for the instances used here)."""
    r = fld.pw  # noqa
    q1 = fld.q - 1
    rr = pow((1 << e) - 1, -1, q1)
    tt = [0] * fld.q
    for t in range(1, fld.q):
        tt[t] = fld.trace(fld.pw(t, (-rr) % q1))
    return tt


def main(n, seed=1):
    t_start = time.time()
    rng = random.Random(seed)
    fld = F(n)
    es = [e for e in range(3, n - 1) if gcd(e, n) == 1][:2]    # identical to cross_relations.py
    A0, A1 = rand_unit_lu(n, rng), rand_unit_lu(n, rng)
    b0, b1 = rng.randrange(fld.q), rng.randrange(fld.q)
    pts = []
    for pt in range(fld.q):
        x0 = apply_mat(A0, pt, n) ^ b0
        x1 = apply_mat(A1, pt, n) ^ b1
        if x0 == 0 or x1 == 0:
            continue
        pts.append(sbox(fld, x0, es[0]) | (sbox(fld, x1, es[1]) << n))
    wt = [n - pow(e, -1, n) for e in es]
    print(f"# params: n={n} seed={seed} exps={es}")
    print(f"n={n} e={es}: Boolean degree of the image criterion Tr(t_i^-r_i) is wt(-r_i) = {wt} "
          f"(a cubic relation is possible iff some wt(-r_i) <= 3)")

    nv = 2 * n
    ms = []
    for d in range(4):
        ms.extend(itertools.combinations(range(nv), d))
    masks = [sum(1 << i for i in m) for m in ms]
    rows = [sum(1 << k for k, msk in enumerate(masks) if p & msk == msk) for p in pts]
    rk, _ = rank_gf2(rows)
    dim = len(ms) - rk
    print(f"   cubic relations on (t0,t1): dim = {dim} of {len(ms)} monomials on {len(pts)} points")

    # is the relation space contained in the span of monomials that touch only t1 (vars n..2n-1)?
    only_t1 = [k for k, m in enumerate(ms) if all(i >= n for i in m)]
    only_t0 = [k for k, m in enumerate(ms) if all(i < n for i in m)]
    # a relation lies in span(only_t1) iff adding the indicator of every other coordinate to the
    # evaluation matrix does not change the relation dimension: test by restricting the columns
    def dim_restricted(cols):
        sub = []
        cmask = 0
        for k in cols:
            cmask |= 1 << k
        for r in rows:
            sub.append(r & cmask)
        rk2, _ = rank_gf2(sub)
        return len(cols) - rk2
    d_t1 = dim_restricted(only_t1)
    d_t0 = dim_restricted(only_t0)
    print(f"   of which supported on t1 only: {d_t1};  on t0 only: {d_t0}")

    # does Tr(t1^{-r1}) actually vanish on the graph, and is it cubic?
    tt = anf_of_trace_criterion(fld, es[1], n)
    vanishes = all(tt[(p >> n) & (fld.q - 1)] == 0 for p in pts)
    print(f"   Tr(t1^-r_1) vanishes on every graph point: {vanishes} "
          f"(Boolean degree wt(-r_1) = {wt[1]})")
    tt0 = anf_of_trace_criterion(fld, es[0], n)
    vanishes0 = all(tt0[p & (fld.q - 1)] == 0 for p in pts)
    print(f"   Tr(t0^-r_0) vanishes on every graph point: {vanishes0} "
          f"(Boolean degree wt(-r_0) = {wt[0]})")
    verdict = ("all cubic relations are univariate image criteria — NO cross-S-box cubic"
               if dim == d_t0 + d_t1 else
               f"{dim - d_t0 - d_t1} cubic relation(s) are NOT univariate — genuine cross-S-box cubics")
    print(f"   => {verdict}")
    record("E12_cross_relations_explained.csv",
           dict(n=n, exps=";".join(map(str, es)), seed=seed, npoints=len(pts), monomials=len(ms),
                dim_cubic=dim, dim_t0_only=d_t0, dim_t1_only=d_t1,
                wt_minus_r0=wt[0], wt_minus_r1=wt[1],
                cross_cubics=dim - d_t0 - d_t1),
           seconds=time.time() - t_start)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 16,
         int(sys.argv[2]) if len(sys.argv) > 2 else 1)
