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

Usage: python3 explain_n16.py [n] [seed] [--deg D] [--backend python|sage]
`--deg 4` runs the same decomposition one degree higher (batch-1 task 5: subtract the univariate
image criteria before reporting any cross-S-box quartic).  `--backend sage` computes the three
ranks with M4RI -- at n=16, deg=4 the evaluation matrix is 65534 x 41449, out of reach in pure
Python.  Both engines are verified to agree by server/jobs/verify_backend.py.
"""
import sys, os, random, time
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3_lab import F, sbox, rand_unit_lu, apply_mat, rank_gf2
from aim3.boolmac_export import rank_sparse_sage
import itertools

SCRATCH = os.environ.get("AIM3_SCRATCH")


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


def main(n, seed=1, deg=3, backend="python"):
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
    print(f"# params: n={n} seed={seed} exps={es} deg={deg} backend={backend}")
    print(f"n={n} e={es}: Boolean degree of the image criterion Tr(t_i^-r_i) is wt(-r_i) = {wt} "
          f"(a degree-{deg} relation is possible iff some wt(-r_i) <= {deg})")

    nv = 2 * n
    ms = []
    for d in range(deg + 1):
        ms.extend(itertools.combinations(range(nv), d))
    masks = [sum(1 << i for i in m) for m in ms]
    only_t1 = [k for k, m in enumerate(ms) if all(i >= n for i in m)]
    only_t0 = [k for k, m in enumerate(ms) if all(i < n for i in m)]

    if backend == "sage":
        # rows as column-index lists (the nonzeros of a point row are the sub-monomials of its
        # support), then one M4RI rank per column subset
        def idx_rows(cols=None):
            keep = None if cols is None else {k: i for i, k in enumerate(cols)}
            for p in pts:
                supp = [i for i in range(nv) if (p >> i) & 1]
                out = []
                for d in range(min(deg, len(supp)) + 1):
                    for c in itertools.combinations(supp, d):
                        k = KIDX[sum(1 << i for i in c)]
                        if keep is None:
                            out.append(k)
                        elif k in keep:
                            out.append(keep[k])
                yield out
        KIDX = {msk: k for k, msk in enumerate(masks)}
        rk = rank_sparse_sage(idx_rows(), len(ms), workdir=SCRATCH, tag=f"exp_n{n}_d{deg}_all")
        dim = len(ms) - rk
        print(f"   degree-<={deg} relations on (t0,t1): dim = {dim} of {len(ms)} monomials "
              f"on {len(pts)} points")
        d_t1 = len(only_t1) - rank_sparse_sage(idx_rows(only_t1), len(only_t1), workdir=SCRATCH,
                                               tag=f"exp_n{n}_d{deg}_t1")
        d_t0 = len(only_t0) - rank_sparse_sage(idx_rows(only_t0), len(only_t0), workdir=SCRATCH,
                                               tag=f"exp_n{n}_d{deg}_t0")
    else:
        rows = [sum(1 << k for k, msk in enumerate(masks) if p & msk == msk) for p in pts]
        rk, _ = rank_gf2(rows)
        dim = len(ms) - rk
        print(f"   degree-<={deg} relations on (t0,t1): dim = {dim} of {len(ms)} monomials "
              f"on {len(pts)} points")

        # a relation lies in span(only_t1) iff restricting the evaluation matrix to those columns
        # does not change the relation dimension
        def dim_restricted(cols):
            cmask = 0
            for k in cols:
                cmask |= 1 << k
            rk2, _ = rank_gf2([r & cmask for r in rows])
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
    verdict = (f"all degree-<={deg} relations are univariate image criteria — NO cross-S-box "
               f"relation of degree <= {deg}"
               if dim == d_t0 + d_t1 else
               f"{dim - d_t0 - d_t1} relation(s) of degree <= {deg} are NOT univariate — genuine "
               f"cross-S-box relations")
    print(f"   => {verdict}")
    record("E12_cross_relations_explained.csv",
           dict(n=n, exps=";".join(map(str, es)), seed=seed, deg=deg, npoints=len(pts),
                monomials=len(ms), dim_cubic=dim, dim_t0_only=d_t0, dim_t1_only=d_t1,
                wt_minus_r0=wt[0], wt_minus_r1=wt[1],
                cross_cubics=dim - d_t0 - d_t1),
           seconds=time.time() - t_start)


if __name__ == "__main__":
    _deg = int(sys.argv[sys.argv.index("--deg") + 1]) if "--deg" in sys.argv else 3
    _backend = sys.argv[sys.argv.index("--backend") + 1] if "--backend" in sys.argv else "python"
    _pos = [a for a in sys.argv[1:] if a.isdigit()]
    if "--deg" in sys.argv:
        _pos = [a for a in _pos if a != sys.argv[sys.argv.index("--deg") + 1]] or _pos
    main(int(_pos[0]) if _pos else 16, int(_pos[1]) if len(_pos) > 1 else 1, _deg, _backend)
