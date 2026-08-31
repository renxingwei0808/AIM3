#!/usr/bin/env python3
"""Symbolic kernel-variable bilinear model of toy AIM3 (l=2), Focus Report §4.3:
  variables U = (pt, t0, t1) [3n bits],  W = (a0, a1, a2) [3n bits]
  B_{e_i}(a_i, A_i pt + b_i) = a_i ,  t_i a_i = a_i^{2^{e_i}}   (i=0,1)
  B_{e_2}(a_2, b_2 + B_0 t_0 + B_1 t_1) = a_2 ,  (pt+ct) a_2 = a_2^{2^{e_2}}
6n Boolean equations, all bilinear in (U,W). Verifies vanishing at the key and reports Macaulay ranks / degree falls,
side by side with S_quad and P1-only on the same instance.
Usage: python3 kernel_bilinear_toy.py n Dmax        (n=5,7 in seconds; n=9 D<=3; larger -> server)
"""
import sys, os, random, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS
from math import gcd
from aim3_lab import F, sbox, rand_unit_lu, apply_mat
from aim3_macaulay import AffExpr, prod_coords, lin_coords, add_polys, build_aim3, macaulay, random_system

def build_kernel_model(n, exps, seed=1):
    ell = len(exps) - 1; assert ell == 2
    fld = F(n); rng = random.Random(seed)
    mats = [rand_unit_lu(n, rng) for _ in range(2 * ell)]
    consts = [rng.randrange(fld.q) for _ in range(ell + 1)]
    while True:
        pt = rng.randrange(1, fld.q)
        xs = [apply_mat(mats[i], pt, n) ^ consts[i] for i in range(ell)]
        if 0 in xs: continue
        ts = [sbox(fld, xs[i], exps[i]) for i in range(ell)]
        xl = consts[ell]
        for i in range(ell): xl ^= apply_mat(mats[ell + i], ts[i], n)
        if xl == 0 or 0 in ts: continue
        tl = sbox(fld, xl, exps[ell])
        if tl == 0: continue
        ct = pt ^ tl; break
    rs = [fld.pw(t, pow((1 << e) - 1, -1, fld.q - 1)) for t, e in zip(ts + [tl], exps)]   # kernel values a_i
    # variable layout: pt 0..n-1 | t0 n..2n-1 | t1 2n..3n-1 | a0 3n.. | a1 4n.. | a2 5n..
    nv = 6 * n
    sol = [0] * nv
    for k in range(n):
        sol[k] = (pt >> k) & 1
        for i in range(2): sol[(1 + i) * n + k] = (ts[i] >> k) & 1
        for i in range(3): sol[(3 + i) * n + k] = (rs[i] >> k) & 1
    def var(block): return AffExpr(fld, {block * n + k: 1 << k for k in range(n)}, 0)
    def colimgs(A): return [apply_mat(A, 1 << j, n) for j in range(n)]
    PT, T0, T1, A0, A1, A2 = (var(b) for b in range(6))
    eqs = []
    for i, (Ai, Ti) in enumerate(((A0, T0), (A1, T1))):
        ci = colimgs(mats[i]); e = exps[i]
        xi = AffExpr(fld, {j: ci[j] for j in range(n)}, consts[i])
        Bi = add_polys(prod_coords(Ai, xi.frob(e)), prod_coords(Ai.frob(e), xi), lin_coords(Ai))   # B(a_i,x_i)+a_i
        Ki = add_polys(prod_coords(Ti, Ai), lin_coords(Ai.frob(e)))                                    # t_i a_i + a_i^{2^e}
        eqs += Bi + Ki
    xl = AffExpr(fld, {}, consts[2])
    for i, Ti in enumerate((T0, T1)):
        ci = colimgs(mats[2 + i])
        for k in range(n):
            xl.coefs[(1 + i) * n + k] = xl.coefs.get((1 + i) * n + k, 0) ^ ci[k]
    e = exps[2]
    tl = AffExpr(fld, {k: 1 << k for k in range(n)}, ct)
    eqs += add_polys(prod_coords(A2, xl.frob(e)), prod_coords(A2.frob(e), xl), lin_coords(A2))
    eqs += add_polys(prod_coords(tl, A2), lin_coords(A2.frob(e)))
    def ev(p):
        s = 0
        for m in p:
            v = 1
            for j in m: v &= sol[j]
            s ^= v
        return s
    assert all(ev(p) == 0 for p in eqs), "kernel model does not vanish at the key"
    return nv, eqs, sol

if __name__ == "__main__":
    n = int(sys.argv[1]); Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    exps = [e for e in range(2, n - 1) if gcd(e, n) == 1][:3]
    nv, eqs, sol = build_kernel_model(n, exps)
    print(f"kernel bilinear model: n={n} exps={exps}: {len(eqs)} eqs in {nv} vars (6n) - vanishes at key: OK")
    nv2, eqsP1, eqsAll, _ = build_aim3(n, exps)
    for name, E, V in [("kernel-bilinear (6n eq / 6n var)", eqs, nv), ("S_quad (9n / 3n)", eqsAll, nv2), ("P1-only (3n / 3n)", eqsP1, nv2)]:
        print(f"  {name}")
        for D in range(2, Dmax + 1):
            rk, low, nrow, ncol = macaulay(E, V, D)
            print(f"     D={D}: rows={nrow:6d} cols={ncol:7d} rank={rk:6d} syz={nrow-rk:5d} dim(span∩deg<={D-1})={low}")
