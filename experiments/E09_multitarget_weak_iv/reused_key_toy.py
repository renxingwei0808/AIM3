#!/usr/bin/env python3
"""E09 — reused-key toy (track T6, attacker model M-reuse): the SAME pt under q different iv.

ROADMAP section 3.6: "q pairs (iv_j, ct_j) sharing pt -> 3nq equations, n + 2nq variables; measure
d_solv for q = 2, 3."  This is the Tier-A part: the D=3 profile.

Layout: pt occupies bits 0..n-1 and is SHARED; instance j contributes its own S-box outputs
t_0^{(j)}..t_{l-1}^{(j)} in the block starting at n + (j*l + i)*n.  Each instance has its own
independent unit-LU matrices and constants (a different iv), and its own ct.  So

    variables  n + q*l*n ,   equations  q * 3*(l+1)*n .

The question is whether sharing pt couples the instances algebraically: if it does not, the D=3
syzygy count should be exactly q*(l+1)*8n -- the direct sum of q independent systems -- and the
number of new elements should be q times a single instance's.  Anything above that is a coupling
gain for the attacker.

Saturation is reported and honoured: the D=3 row space of a heavily over-determined system can fill
the whole degree-<=2 space, in which case the numbers say nothing (same trap as E04 at D=4).

Usage: python3 reused_key_toy.py N [-l ELL] [-q Q] [--seeds 1,2] [--budget SEC]
Appends to results/E09_theory_reused_key.csv.
"""
import sys, os, time, random, argparse
from math import gcd, comb
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3_lab import F, sbox, rand_unit_lu, apply_mat
from aim3_macaulay import AffExpr, prod_coords, lin_coords, add_polys
from aim3.boolmac import macaulay_stats, MacaulayTimeout


def build_reused(n, exps, q, seed=1):
    """q AIM3 instances (l = len(exps)-1) sharing one pt.  Returns (nv, eqs, sol)."""
    ell = len(exps) - 1
    fld = F(n)
    rng = random.Random(seed)
    nv = n + q * ell * n

    def block(j, i):
        return n + (j * ell + i) * n

    # one pt that is a valid key for every instance
    tries = 0
    while True:
        tries += 1
        if tries > 500:
            raise RuntimeError("no shared pt valid for all instances")
        rng2 = random.Random(seed * 7919 + tries)
        mats = [[rand_unit_lu(n, rng2) for _ in range(2 * ell)] for _ in range(q)]
        consts = [[rng2.randrange(fld.q) for _ in range(ell + 1)] for _ in range(q)]
        pt = rng2.randrange(1, fld.q)
        ok = True
        ts, cts = [], []
        for j in range(q):
            xs = [apply_mat(mats[j][i], pt, n) ^ consts[j][i] for i in range(ell)]
            if 0 in xs:
                ok = False
                break
            tj = [sbox(fld, xs[i], exps[i]) for i in range(ell)]
            xl = consts[j][ell]
            for i in range(ell):
                xl ^= apply_mat(mats[j][ell + i], tj[i], n)
            if xl == 0:
                ok = False
                break
            tl = sbox(fld, xl, exps[ell])
            ts.append(tj)
            cts.append(pt ^ tl)
        if ok:
            break

    sol = [0] * nv
    for k in range(n):
        sol[k] = (pt >> k) & 1
    for j in range(q):
        for i in range(ell):
            for k in range(n):
                sol[block(j, i) + k] = (ts[j][i] >> k) & 1

    def colimgs(M):
        return [apply_mat(M, 1 << j, n) for j in range(n)]

    eqs = []
    for j in range(q):
        for i in range(ell):
            ci = colimgs(mats[j][i])
            u = AffExpr(fld, {k: ci[k] for k in range(n)}, consts[j][i])       # x_i = A_i pt + b_i
            t = AffExpr(fld, {block(j, i) + k: 1 << k for k in range(n)}, 0)
            t2 = AffExpr(fld, {block(j, i) + k: fld.pw(1 << k, 2) for k in range(n)}, 0)
            e = exps[i]
            eqs += add_polys(prod_coords(u, t), lin_coords(u.frob(e)), lin_coords(AffExpr(fld, {}, 1)))
            eqs += add_polys(prod_coords(u.frob(1), t), prod_coords(u.frob(e), u), lin_coords(u))
            eqs += add_polys(prod_coords(u, t2), prod_coords(u.frob(e), t), lin_coords(t))
        ul = AffExpr(fld, {}, consts[j][ell])
        for i in range(ell):
            ci = colimgs(mats[j][ell + i])
            for k in range(n):
                ul.coefs[block(j, i) + k] = ul.coefs.get(block(j, i) + k, 0) ^ ci[k]
        tl = AffExpr(fld, {k: 1 << k for k in range(n)}, cts[j])                # t_l = pt + ct_j
        tl2 = AffExpr(fld, {k: fld.pw(1 << k, 2) for k in range(n)}, fld.pw(cts[j], 2))
        e = exps[ell]
        eqs += add_polys(prod_coords(ul, tl), lin_coords(ul.frob(e)), lin_coords(AffExpr(fld, {}, 1)))
        eqs += add_polys(prod_coords(ul.frob(1), tl), prod_coords(ul.frob(e), ul), lin_coords(ul))
        eqs += add_polys(prod_coords(ul, tl2), prod_coords(ul.frob(e), tl), lin_coords(tl))

    for p in eqs:
        s = 0
        for m in p:
            v = 1
            for k in m:
                v &= sol[k]
            s ^= v
        assert s == 0, "reused-key system does not vanish at the shared key"
    return nv, eqs, sol


def profile(n, exps, q, seed, budget):
    t0 = time.time()
    ell = len(exps) - 1
    nv, eqs, _ = build_reused(n, exps, q, seed)
    low_cols = sum(comb(nv, d) for d in range(3))
    st2 = macaulay_stats(eqs, nv, 2)
    st3 = macaulay_stats(eqs, nv, 3, budget=budget)
    new = st3.dim_low - st2.rank
    sat = st3.dim_low / low_cols
    target = q * (ell + 1) * 8 * n
    return dict(n=n, ell=ell, q=q, seed=seed, exps=";".join(map(str, exps)), nvars=nv, neqs=len(eqs),
                mac3_rows=st3.rows, mac3_rank=st3.rank, mac3_syz=st3.syz, syz_direct_sum=target,
                syz_excess=st3.syz - target, mac3_dim_low=st3.dim_low, new_at_D3=new,
                low_cols=low_cols, saturation=round(sat, 4),
                verdict=("saturated" if sat >= 0.99 else
                         ("direct-sum" if st3.syz == target else "coupled")),
                seconds=round(time.time() - t0, 3))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("-l", "--ell", type=int, default=1)
    ap.add_argument("-q", "--q", type=int, default=2)
    ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--budget", type=float, default=3600.0)
    a = ap.parse_args()
    exps = [e for e in range(2, a.n) if gcd(e, a.n) == 1][:a.ell + 1]
    print(f"# params: n={a.n} ell={a.ell} q={a.q} seeds={a.seeds} exps={exps}")
    for seed in [int(s) for s in a.seeds.split(",")]:
        single = profile(a.n, exps, 1, seed, a.budget)
        multi = profile(a.n, exps, a.q, seed, a.budget)
        print(f"  seed {seed}:")
        for tag, r in (("q=1 (single)", single), (f"q={a.q} (shared pt)", multi)):
            print(f"    {tag:>18}: {r['nvars']:3d} vars, {r['neqs']:4d} eqs; D=3 rank {r['mac3_rank']:6d} "
                  f"syz {r['mac3_syz']:5d} (direct sum {r['syz_direct_sum']}) "
                  f"dim_low {r['mac3_dim_low']:5d}/{r['low_cols']:5d} sat {r['saturation']:.3f} "
                  f"new {r['new_at_D3']:5d} -> {r['verdict']}")
            record("E09_theory_reused_key.csv", r, seconds=r.pop("seconds"))
        if multi["saturation"] < 0.99:
            print(f"    => syzygy excess over the direct sum of {a.q} independent instances: "
                  f"{multi['syz_excess']};  new elements {multi['new_at_D3']} vs "
                  f"{a.q} x {single['new_at_D3']} = {a.q * single['new_at_D3']} "
                  f"(difference {multi['new_at_D3'] - a.q * single['new_at_D3']})")
        else:
            print(f"    => q={a.q} run is SATURATED (dim_low/C(nv,<=2) = {multi['saturation']:.4f}); "
                  f"the D=3 numbers are not informative here")
