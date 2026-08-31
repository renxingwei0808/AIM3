#!/usr/bin/env python3
"""
Boolean Macaulay-rank experiment for toy AIM3.

Builds the official quadratic system S_quad of AIM3 (variables x, t_0, ..., t_{l-1};
equations P1,P2,P3 per S-box) symbolically over F_2, and compares its Macaulay
ranks / degree falls against

  * the P1-only sub-system,
  * a random affine Boolean quadratic system with the same shape and a planted
    root.
"""
from __future__ import annotations
import itertools, random
from math import gcd
from aim3_lab import F, sbox, rand_unit_lu, apply_mat

# ---------------------------------------------------------- boolean polys
# monomial = frozenset of variable indices ; polynomial = set of monomials

def pmul_var(poly, v):
    return {m | {v} for m in poly}

def pxor(a, b):
    return a ^ b

# --------------------------------------- affine field expressions over F2

class AffExpr:
    """f = sum_j  var_j * coef_j  +  const   with coef_j, const in F_{2^n}."""
    def __init__(self, fld, coefs, const):
        self.fld = fld; self.coefs = dict(coefs); self.const = const

    def frob(self, k):
        f = self.fld
        return AffExpr(f, {j: f.pw(c, 1 << k) for j, c in self.coefs.items()},
                       f.pw(self.const, 1 << k))

def prod_coords(A: AffExpr, B: AffExpr):
    """coordinate polynomials of the field product A*B (list of n boolean polys)"""
    f = A.fld; n = f.n
    acc = {}                       # monomial -> field coefficient
    def add(mono, c):
        if c == 0: return
        acc[mono] = acc.get(mono, 0) ^ c
    for j, cj in A.coefs.items():
        for k, ck in B.coefs.items():
            add(frozenset((j, k)), f.mul(cj, ck))
        add(frozenset((j,)), f.mul(cj, B.const))
    for k, ck in B.coefs.items():
        add(frozenset((k,)), f.mul(A.const, ck))
    add(frozenset(), f.mul(A.const, B.const))
    out = [set() for _ in range(n)]
    for mono, c in acc.items():
        for i in range(n):
            if (c >> i) & 1:
                out[i] ^= {mono}
    return out

def lin_coords(A: AffExpr):
    f = A.fld; n = f.n
    acc = {}
    for j, cj in A.coefs.items():
        acc[frozenset((j,))] = acc.get(frozenset((j,)), 0) ^ cj
    acc[frozenset()] = acc.get(frozenset(), 0) ^ A.const
    out = [set() for _ in range(n)]
    for mono, c in acc.items():
        for i in range(n):
            if (c >> i) & 1:
                out[i] ^= {mono}
    return out

def add_polys(*lists):
    n = len(lists[0])
    out = [set() for _ in range(n)]
    for L in lists:
        for i in range(n):
            out[i] ^= L[i]
    return out

# ------------------------------------------------------- toy AIM3 system

def build_aim3(n, exps, seed=1):
    """returns (nvars, eqs_P1, eqs_all, solution_dict)"""
    ell = len(exps) - 1
    fld = F(n)
    rng = random.Random(seed)
    nv = (ell + 1) * n                       # x , t_0..t_{ell-1}
    # variable index: x bits 0..n-1 ; t_i bits (i+1)*n + k
    mats = [rand_unit_lu(n, rng) for _ in range(2 * ell)]
    consts = [rng.randrange(fld.q) for _ in range(ell + 1)]
    # find a valid key
    while True:
        pt = rng.randrange(1, fld.q)
        xs = [apply_mat(mats[i], pt, n) ^ consts[i] for i in range(ell)]
        if 0 in xs: continue
        ts = [sbox(fld, xs[i], exps[i]) for i in range(ell)]
        xl = consts[ell]
        for i in range(ell):
            xl ^= apply_mat(mats[ell + i], ts[i], n)
        if xl == 0: continue
        tl = sbox(fld, xl, exps[ell])
        ct = pt ^ tl
        break
    sol = [0] * nv
    for k in range(n): sol[k] = (pt >> k) & 1
    for i in range(ell):
        for k in range(n): sol[(i + 1) * n + k] = (ts[i] >> k) & 1

    # column images of the matrices as field elements
    def colimgs(A):
        return [apply_mat(A, 1 << j, n) for j in range(n)]

    eqsP1, eqsAll = [], []
    # first-round S-boxes
    for i in range(ell):
        ci = colimgs(mats[i])
        u = AffExpr(fld, {j: ci[j] for j in range(n)}, consts[i])   # in x-bits
        t = AffExpr(fld, {(i + 1) * n + k: 1 << k for k in range(n)}, 0)
        e = exps[i]
        P1 = add_polys(prod_coords(u, t), lin_coords(u.frob(e)),
                       lin_coords(AffExpr(fld, {}, 1)))
        u2 = u.frob(1)
        P2 = add_polys(prod_coords(u2, t), prod_coords(u.frob(e), u), lin_coords(u))
        # P3 : u t^2 + u^{2^e} t + t
        t2 = AffExpr(fld, {(i + 1) * n + k: fld.mul(1 << k, 1 << k) for k in range(n)}, 0)
        # careful: t^2 as field element is F2-linear in t bits with coef beta_k^2
        t2 = AffExpr(fld, {(i + 1) * n + k: fld.pw(1 << k, 2) for k in range(n)}, 0)
        P3 = add_polys(prod_coords(u, t2), prod_coords(u.frob(e), t), lin_coords(t))
        eqsP1 += P1
        eqsAll += P1 + P2 + P3
    # last S-box
    ul = AffExpr(fld, {}, consts[ell])
    for i in range(ell):
        ci = colimgs(mats[ell + i])
        for k in range(n):
            ul.coefs[(i + 1) * n + k] = ul.coefs.get((i + 1) * n + k, 0) ^ ci[k]
    tl = AffExpr(fld, {k: 1 << k for k in range(n)}, ct)   # t_l = x + ct
    tl2 = AffExpr(fld, {k: fld.pw(1 << k, 2) for k in range(n)}, fld.pw(ct, 2))
    e = exps[ell]
    P1 = add_polys(prod_coords(ul, tl), lin_coords(ul.frob(e)),
                   lin_coords(AffExpr(fld, {}, 1)))
    P2 = add_polys(prod_coords(ul.frob(1), tl), prod_coords(ul.frob(e), ul), lin_coords(ul))
    P3 = add_polys(prod_coords(ul, tl2), prod_coords(ul.frob(e), tl), lin_coords(tl))
    eqsP1 += P1
    eqsAll += P1 + P2 + P3

    # sanity: all equations vanish at the solution
    def ev(poly):
        s = 0
        for m in poly:
            v = 1
            for j in m: v &= sol[j]
            s ^= v
        return s
    assert all(ev(p) == 0 for p in eqsAll), "system does not vanish at the key"
    return nv, eqsP1, eqsAll, sol


def random_system(nv, neq, seed=2):
    rng = random.Random(seed)
    sol = [rng.getrandbits(1) for _ in range(nv)]
    eqs = []
    for _ in range(neq):
        p = set()
        for i in range(nv):
            for j in range(i + 1, nv):
                if rng.getrandbits(1): p ^= {frozenset((i, j))}
        for i in range(nv):
            if rng.getrandbits(1): p ^= {frozenset((i,))}
        # force vanishing at sol
        s = 0
        for m in p:
            v = 1
            for j in m: v &= sol[j]
            s ^= v
        if s: p ^= {frozenset()}
        eqs.append(p)
    return eqs, sol

# ------------------------------------------------------ Macaulay machinery

def monomial_index(nv, D):
    idx = {}; by_deg = [[] for _ in range(D + 1)]
    c = 0
    for d in range(D + 1):
        for m in itertools.combinations(range(nv), d):
            idx[frozenset(m)] = c; by_deg[d].append(c); c += 1
    return idx, c, by_deg

def poly_to_row(poly, idx):
    r = 0
    for m in poly:
        r |= 1 << idx[m]
    return r

def rank_bits(rows):
    piv = {}; rk = 0
    for r in rows:
        while r:
            b = r.bit_length() - 1
            if b in piv: r ^= piv[b]
            else:
                piv[b] = r; rk += 1; break
    return rk

def macaulay(eqs, nv, D):
    """returns (rank_D, dim of (rowspace ∩ deg<=D-1))"""
    idx, ncol, by_deg = monomial_index(nv, D)
    topcols = set(by_deg[D])
    mults = [frozenset(m) for d in range(D - 1)
             for m in itertools.combinations(range(nv), d)]
    rows, rows_top = [], []
    topmask = 0
    for c in topcols: topmask |= 1 << c
    for f in eqs:
        for mu in mults:
            g = set()
            for m in f:
                g ^= {m | mu}
            if not g: continue
            if max(len(m) for m in g) > D: continue
            r = poly_to_row(g, idx)
            rows.append(r)
    rk = rank_bits(rows)
    rk_top = rank_bits([r & topmask for r in rows])
    return rk, rk - rk_top, len(rows), ncol


def report(name, eqs, nv, Dmax):
    print(f"  {name}: {len(eqs)} eqs, {nv} vars")
    prev = None
    for D in range(2, Dmax + 1):
        rk, low, nrow, ncol = macaulay(eqs, nv, D)
        # dimension of degree<=D-1 part obtained one degree earlier
        note = ""
        if prev is not None and low > prev:
            note = f"   <-- degree fall: {low - prev} new elts of deg<={D-1}"
        print(f"     D={D}: rows={nrow:6d} cols={ncol:6d} rank={rk:6d} "
              f"dim(span ∩ deg<={D-1})={low:6d}{note}")
        prev = low if prev is None else max(prev, low)
    print()


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    cand = [e for e in range(2, n - 1) if gcd(e, n) == 1]
    exps = cand[:3] if len(cand) >= 3 else cand
    print(f"toy AIM3: n={n}, l=2, exponents={exps}, {3*n} variables")
    nv, eqsP1, eqsAll, sol = build_aim3(n, exps)
    print("=" * 74)
    report("AIM3  P1+P2+P3 (official S_quad)", eqsAll, nv, Dmax)
    report("AIM3  P1 only (bilinear core)   ", eqsP1, nv, Dmax)
    reqs, _ = random_system(nv, len(eqsAll))
    report("random quadratic, same shape    ", reqs, nv, Dmax)
    reqs2, _ = random_system(nv, len(eqsP1), seed=5)
    report("random quadratic, P1 shape      ", reqs2, nv, Dmax)
