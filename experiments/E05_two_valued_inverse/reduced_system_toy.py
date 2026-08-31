#!/usr/bin/env python3
"""E05 (Tier A) — explicit two-valued-inverse REDUCED-VARIABLE system for toy AIM3 (l=2).

Track T3 / spec claim C2.  Two field unknowns (t0,t1) instead of the official (l+1)n = 3n bits.

Pipeline (standard library only; field arithmetic from src/aim3/gf2n.py):

 1. H_e as a 2-linearised polynomial.  `SboxStructures.H` is an F_2-linear section of
    L(s)=s^{2^e}+s on the trace-zero hyperplane.  We pick a basis {b_k} of F_{2^n}/F_2 made of
    n-1 trace-zero vectors plus one trace-one vector g, set H(g):=0 (this pins the one degree of
    freedom: linearised polys vanishing on the trace-zero hyperplane are exactly c*Tr(u)), and
    solve the n x n Moore system  sum_j h_j b_k^{2^j} = H(b_k)  by Gaussian elimination over
    GF(2^n).  Checked on 100 random trace-zero u.

 2. G_{e,eps}(t) = sum_j h_j t^{r(1-2^j)} + eps t^r + 1   as a sparse exponent->coefficient dict.
    Checked: for random x != 1 and t = S_e(x),  {G_{e,0}(t), G_{e,1}(t)} = {x, x+a}, a = t^r.

 3. The reduced system in F_{2^n}[t0,t1].  With x_i = A_i pt + b_i, t_i = S_{e_i}(x_i),
    x_2 = b_2 + B_0 t_0 + B_1 t_1, t_2 = S_{e_2}(x_2), ct = pt + t_2:

      Eq1:  (A_1 A_0^-1)(G_{e0,eps0}(t0) + b_0) + b_1 + G_{e1,eps1}(t1) = 0
      Eq2:  x_2 * t_2 + x_2^{2^{e_2}} + 1 = 0,
            with  t_2 = A_0^-1(G_{e0,eps0}(t0) + b_0) + ct   and  x_2 = b_2 + B_0 t_0 + B_1 t_1.

    Only "linearised polynomial o sparse polynomial", sparse products and additions are used, so
    both equations stay sparse.  IMPORTANT (reported, not hidden): the ROADMAP form of the second
    equation,  b_2 + B_0 t0 + B_1 t1 = G_{e2,eps2}(A_0^-1(G_{e0,eps0}(t0)+b_0) + ct),  CANNOT be
    expanded into a sparse polynomial -- it is a non-linearised sparse polynomial composed with a
    non-linearised sparse polynomial, whose expansion has ~2^n terms.  Eq2 above is the same
    variety written through the forward S-box relation P1 of the last S-box; it is equivalent for
    x_2 != 0 and it removes the third branch bit, so the system carries 2 branch bits, not 3.
    The ROADMAP form is still verified POINTWISE at the true key for the correct (eps0,eps1,eps2).

 4. Statistics per equation: #terms, ordinary degree in t0 and in t1, Boolean-degree (Hamming
    weight) distribution of the terms, and the Bezout bound deg_t0 * deg_t1.  Plus the exact number
    of roots of the sparse system per branch (see `count_roots`) -- the pseudo-solution count that
    COMPLEXITY_ACCOUNTING section 5 requires.  No resultant here (Tier C, see README).

Usage: python3 reduced_system_toy.py n [seed] [--no-roots]
       python3 reduced_system_toy.py --sweep 8,9,10,11,12,13,14 [seed]
"""
import sys, os, time, random
from math import gcd
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3.gf2n import GF2n
from aim3.structures import SboxStructures
from aim3_lab import rand_unit_lu, apply_mat


# ---------------------------------------------------------------- GF(2^n) linear algebra
def solve_gf2n(F, M, v):
    """Solve M h = v over GF(2^n) by Gaussian elimination. M is a list of rows (lists), v a list."""
    n = len(M)
    A = [list(row) + [v[i]] for i, row in enumerate(M)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c]), None)
        if p is None:
            raise ValueError("singular Moore matrix (the b_k are not an F_2-basis)")
        A[c], A[p] = A[p], A[c]
        ic = F.inv(A[c][c])
        A[c] = [F.mul(ic, x) for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [x ^ F.mul(f, y) for x, y in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]


def linearise(F, basis, values):
    """coefficients h_j of the linearised polynomial sum_j h_j u^{2^j} agreeing with
    basis[k] |-> values[k] (Moore matrix M[k][j] = basis[k]^{2^j})."""
    n = F.n
    M = [[F.frob(b, j) for j in range(n)] for b in basis]
    return solve_gf2n(F, M, list(values))


def linear_map_coeffs(F, mat):
    """linearised-polynomial coefficients of the F_2-linear map given by a bit matrix."""
    n = F.n
    basis = [1 << k for k in range(n)]
    return linearise(F, basis, [apply_mat(mat, b, n) for b in basis])


def invert_bitmat(n, mat):
    """inverse of an n x n bit matrix (rows are ints in the apply_mat convention: y_i = sum_j A[i][j] x_j)."""
    A = [sum(mat[i][j] << j for j in range(n)) | (1 << (n + i)) for i in range(n)]
    piv = {}
    for i in range(n):
        r = A[i]
        for c in range(n):
            if (r >> c) & 1:
                if c in piv:
                    r ^= piv[c]
                else:
                    piv[c] = r
                    break
        else:
            raise ValueError("singular matrix")
    # back-substitute to reduced form
    for c in sorted(piv, reverse=True):
        for c2 in sorted(piv):
            if c2 != c and (piv[c2] >> c) & 1:
                piv[c2] ^= piv[c]
    inv = [[0] * n for _ in range(n)]
    for c in range(n):                      # after Gauss-Jordan, piv[c] = e_c | (A^-1 row c) << n
        row = piv[c] >> n
        for k in range(n):
            inv[c][k] = (row >> k) & 1
    for k in range(n):                      # self-check A * A^-1 = I on the basis
        assert apply_mat(mat, apply_mat(inv, 1 << k, n), n) == 1 << k, "matrix inversion failed"
    return inv


def eval_lin(F, coeffs, u):
    """sum_j coeffs[j] u^{2^j}"""
    s = 0
    for j, c in enumerate(coeffs):
        if c:
            s ^= F.mul(c, F.frob(u, j))
    return s


# ---------------------------------------------------------------- sparse F_{2^n}[t0,t1]
class SPoly:
    """sparse polynomial: {(i, j): coefficient}, i = exponent of t0, j of t1, taken mod 2^n-1."""
    __slots__ = ("F", "t", "q1")

    def __init__(self, F, terms=None):
        self.F, self.q1 = F, F.q - 1
        self.t = {}
        if terms:
            for k, c in terms.items():
                self.add_term(k, c)

    def add_term(self, key, c):
        if not c:
            return
        i, j = key[0] % self.q1, key[1] % self.q1
        v = self.t.get((i, j), 0) ^ c
        if v:
            self.t[(i, j)] = v
        else:
            self.t.pop((i, j), None)

    def __add__(self, o):
        r = SPoly(self.F, self.t)
        for k, c in o.t.items():
            r.add_term(k, c)
        return r

    def __mul__(self, o):
        r = SPoly(self.F)
        for (i1, j1), c1 in self.t.items():
            for (i2, j2), c2 in o.t.items():
                r.add_term((i1 + i2, j1 + j2), self.F.mul(c1, c2))
        return r

    def frob(self, k):
        """self^{2^k}: exponents doubled k times, coefficients Frobenius-twisted."""
        r = SPoly(self.F)
        f = 1 << k
        for (i, j), c in self.t.items():
            r.add_term((i * f, j * f), self.F.frob(c, k))
        return r

    def compose_linear(self, coeffs):
        """sum_j coeffs[j] * self^{2^j}  -- linearised polynomial applied to a sparse polynomial."""
        r = SPoly(self.F)
        for j, c in enumerate(coeffs):
            if not c:
                continue
            for (a, b), cc in self.frob(j).t.items():
                r.add_term((a, b), self.F.mul(c, cc))
        return r

    def eval(self, t0, t1):
        s = 0
        for (i, j), c in self.t.items():
            v = c
            if i:
                v = self.F.mul(v, self.F.pow(t0, i))
            if j:
                v = self.F.mul(v, self.F.pow(t1, j))
            s ^= v
        return s

    def stats(self):
        di = max((i for i, _ in self.t), default=0)
        dj = max((j for _, j in self.t), default=0)
        bool_deg = Counter(bin(i).count("1") + bin(j).count("1") for i, j in self.t)
        return dict(terms=len(self.t), deg_t0=di, deg_t1=dj,
                    bool_deg_max=max(bool_deg, default=0), bool_deg_hist=bool_deg)


def const(F, c):
    return SPoly(F, {(0, 0): c})


def var(F, which):
    return SPoly(F, {(1, 0) if which == 0 else (0, 1): 1})


# ---------------------------------------------------------------- H_e and G_{e,eps}
def H_linearised(F, S):
    """linearised coefficients of the section H_e, using a basis of n-1 trace-zero vectors + one
    trace-one vector g with H(g) := 0."""
    n = F.n
    raw = [1 << k for k in range(n)]
    g = next((b for b in raw if F.trace(b) == 1), None)
    if g is None:                                  # cannot happen: Tr is onto
        raise RuntimeError("no trace-one basis vector")
    basis = [g] + [(b ^ g if F.trace(b) == 1 else b) for b in raw if b != g]
    values = [0] + [S.H(b) for b in basis[1:]]
    return linearise(F, basis, values)


def check_H(F, S, h, rng, trials=100):
    for _ in range(trials):
        while True:
            u = rng.randrange(1, F.q)
            if F.trace(u) == 0:
                break
        assert eval_lin(F, h, u) == S.H(u), f"H mismatch at u={u}"
    return True


def G_poly(F, S, h, eps, which):
    """G_{e,eps} as a sparse polynomial in t0 (which=0) or t1 (which=1)."""
    r = S.r
    P = SPoly(F)
    for j, c in enumerate(h):
        k = (r * (1 - (1 << j)))
        P.add_term((k, 0) if which == 0 else (0, k), c)
    if eps:
        P.add_term((r, 0) if which == 0 else (0, r), 1)
    P.add_term((0, 0), 1)
    return P


def check_G(F, S, h, rng, trials=100):
    """{G_{e,0}(t), G_{e,1}(t)} = {x, x+a} for t = S_e(x), x != 1."""
    e = S.e
    G0, G1 = G_poly(F, S, h, 0, 0), G_poly(F, S, h, 1, 0)
    for _ in range(trials):
        x = rng.randrange(1, F.q)
        if x == 1:
            continue
        t = F.sbox(x, e)
        a = F.pow(t, S.r)
        assert {G0.eval(t, 0), G1.eval(t, 0)} == {x, x ^ a}, f"G mismatch at x={x}"
    return True


# ---------------------------------------------------------------- the toy AIM3 instance
def build_instance(n, exps, seed):
    """legacy unit-LU matrices + a valid key, exactly the shape of legacy/aim3_macaulay.build_aim3."""
    F = GF2n(n)
    rng = random.Random(seed)
    mats = [rand_unit_lu(n, rng) for _ in range(4)]          # A0, A1, B0, B1
    consts = [rng.randrange(F.q) for _ in range(3)]          # b0, b1, b2
    while True:
        pt = rng.randrange(1, F.q)
        xs = [apply_mat(mats[i], pt, n) ^ consts[i] for i in range(2)]
        if 0 in xs:
            continue
        ts = [F.sbox(xs[i], exps[i]) for i in range(2)]
        if 0 in ts or 1 in xs:
            continue
        x2 = consts[2] ^ apply_mat(mats[2], ts[0], n) ^ apply_mat(mats[3], ts[1], n)
        if x2 == 0:
            continue
        t2 = F.sbox(x2, exps[2])
        if t2 == 0:
            continue
        return F, mats, consts, pt, xs, ts, x2, t2, pt ^ t2


def build_reduced_system(n, exps, seed, rng):
    F, mats, consts, pt, xs, ts, x2, t2, ct = build_instance(n, exps, seed)
    S = [SboxStructures(n, e, F) for e in exps]
    h = [H_linearised(F, s) for s in S]
    for s, hh in zip(S, h):
        check_H(F, s, hh, rng)
        check_G(F, s, hh, rng)

    A0inv = invert_bitmat(n, mats[0])
    L_A0inv = linear_map_coeffs(F, A0inv)
    L_A1A0inv = linearise(F, [1 << k for k in range(n)],
                          [apply_mat(mats[1], apply_mat(A0inv, 1 << k, n), n) for k in range(n)])
    L_B0, L_B1 = linear_map_coeffs(F, mats[2]), linear_map_coeffs(F, mats[3])

    # which branch bit is the true one?
    eps = []
    for i in (0, 1):
        G0 = G_poly(F, S[i], h[i], 0, i)
        eps.append(0 if G0.eval(ts[0] if i == 0 else 0, ts[1] if i == 1 else 0) == xs[i] else 1)
    eps2 = 0 if G_poly(F, S[2], h[2], 0, 0).eval(t2, 0) == x2 else 1

    G0 = G_poly(F, S[0], h[0], eps[0], 0)
    G1 = G_poly(F, S[1], h[1], eps[1], 1)
    inner = G0 + const(F, consts[0])                      # G_{e0,eps0}(t0) + b0   (= pt image under A0)
    eq1 = inner.compose_linear(L_A1A0inv) + const(F, consts[1]) + G1
    pt_poly = inner.compose_linear(L_A0inv)               # = pt
    t2_poly = pt_poly + const(F, ct)                      # = t2
    x2_poly = (var(F, 0).compose_linear(L_B0) + var(F, 1).compose_linear(L_B1) + const(F, consts[2]))
    eq2 = x2_poly * t2_poly + x2_poly.frob(exps[2]) + const(F, 1)
    ctx = dict(F=F, S=S, h=h, mats=mats, consts=consts, pt=pt, xs=xs, ts=ts, x2=x2, t2=t2, ct=ct,
               eps=eps, eps2=eps2, A0inv=A0inv, L_A0inv=L_A0inv, L_A1A0inv=L_A1A0inv,
               pt_poly=pt_poly, t2_poly=t2_poly, x2_poly=x2_poly, G0=G0, G1=G1)
    return eq1, eq2, ctx


def verify(eq1, eq2, ctx):
    """the true key must be a root of the sparse system, and of the ROADMAP form pointwise."""
    F, ts = ctx["F"], ctx["ts"]
    t0, t1 = ts
    assert eq1.eval(t0, t1) == 0, "Eq1 does not vanish at the true key"
    assert eq2.eval(t0, t1) == 0, "Eq2 does not vanish at the true key"
    assert ctx["pt_poly"].eval(t0, t1) == ctx["pt"]
    assert ctx["t2_poly"].eval(t0, t1) == ctx["t2"]
    assert ctx["x2_poly"].eval(t0, t1) == ctx["x2"]
    # ROADMAP form of the second equation, evaluated (it cannot be expanded, see the module docstring)
    G2 = G_poly(F, ctx["S"][2], ctx["h"][2], ctx["eps2"], 0)
    assert G2.eval(ctx["t2"], 0) == ctx["x2"], "ROADMAP Eq2 fails at the true key"
    return True


def count_roots(eq1, eq2, ctx):
    """exact number of (t0,t1) in (F*)^2 satisfying BOTH sparse equations, for the true branch.

    Eq1 is  W(t0) = V(t1)  with W(t0) = (A1 A0^-1)(G0(t0)+b0) + b1 and V(t1) = G1(t1); both are
    tabulated in O(2^n), so the whole Eq1 variety is enumerated without a 2^{2n} scan.  Each of its
    points is then tested against Eq2 by field arithmetic.
    """
    F, ctx_G0, ctx_G1 = ctx["F"], ctx["G0"], ctx["G1"]
    n = F.n
    by_v = {}
    for t1 in range(1, F.q):
        by_v.setdefault(ctx_G1.eval(0, t1), []).append(t1)
    L, b1 = ctx["L_A1A0inv"], ctx["consts"][1]
    L0, ct, b2 = ctx["L_A0inv"], ctx["ct"], ctx["consts"][2]
    mats, e2 = ctx["mats"], ctx["S"][2].e
    roots = []
    for t0 in range(1, F.q):
        w = eval_lin(F, L, ctx_G0.eval(t0, 0) ^ ctx["consts"][0]) ^ b1
        cands = by_v.get(w)
        if not cands:
            continue
        t2v = eval_lin(F, L0, ctx_G0.eval(t0, 0) ^ ctx["consts"][0]) ^ ct
        base = b2 ^ apply_mat(mats[2], t0, n)
        for t1 in cands:
            x2v = base ^ apply_mat(mats[3], t1, n)
            if F.mul(x2v, t2v) ^ F.frob(x2v, e2) ^ 1 == 0:
                roots.append((t0, t1))
    return roots


def run(n, seed=1, do_roots=True):
    t0 = time.time()
    rng = random.Random(1000 + seed)
    # e ranges over 2..n-1 (the spec only requires gcd(e,n)=1); e=n-1 is needed to get three
    # distinct admissible exponents at n=8,10,12.  E03 uses 2..n-2, so the two exponent lists differ
    # only at those n -- stated here so the CSVs are not compared blindly.
    exps = [e for e in range(2, n) if gcd(e, n) == 1][:3]
    if len(exps) < 3:
        raise SystemExit(f"n={n}: fewer than 3 admissible exponents")
    eq1, eq2, ctx = build_reduced_system(n, exps, seed, rng)
    verify(eq1, eq2, ctx)
    s1, s2 = eq1.stats(), eq2.stats()
    F = ctx["F"]
    print(f"n={n} exps={exps} seed={seed}: H/G identities OK on 100 random points each; "
          f"true key (t0,t1)=({ctx['ts'][0]},{ctx['ts'][1]}) is a root for "
          f"(eps0,eps1,eps2)=({ctx['eps'][0]},{ctx['eps'][1]},{ctx['eps2']})")
    for nm, s in (("Eq1", s1), ("Eq2", s2)):
        hist = " ".join(f"{d}:{c}" for d, c in sorted(s["bool_deg_hist"].items()))
        print(f"   {nm}: {s['terms']:6d} terms, deg_t0={s['deg_t0']:6d} deg_t1={s['deg_t1']:6d} "
              f"(2^n-1={F.q-1}), max Boolean degree {s['bool_deg_max']}; Boolean-degree histogram {hist}")
    bez = s1["deg_t0"] * s1["deg_t1"] if s1["deg_t1"] else 0
    bez2 = s2["deg_t0"] * s2["deg_t1"] if s2["deg_t1"] else 0
    print(f"   Bezout bound deg_t0*deg_t1: Eq1 {bez} (2^2n={4**n}), Eq2 {bez2}; "
          f"product of total degrees = {(s1['deg_t0']+s1['deg_t1'])*(s2['deg_t0']+s2['deg_t1'])}")
    nroots, spurious = "", ""
    if do_roots:
        roots = count_roots(eq1, eq2, ctx)
        nroots = len(roots)
        spurious = nroots - 1
        assert tuple(ctx["ts"]) in roots, "true key missing from the enumerated root set"
        print(f"   roots of the sparse system in (F*)^2 for this branch: {nroots} "
              f"(true key + {spurious} spurious)")
    print(f"   {time.time()-t0:.1f}s")
    rec = dict(n=n, exps=";".join(map(str, exps)), seed=seed,
               eps0=ctx["eps"][0], eps1=ctx["eps"][1], eps2=ctx["eps2"],
               eq1_terms=s1["terms"], eq1_deg_t0=s1["deg_t0"], eq1_deg_t1=s1["deg_t1"],
               eq1_bool_deg_max=s1["bool_deg_max"],
               eq1_bool_deg_hist=";".join(f"{d}:{c}" for d, c in sorted(s1["bool_deg_hist"].items())),
               eq2_terms=s2["terms"], eq2_deg_t0=s2["deg_t0"], eq2_deg_t1=s2["deg_t1"],
               eq2_bool_deg_max=s2["bool_deg_max"],
               eq2_bool_deg_hist=";".join(f"{d}:{c}" for d, c in sorted(s2["bool_deg_hist"].items())),
               bezout_eq1=bez, bezout_eq2=bez2, roots=nroots, spurious=spurious)
    record("E05_reduced_system_stats.csv", rec, seconds=time.time() - t0)
    return rec


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    do_roots = "--no-roots" not in sys.argv
    if "--sweep" in sys.argv:
        ns = [int(x) for x in args[0].split(",")]
        seed = int(args[1]) if len(args) > 1 else 1
    else:
        ns = [int(args[0])]
        seed = int(args[1]) if len(args) > 1 else 1
    print(f"# params: n in {ns} seed={seed} roots={do_roots}")
    for n in ns:
        run(n, seed, do_roots)
