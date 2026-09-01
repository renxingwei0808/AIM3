# Server batch 2, task 4 -- independent cross-validation and scale-up of the Tier-A T4 measurement.
#
# The Tier-A probe (`../skew_elimination_probe.py`, results/E07_theory_skew_degrees.csv) computes in a
# hand-written Ore ring and reaches n = 14.  This script recomputes the same quantities with Sage's
# `OrePolynomialRing(K, K.frobenius_endomorphism())` -- a completely independent implementation -- and
# then goes past n = 14 in two ways:
#
#   (a) RING LEVEL (cheap at every n): deg Lambda_i, deg gcrd and deg lclm at the true (t0,t1) and at
#       random ones, dim(ker Lambda_0 cap ker Lambda_1), and the affine solvability predicate.
#
#   (b) THE DEGREE IN t0, D(n) -- measured twice, by two instruments that must agree:
#         * SYMBOLIC: build the same Ore ring over L = K(T) with the twist a -> a^2, T -> T^2 (the
#           Frobenius of L; T is the operator "multiply by t0", and phi.(t0) = t0^2.phi forces T -> T^2).
#           Then Lambda_0(T) = (phi^{e0} + T)*A_0 is an exact object and every coefficient of the monic
#           lclm is an exact element of K(T).  D_sym = max degree of a numerator.  Cost is independent
#           of 2^n, so this runs to n = 24.
#         * EVALUATION-INTERPOLATION (what the brief asks for): t0 ranges over ALL of F_{2^n}, the monic
#           lclm is recomputed for each value, and each coefficient is interpolated exactly
#           (c_j = sum_{a != 0} f(a) a^{-j}, scanned from the top; plus the ANF degree by Moebius).
#           Cost 2^n lclm computations, so it stops where the budget stops -- and a stopped run is
#           reported as censored, never as a value (CLAUDE.md rule 4).
#       The two are different quantities and are reported separately: D_sym is the degree of a formal
#       polynomial in the transcendental T, D_interp the degree of the induced FUNCTION on F_{2^n},
#       i.e. D_sym reduced modulo T^{2^n} - T.  They coincide exactly when D_sym < 2^n.
#
#   (c) The formal degree of the eliminated condition E, without expanding it.  With
#       L = u*Lambda_0 = v*Lambda_1 the eliminated scalar condition is E(t0) = u(c_0(t0)) + v(c_1),
#       and since c_0(t0) = beta + b_0*t0 is affine, u(c_0) = sum_j u_j*(beta^{2^j} + b_0^{2^j}*T^{2^j}).
#       So deg_T E = max_j (deg u_j + 2^j) BEFORE the reduction T^{2^n} = T -- computable from the
#       cofactor degrees alone.  This is the number that decides whether a resultant in t0 is
#       affordable, and it is what E05's censored resultants were running into.
#
# Usage (driven by server/jobs/t4_ore_scale.py, which exports the instance):
#   sage ore_scale.sage INSTANCE.json --out OUT.json [--trials 32] [--no-interp] [--no-symbolic]
#                                     [--interp-budget SECONDS]
import json, sys, time

try:
    sys.stdout.reconfigure(line_buffering=True)   # a buffered log is unpollable on a long run
except Exception:
    pass

argv = sys.argv[1:]
path = argv[0]
opts = argv[1:]


def opt(name, default=None):
    if name in opts:
        i = opts.index(name)
        return opts[i + 1] if i + 1 < len(opts) else default
    return default


outpath = opt("--out", "ore_scale_out.json")
trials = int(opt("--trials", "32"))
do_interp = "--no-interp" not in opts
do_sym = "--no-symbolic" not in opts
interp_budget = float(opt("--interp-budget", "3600"))
sym_budget = float(opt("--symbolic-budget", "3600"))

D = json.load(open(path))
n, mod = int(D["n"]), int(D["modulus"])
e0, e1 = int(D["e"][0]), int(D["e"][1])
q = int(2) ** n
print("# params: instance=%s n=%d e=(%d,%d) trials=%d interp=%s symbolic=%s budgets=(%.0f,%.0f)s"
      % (path, n, e0, e1, trials, do_interp, do_sym, interp_budget, sym_budget))

P = PolynomialRing(GF(2), "X")
X = P.gen()
f = P([(mod >> i) & 1 for i in range(n + 1)])
K = GF(2 ** n, "w", modulus=f)
R = OrePolynomialRing(K, K.frobenius_endomorphism(), "phi")
phi = R.gen()


def elt(v):
    v = int(v)
    return K([(v >> i) & 1 for i in range(n)])


def num(a):
    lst = a.polynomial().list()
    return sum(int(c) << i for i, c in enumerate(lst))


A = [R([elt(c) for c in D["A_skew_coeffs"][i]]) for i in (0, 1)]
b = [elt(v) for v in D["b"]]
t = [elt(v) for v in D["t"]]
c = [elt(v) for v in D["c"]]
pt = elt(D["pt"])
ee = (e0, e1)

res = {"n": n, "e0": e0, "e1": e1, "instance": path}

# ------------------------------------------------------------------ (a) the ring level
t_a = time.time()
Lam = [(phi ** ee[i] + t[i]) * A[i] for i in (0, 1)]
res["lambda_applies"] = [bool(Lam[i](pt) == c[i]) for i in (0, 1)]
if not all(res["lambda_applies"]):
    print("!! Lambda_i(pt) != c_i -- the instance or the convention is wrong; refusing to measure")
    json.dump(res, open(outpath, "w"), indent=1)
    sys.exit(int(2))
res["deg_Lambda"] = [int(Lam[i].degree()) for i in (0, 1)]
g_true = Lam[0].right_gcd(Lam[1])
l_true = Lam[0].left_lcm(Lam[1])
print("   [%.1fs] gcrd/lclm at the true key" % (time.time() - t_a))
res["deg_gcrd_true"] = int(g_true.degree())
res["deg_lclm_true"] = int(l_true.degree())
res["deg_sum"] = int(Lam[0].degree() + Lam[1].degree())


def root_space_dim(F):
    """dim_F2 ker F acting on K (the action factors through phi^n = 1, so it is a bit-matrix rank)."""
    M = matrix(GF(2), n, n)
    for k in range(n):
        v = K([1 if j == k else 0 for j in range(n)])
        col = F(v).polynomial().list()
        for i in range(len(col)):
            M[i, k] = col[i]
    return int(n - M.rank())


def ker_intersection_dim():
    M = matrix(GF(2), 2 * n, n)
    for k in range(n):
        v = K([1 if j == k else 0 for j in range(n)])
        for r, val in ((0, Lam[0](v)), (n, Lam[1](v))):
            col = val.polynomial().list()
            for i in range(len(col)):
                M[r + i, k] = col[i]
    return int(n - M.rank())


res["kerint_true"] = ker_intersection_dim()
print("   [%.1fs] kernel dimensions" % (time.time() - t_a))
res["dim_ker_Lambda"] = [root_space_dim(Lam[i]) for i in (0, 1)]


def rhs(i, ti):
    return K(1) + b[i] ** (2 ** ee[i]) + ti * b[i]


def affine_solvable(t0, t1):
    L0 = (phi ** e0 + t0) * A[0]
    L1 = (phi ** e1 + t1) * A[1]
    M = matrix(GF(2), 2 * n, n)
    for k in range(n):
        v = K([1 if j == k else 0 for j in range(n)])
        for r, val in ((0, L0(v)), (n, L1(v))):
            col = val.polynomial().list()
            for i in range(len(col)):
                M[r + i, k] = col[i]
    r0 = rhs(0, t0).polynomial().list() + [0] * n
    r1 = rhs(1, t1).polynomial().list() + [0] * n
    rv = vector(GF(2), [GF(2)(r0[i]) for i in range(n)] + [GF(2)(r1[i]) for i in range(n)])
    try:
        M.solve_right(rv)
        return True
    except ValueError:
        return False


res["affine_solvable_true"] = bool(affine_solvable(t[0], t[1]))
print("   [%.1fs] affine solvability at the true key" % (time.time() - t_a))

import random as _random
rng = _random.Random(int(20260901) + n)
gr, solv = [], 0
for _ in range(trials):
    r0 = elt(rng.randrange(1, q))
    r1 = elt(rng.randrange(1, q))
    L0 = (phi ** e0 + r0) * A[0]
    L1 = (phi ** e1 + r1) * A[1]
    gr.append(int(L0.right_gcd(L1).degree()))
    if affine_solvable(r0, r1):
        solv += 1
res["deg_gcrd_rand_min"] = int(min(gr))
res["deg_gcrd_rand_max"] = int(max(gr))
res["affine_solvable_rand"] = int(solv)
res["trials"] = int(trials)
res["seconds_ring"] = float(time.time() - t_a)
print("deg Lambda = %s ; deg gcrd(true) = %d ; deg lclm(true) = %d (sum %d)"
      % (res["deg_Lambda"], res["deg_gcrd_true"], res["deg_lclm_true"], res["deg_sum"]))
print("dim(ker L0 cap ker L1) = %d ; affine solvable at the true key: %s ; at random: %d/%d"
      % (res["kerint_true"], res["affine_solvable_true"], solv, trials))
print("deg gcrd over %d random (t0,t1): %d..%d   [%.1fs]"
      % (trials, res["deg_gcrd_rand_min"], res["deg_gcrd_rand_max"], res["seconds_ring"]))

# ------------------------------------------------------------------ (b1) symbolic in t0
if do_sym:
    t_s = time.time()
    PT = PolynomialRing(K, "T")
    T = PT.gen()
    LF = PT.fraction_field()
    sig = LF.hom([LF(T) ** 2], base_map=K.frobenius_endomorphism())
    S = OrePolynomialRing(LF, sig, "psi")
    psi = S.gen()
    A0s = S([LF(a) for a in A[0].coefficients(sparse=False)])
    A1s = S([LF(a) for a in A[1].coefficients(sparse=False)])
    Lam0s = (psi ** e0 + LF(T)) * A0s
    Lam1s = (psi ** e1 + LF(t[1])) * A1s
    gs = Lam0s.right_gcd(Lam1s)
    ls, us, vs = Lam0s.left_xlcm(Lam1s)
    ok = bool(us * Lam0s == ls) and bool(vs * Lam1s == ls)

    def degs(F):
        out = []
        for co in F.coefficients(sparse=False):
            r = LF(co)
            out.append((int(r.numerator().degree()), int(r.denominator().degree())))
        return out

    dl, du, dv = degs(ls), degs(us), degs(vs)
    res["symbolic"] = {
        "deg_gcrd_sym": int(gs.degree()), "deg_lclm_sym": int(ls.degree()),
        "deg_u": int(us.degree()), "deg_v": int(vs.degree()), "cofactors_check": ok,
        "D_sym_lclm": int(max(d for d, _ in dl)),
        "D_sym_lclm_den": int(max(d for _, d in dl)),
        "D_sym_u": int(max(d for d, _ in du)),
        "D_sym_v": int(max(d for d, _ in dv)),
        "lclm_num_degrees": [d for d, _ in dl],
        "u_num_degrees": [d for d, _ in du],
        # (c) formal degree of E = u(c_0(T)) + v(c_1), before the reduction T^{2^n} = T
        "log2_deg_E_formal": float(max([RR(d + 2 ** j).log(2) for j, (d, _) in enumerate(du)])),
        "deg_E_formal_terms": [[j, d] for j, (d, _) in enumerate(du)],
        "seconds": float(time.time() - t_s),
    }
    print("symbolic: deg gcrd = %d, deg lclm = %d, cofactor degs (%d,%d), check=%s"
          % (res["symbolic"]["deg_gcrd_sym"], res["symbolic"]["deg_lclm_sym"],
             res["symbolic"]["deg_u"], res["symbolic"]["deg_v"], ok))
    print("symbolic: D_sym(lclm coefficients in t0) = %d (denominator degree %d); D_sym(u) = %d"
          % (res["symbolic"]["D_sym_lclm"], res["symbolic"]["D_sym_lclm_den"],
             res["symbolic"]["D_sym_u"]))
    print("symbolic: log2 of the FORMAL degree of E in t0 = %.1f  (before T^{2^n} = T)   [%.1fs]"
          % (res["symbolic"]["log2_deg_E_formal"], res["symbolic"]["seconds"]))

# ------------------------------------------------------------------ (b2) evaluation-interpolation
if do_interp:
    t_i = time.time()
    Kelts = list(K)                                  # index of an element = its integer encoding
    order = [None] * q
    for a in Kelts:
        order[num(a)] = a
    print("   interpolation: %d values of t0" % q)
    maxdeg = int(res["deg_lclm_true"])
    vals = [[0] * q for _ in range(maxdeg + 1)]      # vals[j][int(t0)] = coefficient j of the lclm
    degdrop, censored = 0, 0
    for i in range(q):
        if time.time() - t_i > interp_budget:
            censored = 1
            print("!! interpolation budget %.0fs exhausted after %d/%d values -- CENSORED"
                  % (interp_budget, i, q))
            break
        a = order[i]
        L0 = (phi ** e0 + a) * A[0]
        lc = L0.left_lcm(Lam[1])
        co = lc.coefficients(sparse=False)
        if i and (i % 4096 == 0):
            print("   ... %d/%d values, %.0fs" % (i, q, time.time() - t_i))
        if len(co) != maxdeg + 1:
            degdrop += 1
        for j in range(min(len(co), maxdeg + 1)):
            vals[j][i] = num(co[j])
    if not censored:
        # discrete-log tables: sums in K are XORs of the integer encodings, so the interpolation
        # coefficient c_j = sum_{a != 0} f(a) a^{-j} is one integer XOR-loop per j.
        gK = K.multiplicative_generator()
        expt = [0] * (q - 1)
        logt = [-1] * q
        acc = K(1)
        for k in range(q - 1):
            ai = num(acc)
            expt[k] = ai
            logt[ai] = k
            acc = acc * gK

        def anf_degree(v):
            arr = list(v)
            step = 1
            while step < q:
                for i2 in range(q):
                    if i2 & step:
                        arr[i2] ^= arr[i2 ^ step]
                step <<= 1
            d = -1
            for i2, cc in enumerate(arr):
                if cc:
                    wt = bin(i2).count("1")
                    if wt > d:
                        d = wt
            return int(d)

        def univ_degree(v, anf):
            """EXACT degree of the interpolating polynomial of degree <= q-1.  The ANF degree bounds
            the weight of any exponent that can occur, so only those j are scanned; the scan is from
            the top and returns at the first nonzero coefficient."""
            f0 = v[0]
            if anf < 0:
                return int(-1)
            nz = [(logt[v[a]], logt[a]) for a in range(1, q) if v[a]]
            qm = q - 1
            for j in range(q - 1, 0, -1):
                if bin(j).count("1") > anf:
                    continue
                s = 0
                for lf, la in nz:
                    s ^= expt[(lf - j * la) % qm]
                if j == q - 1:
                    s ^= f0
                if s:
                    return int(j)
            return int(0) if f0 else int(-1)

        ad = [anf_degree(vals[j]) for j in range(maxdeg + 1)]
        ud = [univ_degree(vals[j], ad[j]) for j in range(maxdeg + 1)]
        res["interp"] = {"D_interp_lclm": int(max(ud)), "anf_deg_lclm": int(max(ad)),
                         "per_coefficient_univariate": ud, "per_coefficient_anf": ad,
                         "degree_drops": int(degdrop), "values": int(q),
                         "seconds": float(time.time() - t_i)}
        print("interpolation: D_interp(lclm coefficients) = %d of a possible %d ; ANF degree = %d of %d"
              % (res["interp"]["D_interp_lclm"], q - 1, res["interp"]["anf_deg_lclm"], n))
        print("interpolation: %d of %d values of t0 gave a lower-degree lclm   [%.1fs]"
              % (degdrop, q, res["interp"]["seconds"]))
    else:
        res["interp"] = {"censored": 1, "values_done": int(i), "values": int(q),
                         "seconds": float(time.time() - t_i)}

json.dump(res, open(outpath, "w"), indent=1)
print("# wrote %s" % outpath)
