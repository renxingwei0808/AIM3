#!/usr/bin/env python3
"""E07 / track T4 -- eliminating pt between two first-layer AIM3 relations in the Ore ring
R = F_{2^n}<phi>, and measuring the degree of what is left.

SET-UP.  With x_i = A_i(pt) + b_i and t_i = y_i = S[e_i](x_i), the S-box relation x_i y_i =
x_i^{2^{e_i}} + 1 becomes, in R,

    Lambda_i(t_i) := (phi^{e_i} + t_i) * A_i,        Eq_i:  Lambda_i(t_i)(pt) = c_i(t_i),
    c_i(t_i) = 1 + b_i^{2^{e_i}} + t_i * b_i.

Eq_i is AFFINE in pt (not linear): that is the whole point, and it is why the naive "common root"
formulation has to be corrected -- see MEASUREMENT 1 below.

Eliminating pt in R: the LCLM cofactors L = u*Lambda_0 = v*Lambda_1 give the scalar condition
    E(t_0, t_1) := u(c_0(t_0)) + v(c_1(t_1)) = 0,
which no longer mentions pt.  T4.1 asks how the degree of E in t_0 grows with n; T4.3 makes that
degree the switch that decides whether the attack line is open.  This script measures it EXACTLY
(no extrapolation, no interpolation guesswork) in two independent ways:

  * the Boolean (ANF) degree of E as a function of the n bits of t_0, by a Moebius transform over
    ALL 2^n values of t_0 -- exact by exhaustion;
  * the univariate degree of E as a polynomial over F_{2^n}, from the interpolation coefficients
    c_j = sum_{a != 0} E(a) a^{-j} -- also exact, since every function F_{2^n} -> F_{2^n} is a
    polynomial of degree <= 2^n - 1.

The two are linked by a proved fact (Boolean ANF degree = max weight of a contributing exponent),
so each checks the other.

MEASUREMENTS
  1 (--probe)   numeric probe at the true (t_0,t_1) versus random ones: deg gcrd, deg lclm in R,
                dim(ker Lambda_0 cap ker Lambda_1) and -- the corrected statistic -- whether the
                AFFINE system {Lambda_0(v)=c_0, Lambda_1(v)=c_1} is solvable.
  2 (--degree)  the core measurement: t_1 fixed at its true value, t_0 ranging over all of F_{2^n};
                exact Boolean and univariate degree of E and of the solvability predicate C.
  3 (--joint)   Boolean degree of C in the 2n bits of (t_0,t_1) jointly (small n only).
  4 (--export)  dump the instance (A_0, A_1, b, true key) as JSON for the Sage cross-check.

Usage:
    python3 skew_elimination_probe.py --probe --n 6 7 8 9 10 11 12 13 14
    python3 skew_elimination_probe.py --degree --n 6 7 8 9 10 11 12 --timeout 900
    python3 skew_elimination_probe.py --joint --n 6 7 8 --timeout 600
    python3 skew_elimination_probe.py --export instance.json --n 10
"""
import argparse, json, os, random, sys, time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from _common import record, ROOT                                          # noqa: E402
from aim3.gf2n import GF2n                                               # noqa: E402
from aim3.ore import SkewPoly, from_bit_matrix                           # noqa: E402
from aim3.structures import degenerate_class                             # noqa: E402
from aim3_lab import rand_unit_lu                                        # noqa: E402

CSV = "E07_theory_skew_degrees.csv"
FIELDS = ["n", "e0", "e1", "measurement", "trials", "value", "unit", "exact", "detail"]


class Budget(Exception):
    pass


class Clock:
    """`timeout` in the sense of CLAUDE.md rule 4: a run that hits it is reported as a LOWER BOUND."""
    def __init__(self, seconds):
        self.t0, self.limit = time.time(), seconds

    def check(self):
        if self.limit and time.time() - self.t0 > self.limit:
            raise Budget(f"budget {self.limit}s exhausted")

    @property
    def elapsed(self):
        return time.time() - self.t0


def admissible(n):
    """exponents in the GENERIC class only (COMPLEXITY_ACCOUNTING rule 5)."""
    from math import gcd
    return [e for e in range(2, n) if gcd(e, n) == 1 and degenerate_class(n, e) == []]


# --------------------------------------------------------------------- the instance
class Instance:
    """one random first-layer pair (A_0,b_0,e_0), (A_1,b_1,e_1) with a planted key pt."""

    def __init__(self, n, e0, e1, rng):
        self.n, self.e = n, (e0, e1)
        self.F = F = GF2n(n)
        while True:
            self.mats = [rand_unit_lu(n, rng), rand_unit_lu(n, rng)]
            self.A = [from_bit_matrix(F, m) for m in self.mats]
            self.b = [rng.randrange(F.q), rng.randrange(F.q)]
            self.pt = rng.randrange(1, F.q)
            self.x = [self.A[i].apply(self.pt) ^ self.b[i] for i in (0, 1)]
            if all(self.x):                              # NoZeroInpSB
                break
        self.t = [F.sbox(self.x[i], self.e[i]) for i in (0, 1)]
        # image of the standard basis under A_i, so that Lambda_i(t) is a scalar away from it
        self.g = [[self.A[i].apply(1 << k) for k in range(n)] for i in (0, 1)]
        for i in (0, 1):                                 # the relation really does hold
            assert self.lam(i, self.t[i]).apply(self.pt) == self.rhs(i, self.t[i])

    def lam(self, i, t):
        """Lambda_i(t) = (phi^{e_i} + t) * A_i in R."""
        F, e = self.F, self.e[i]
        return (SkewPoly.monomial(F, 1, e) + SkewPoly(F, [t])) * self.A[i]

    def rhs(self, i, t):
        """c_i(t) = 1 + b_i^{2^{e_i}} + t*b_i."""
        F = self.F
        return 1 ^ F.frob(self.b[i], self.e[i]) ^ F.mul(t, self.b[i])

    def columns(self, i, t):
        """the n images Lambda_i(t)(2^k), computed from the cached A_i images (n mults, not a
        skew-polynomial multiplication)."""
        F, e = self.F, self.e[i]
        return [F.frob(g, e) ^ F.mul(t, g) for g in self.g[i]]


# --------------------------------------------------------------------- F_2 linear algebra
def solve_affine(n, cols, rhs):
    """(solvable, dim_kernel) for sum_k v_k cols[k] = rhs over F_2 (cols are bitmask vectors)."""
    piv, rank = {}, 0
    for c in cols:
        v = c
        while v:
            b = v.bit_length() - 1
            if b in piv:
                v ^= piv[b]
            else:
                piv[b] = v
                rank += 1
                break
    r = rhs
    while r:
        b = r.bit_length() - 1
        if b not in piv:
            return False, n - rank
        r ^= piv[b]
    return True, n - rank


def kernel_dim(n, cols):
    return solve_affine(n, cols, 0)[1]


def preimage(n, cols, rhs):
    """some v with sum_k v_k cols[k] = rhs, or None."""
    piv = {}
    for k, c in enumerate(cols):
        v, p = c, 1 << k
        while v:
            b = v.bit_length() - 1
            if b in piv:
                v ^= piv[b][0]
                p ^= piv[b][1]
            else:
                piv[b] = (v, p)
                break
    r, sol = rhs, 0
    while r:
        b = r.bit_length() - 1
        if b not in piv:
            return None
        r ^= piv[b][0]
        sol ^= piv[b][1]
    return sol


# --------------------------------------------------------------------- exact degree tools
def anf_degree(values, m):
    """EXACT Boolean (ANF) degree of f : {0,1}^m -> Z_2^k given by `values` (a list of 2^m ints,
    each a bit-packed output).  Moebius transform in place, then the largest popcount of a mask
    with a nonzero coefficient."""
    a = list(values)
    step = 1
    while step < (1 << m):
        for i in range(1 << m):
            if i & step:
                a[i] ^= a[i ^ step]
        step <<= 1
    d = -1
    for i, c in enumerate(a):
        if c:
            w = bin(i).count("1")
            if w > d:
                d = w
    return d


def univariate_degree(F, values):
    """EXACT degree of the unique polynomial of degree <= q-1 over F_{2^n} interpolating
    `values` (indexed by the field element).  c_j = sum_{a!=0} f(a) a^{-j} for 1 <= j <= q-1,
    plus the f(0) correction at j = q-1 and j = 0."""
    q = F.q
    f0 = values[0]
    for j in range(q - 1, 0, -1):
        s = 0
        for a in range(1, q):
            if values[a]:
                s ^= F.mul(values[a], F.pow(a, -j))
        if j == q - 1:
            s ^= f0
        if s:
            return j
    return 0 if f0 else -1


# --------------------------------------------------------------------- measurement 1
def measurement_probe(n, e0, e1, trials, rng, clock):
    """deg gcrd / deg lclm in R, kernel intersection and affine solvability, at the true (t_0,t_1)
    and at random ones."""
    inst = Instance(n, e0, e1, rng)
    F = inst.F
    out = {"true": {}, "rand": {}}

    def one(t0, t1):
        L0, L1 = inst.lam(0, t0), inst.lam(1, t1)
        g = L0.gcrd(L1)
        L = L0.lclm(L1)
        c0, c1 = inst.columns(0, t0), inst.columns(1, t1)
        cols = [c0[k] | (c1[k] << n) for k in range(n)]
        rhs = inst.rhs(0, t0) | (inst.rhs(1, t1) << n)
        ok, _ = solve_affine(n, cols, rhs)
        return dict(deg_gcrd=g.degree(), deg_lclm=L.degree(),
                    kerint=kernel_dim(n, cols), solvable=int(ok))

    out["true"] = one(inst.t[0], inst.t[1])
    acc = {}
    for _ in range(trials):
        clock.check()
        r = one(rng.randrange(F.q), rng.randrange(F.q))
        for k, v in r.items():
            acc.setdefault(k, []).append(v)
    out["rand"] = {k: (min(v), max(v), sum(v)) for k, v in acc.items()}
    out["deg_lambda"] = [inst.lam(i, inst.t[i]).degree() for i in (0, 1)]
    return inst, out


# --------------------------------------------------------------------- measurement 2
def selftest():
    """`anf_degree` and `univariate_degree` on functions whose degree is known in closed form.
    Both are the instruments of the core measurement, so they are checked before it is trusted."""
    for m in (4, 6, 8):
        v = [0] * (1 << m); v[5] = 1
        assert anf_degree(v, m) == m, "indicator of one point must have degree m"
        v = [0] * (1 << m); v[5] = v[6] = 1
        assert anf_degree(v, m) == m - 1, "indicator of a 2-set must have degree m-1"
        assert anf_degree([bin(i).count("1") & 1 for i in range(1 << m)], m) == 1
        assert anf_degree([1] * (1 << m), m) == 0
    F = GF2n(6)
    for k in (1, 2, 5, 17, F.q - 2):
        vals = [F.pow(a, k) if a else 0 for a in range(F.q)]
        assert univariate_degree(F, vals) == k, f"univariate_degree wrong on t^{k}"
    vals = [1] * F.q
    assert univariate_degree(F, vals) == 0
    print("selftest OK: anf_degree and univariate_degree agree with the closed-form degrees")


def measurement_degree(n, e0, e1, rng, clock, inst=None, t1_override=None):
    """t_1 = true, t_0 over ALL of F_{2^n}: the eliminated scalar condition E(t_0) and the
    linear-algebra solvability predicate C(t_0), both measured exactly."""
    inst = inst or Instance(n, e0, e1, rng)
    F, q = inst.F, inst.F.q
    t1 = inst.t[1] if t1_override is None else t1_override
    L1 = inst.lam(1, t1)
    c1 = inst.rhs(1, t1)
    cols1 = inst.columns(1, t1)
    rhs1 = c1 << n

    Evals, Cvals = [0] * q, [0] * q
    degenerate = 0
    for t0 in range(q):
        clock.check()
        # (a) the eliminated condition, from the LCLM cofactors in R
        L0 = inst.lam(0, t0)
        L, u, v = L0.lclm_cofactors(L1)
        Evals[t0] = u.apply(inst.rhs(0, t0)) ^ v.apply(c1)
        if L.degree() != L0.degree() + L1.degree():
            degenerate += 1               # gcrd(Lambda_0, Lambda_1) nontrivial in R
        # (b) the affine system, by F_2 linear algebra on the action
        c0 = inst.columns(0, t0)
        cols = [c0[k] | (cols1[k] << n) for k in range(n)]
        Cvals[t0] = int(solve_affine(n, cols, inst.rhs(0, t0) | rhs1)[0])

    zeroE = [t for t in range(q) if Evals[t] == 0]
    solvC = [t for t in range(q) if Cvals[t]]

    # Structural prediction for solvC.  Eq_1 alone (t_1 true) has solution set {pt, pt + w} with
    # A_1(w) = a_1 = t_1^{r_1} the kernel variable of the second S-box; each of the two candidate
    # pt values forces one value of t_0.  So |C| should be exactly 2 and its second element should
    # be S[e_0](A_0(pt+w) + b_0).
    r1 = F.r_inv(inst.e[1])
    a1 = F.pow(t1, r1)
    w = preimage(n, [inst.A[1].apply(1 << k) for k in range(n)], a1)
    pt2 = inst.pt ^ w if w is not None else None
    x2 = inst.A[0].apply(pt2) ^ inst.b[0] if pt2 is not None else 0
    t0_second = F.sbox(x2, inst.e[0]) if x2 else None
    predicted = sorted({inst.t[0]} | ({t0_second} if t0_second is not None else set())) \
        if t1_override is None else None
    return inst, dict(
        bool_deg_E=anf_degree(Evals, n),
        bool_deg_C=anf_degree(Cvals, n),
        univ_deg_E=univariate_degree(F, Evals),
        n_zero_E=len(zeroE), n_solv_C=len(solvC),
        true_in_zeroE=int(inst.t[0] in zeroE), true_in_solvC=int(inst.t[0] in solvC),
        C_subset_ZE=int(all(Evals[t] == 0 for t in solvC)),
        ZE_subset_C=int(all(Cvals[t] for t in zeroE)),
        degenerate_t0=degenerate,
        second_solution_explained=(int(sorted(solvC) == predicted) if predicted is not None else ""),
        predicted_C=predicted,
        deg_lclm=inst.lam(0, inst.t[0]).degree() + L1.degree(),
    )


# --------------------------------------------------------------------- measurement 3
def measurement_joint(n, e0, e1, rng, clock, inst=None):
    """Boolean degree of the solvability predicate C(t_0,t_1) in all 2n bits."""
    inst = inst or Instance(n, e0, e1, rng)
    q = inst.F.q
    cols0 = [inst.columns(0, t0) for t0 in range(q)]
    rhs0 = [inst.rhs(0, t0) for t0 in range(q)]
    vals = [0] * (q * q)
    solv = 0
    for t1 in range(q):
        clock.check()
        c1, r1 = inst.columns(1, t1), inst.rhs(1, t1) << n
        for t0 in range(q):
            c0 = cols0[t0]
            cols = [c0[k] | (c1[k] << n) for k in range(n)]
            s = int(solve_affine(n, cols, rhs0[t0] | r1)[0])
            vals[t1 * q + t0] = s
            solv += s
    return inst, dict(bool_deg_C_joint=anf_degree(vals, 2 * n), n_solv_joint=solv,
                      true_solv=int(vals[inst.t[1] * q + inst.t[0]]))


# --------------------------------------------------------------------- JSON export
def export_instance(inst, path):
    d = dict(n=inst.n, modulus=inst.F.mod, e=list(inst.e),
             A_bitmatrix=[[list(map(int, row)) for row in m] for m in inst.mats],
             A_skew_coeffs=[[int(c) for c in inst.A[i].c] for i in (0, 1)],
             b=[int(x) for x in inst.b], pt=int(inst.pt),
             x=[int(x) for x in inst.x], t=[int(t) for t in inst.t],
             c=[int(inst.rhs(i, inst.t[i])) for i in (0, 1)],
             convention=("A_skew_coeffs[i] are the coefficients of sum_j a_j phi^j, i.e. "
                         "A_i(v) = sum_j a_j v^{2^j}; Lambda_i = (phi^{e_i} + t_i)*A_i; "
                         "Eq_i: Lambda_i(pt) = c_i.  Bit i of a field element is the coefficient "
                         "of X^i modulo `modulus` (AIMer v3.0 section 4.1.1)."))
    with open(path, "w") as f:
        json.dump(d, f, indent=1)
    return path


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=[6, 7, 8, 9, 10])
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--degree", action="store_true")
    ap.add_argument("--joint", action="store_true")
    ap.add_argument("--control", action="store_true",
                    help="repeat the core measurement with a RANDOM (wrong) t_1")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--export", default=None, help="write the instance of the FIRST n as JSON")
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--timeout", type=float, default=900.0, help="per (n, measurement) budget")
    ap.add_argument("--seed", type=int, default=20260901)
    ap.add_argument("--no-csv", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    print(f"# params: n={a.n} probe={a.probe} degree={a.degree} joint={a.joint} "
          f"trials={a.trials} timeout={a.timeout} seed={a.seed}")

    def row(n, e0, e1, meas, value, unit, exact, detail, seconds, censored=0):
        print(f"n={n} e=({e0},{e1}) {meas:<18} {value:<8} {unit:<12} {exact:<11} {detail}")
        if not a.no_csv:
            record(CSV, dict(n=n, e0=e0, e1=e1, measurement=meas, trials=a.trials,
                             value=value, unit=unit, exact=exact, detail=detail),
                   seconds=seconds, censored=censored)

    for n in a.n:
        exps = admissible(n)
        if not exps:
            print(f"n={n}: no generic exponent -- skipped")
            continue
        # COMPLEXITY_ACCOUNTING rule 5: only exponents of the GENERIC class.  For n <= 10 and
        # n = 12 the only generic exponent is n-1, so both S-boxes then get the same e; that is
        # not an AIM3 parameter shape (the spec uses distinct e_i) but it keeps the instance out
        # of the degenerate classes, which is what matters for a degree measurement.
        e0 = exps[0]
        e1 = exps[1] if len(exps) > 1 else exps[0]
        rng = random.Random(a.seed + n)

        if a.probe:
            clock = Clock(a.timeout)
            try:
                inst, r = measurement_probe(n, e0, e1, a.trials, rng, clock)
                t, rd = r["true"], r["rand"]
                row(n, e0, e1, "deg_lclm_true", t["deg_lclm"], "phi-degree", "exact",
                    f"deg Lambda = {r['deg_lambda']}; deg gcrd = {t['deg_gcrd']}", clock.elapsed)
                row(n, e0, e1, "deg_lclm_rand", f"{rd['deg_lclm'][0]}..{rd['deg_lclm'][1]}",
                    "phi-degree", "exact",
                    f"deg gcrd in {rd['deg_gcrd'][0]}..{rd['deg_gcrd'][1]} over {a.trials} random (t0,t1)",
                    clock.elapsed)
                row(n, e0, e1, "kerint_true", t["kerint"], "F2-dim", "exact",
                    f"random: {rd['kerint'][0]}..{rd['kerint'][1]}", clock.elapsed)
                row(n, e0, e1, "affine_solvable", t["solvable"], "0/1", "exact",
                    f"random (t0,t1): {rd['solvable'][2]}/{a.trials} solvable", clock.elapsed)
            except Budget as ex:
                row(n, e0, e1, "probe", "", "", "lower_bound", str(ex), a.timeout, censored=1)

        if a.degree:
            clock = Clock(a.timeout)
            try:
                inst, r = measurement_degree(n, e0, e1, rng, clock)
                d = (f"|Z(E)|={r['n_zero_E']} |C|={r['n_solv_C']} "
                     f"C = predicted 2-set: {r['second_solution_explained']} "
                     f"true in Z(E)/C = {r['true_in_zeroE']}/{r['true_in_solvC']} "
                     f"C<=Z(E)={r['C_subset_ZE']} Z(E)<=C={r['ZE_subset_C']} "
                     f"degenerate t0 = {r['degenerate_t0']}")
                row(n, e0, e1, "bool_deg_E", r["bool_deg_E"], f"of n={n} bits", "exact", d, clock.elapsed)
                row(n, e0, e1, "univ_deg_E", r["univ_deg_E"], f"of q-1={inst.F.q-1}", "exact",
                    f"deg lclm = {r['deg_lclm']}", clock.elapsed)
                row(n, e0, e1, "bool_deg_C", r["bool_deg_C"], f"of n={n} bits", "exact", d, clock.elapsed)
                if a.export and n == a.n[-1]:
                    print("   exported:", export_instance(inst, a.export))
            except Budget as ex:
                row(n, e0, e1, "bool_deg_E", "", f"of n={n} bits", "lower_bound", str(ex),
                    a.timeout, censored=1)

        if a.control:
            clock = Clock(a.timeout)
            F = GF2n(n)
            try:
                inst0 = Instance(n, e0, e1, random.Random(a.seed + n))
                bad = random.Random(a.seed + 7777 + n).randrange(1, F.q)
                while bad == inst0.t[1]:
                    bad = random.Random(a.seed + 7778 + n).randrange(1, F.q)
                _, r = measurement_degree(n, e0, e1, rng, clock, inst=inst0, t1_override=bad)
                row(n, e0, e1, "control_bool_deg_C", r["bool_deg_C"], f"of n={n} bits", "exact",
                    f"WRONG t_1: |C|={r['n_solv_C']} |Z(E)|={r['n_zero_E']} "
                    f"true t_0 in C = {r['true_in_solvC']}", clock.elapsed)
            except Budget as ex:
                row(n, e0, e1, "control_bool_deg_C", "", f"of n={n} bits", "lower_bound",
                    str(ex), a.timeout, censored=1)

        if a.joint:
            clock = Clock(a.timeout)
            try:
                inst, r = measurement_joint(n, e0, e1, rng, clock)
                row(n, e0, e1, "bool_deg_C_joint", r["bool_deg_C_joint"], f"of 2n={2*n} bits",
                    "exact", f"|C|={r['n_solv_joint']} of {inst.F.q**2}; true pair solvable="
                             f"{r['true_solv']}", clock.elapsed)
            except Budget as ex:
                row(n, e0, e1, "bool_deg_C_joint", "", f"of 2n={2*n} bits", "lower_bound",
                    str(ex), a.timeout, censored=1)

        if a.export and not a.degree:
            inst = Instance(n, e0, e1, random.Random(a.seed + n))
            print("   exported:", export_instance(inst, a.export))
            break


if __name__ == "__main__":
    main()
