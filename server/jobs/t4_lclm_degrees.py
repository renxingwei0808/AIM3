#!/usr/bin/env python3
"""Batch-2 task 4, the scale-up half: D(n), the degree in t₀ of the eliminated (lclm) coefficients.

`server/jobs/t4_ore_scale.py` runs Sage's `OrePolynomialRing` and is the *independent* instrument --
that is what cross-validates the Tier-A branch.  It is however about **270× slower per lclm** than the
hand-written `src/aim3/ore.py` (measured: ~1 s versus ~3.7 ms at n = 8..14), which puts every
2^n-evaluation measurement out of its reach above n ≈ 12.  So the division of labour is:

    Sage            -- cross-validation at every n it can reach, ring level to n = 24;
    src/aim3/ore.py -- the 2^n evaluation-interpolation, pushed as far as the budget goes.

Both are run on the SAME instances (built by the Tier-A probe, same seed), and every n where both
produce a number is compared; a disagreement would be a real finding and is reported as one.

WHAT IS MEASURED.  t₁ is fixed at its true value and t₀ ranges over all of F_{2^n}.  For each t₀ the
monic lclm L(t₀) = lclm((φ^{e₀}+t₀)A₀, Λ₁) is formed and its coefficients recorded.  Each coefficient
is then a *function* F_{2^n} → F_{2^n}, and

    D(n)   = max over coefficients of its exact univariate degree  (c_j = Σ_{a≠0} f(a)a^{-j},
             computed with discrete-log tables so one candidate j costs O(2^n) integer XORs),
    A(n)   = max over coefficients of its exact Boolean (ANF) degree (Möbius transform; n ≤ --anf-max
             only, since it costs 2^n·n word operations per coefficient).

D(n) is the quantity that decides the cost of eliminating t₀ by a resultant: a degree-D condition in
t₀ costs at least D field operations to even write down.

Usage:
  python3 server/jobs/t4_lclm_degrees.py --ns 8,9,...,18 [--budget 7200] [--anf-max 16]
                                         [--csv E07_server_lclm_degrees.csv]
"""
import argparse
import os
import sys
import time
from array import array

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP,
          os.path.join(_EXP, "E07_ore_skew_elimination")]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record                                          # noqa: E402
from aim3.ore import SkewPoly                                       # noqa: E402
import skew_elimination_probe as probe                              # noqa: E402


def dlog_tables(F, q):
    """(exp, log) for a multiplicative generator of F_{2^n}; sums in the field are integer XORs."""
    for g in range(2, q):
        seen, x, ok = [0] * q, 1, True
        exp = array("l", [0]) * 0
        exp = array("l", [0] * (q - 1))
        log = array("l", [-1] * q)
        for k in range(q - 1):
            if seen[x]:
                ok = False
                break
            seen[x] = 1
            exp[k] = x
            log[x] = k
            x = F.mul(x, g)
        if ok and x == 1:
            return exp, log
    raise RuntimeError("no generator found")


def anf_degree(vals, q, n):
    a = list(vals)
    step = 1
    while step < q:
        for i in range(q):
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


def univariate_degree(vals, q, exp, log, cap):
    """EXACT degree of the interpolating polynomial, scanning from the top.  Returns
    (degree, exhausted): `exhausted` means the scan hit `cap` candidates without a nonzero
    coefficient, so the answer is a lower bound of −1 rather than a value (never reported as one)."""
    f0 = vals[0]
    nz = [(log[vals[a]], log[a]) for a in range(1, q) if vals[a]]
    if not nz and not f0:
        return -1, False
    qm = q - 1
    tried = 0
    for j in range(q - 1, 0, -1):
        s = 0
        for lf, la in nz:
            s ^= exp[(lf - j * la) % qm]
        if j == q - 1:
            s ^= f0
        if s:
            return j, False
        tried += 1
        if tried >= cap:
            return -1, True
    return (0 if f0 else -1), False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="8,9,10,11,12,13,14,15,16,17,18")
    ap.add_argument("--budget", type=float, default=7200.0, help="per-n wall-clock cap")
    ap.add_argument("--anf-max", type=int, default=16)
    ap.add_argument("--cap", type=int, default=4096, help="candidate exponents scanned from the top")
    ap.add_argument("--seed", type=int, default=20260901)
    ap.add_argument("--csv", default="E07_server_lclm_degrees.csv")
    a = ap.parse_args()
    ns = [int(v) for v in a.ns.split(",")]
    print(f"# params: ns={ns} budget={a.budget}s anf_max={a.anf_max} cap={a.cap} seed={a.seed}")
    print("# instrument: src/aim3/ore.py (the fast one); Sage cross-validates the same instances "
          "in server/jobs/t4_ore_scale.py")

    import random
    for n in ns:
        adm = probe.admissible(n)
        if len(adm) < 1:
            print(f"n={n}: no generic exponent, skipped")
            continue
        rng = random.Random(a.seed + n)
        inst = probe.Instance(n, adm[0], adm[1] if len(adm) > 1 else adm[0], rng)
        F, q = inst.F, 1 << n
        e0, e1 = inst.e
        L1 = (SkewPoly.monomial(F, 1, e1) + SkewPoly(F, [inst.t[1]])) * inst.A[1]
        base0 = SkewPoly.monomial(F, 1, e0)
        Lt = (base0 + SkewPoly(F, [inst.t[0]])) * inst.A[0]
        assert Lt.apply(inst.pt) == inst.rhs(0, inst.t[0]), "instance is inconsistent"
        maxdeg = Lt.lclm(L1).degree()
        print(f"\n== n={n} e=({e0},{e1}) deg Lambda=({Lt.degree()},{L1.degree()}) "
              f"deg lclm={maxdeg}  ({q} values of t0)", flush=True)

        t0 = time.time()
        cols = [array("Q", [0]) * 0 for _ in range(maxdeg + 1)]
        cols = [array("Q", bytes(8 * q)) for _ in range(maxdeg + 1)]
        drops, censored, done = 0, 0, 0
        for v in range(q):
            if (v & 0x3FFF) == 0 and time.time() - t0 > a.budget:
                censored = 1
                break
            lc = (base0 + SkewPoly(F, [v])) * inst.A[0]
            L = lc.lclm(L1).monic()
            cc = L.c
            if len(cc) != maxdeg + 1:
                drops += 1
            for j in range(min(len(cc), maxdeg + 1)):
                cols[j][v] = cc[j]
            done = v + 1
        t_lclm = time.time() - t0
        if censored:
            print(f"   !! budget {a.budget:.0f}s exhausted after {done}/{q} values -- CENSORED, "
                  f"no degree reported ({t_lclm:.0f}s)")
            record(a.csv, dict(n=n, e0=e0, e1=e1, nvalues=q, values_done=done,
                               deg_lclm=maxdeg, D_interp="", anf="", degree_drops=drops,
                               tool="aim3.ore", detail="budget exhausted"),
                   seconds=t_lclm, censored=1, tier="B")
            continue
        print(f"   {q} lclm computations in {t_lclm:.1f}s ({1000*t_lclm/q:.2f} ms each); "
              f"{drops} values of t0 gave a lower-degree lclm", flush=True)

        exp, log = dlog_tables(F, q)
        ud, exhausted = [], 0
        for j in range(maxdeg + 1):
            d, ex = univariate_degree(cols[j], q, exp, log, a.cap)
            ud.append(d)
            exhausted += int(ex)
        ad = []
        if n <= a.anf_max:
            ad = [anf_degree(cols[j], q, n) for j in range(maxdeg + 1)]
        secs = time.time() - t0
        D = max(ud)
        print(f"   D_interp = {D} of a possible {q-1}"
              + (f";  ANF = {max(ad)} of {n}" if ad else ";  ANF not attempted (n > anf-max)")
              + f"   [{secs:.1f}s total]", flush=True)
        print(f"   per-coefficient univariate degrees: {ud}")
        if exhausted:
            print(f"   !! {exhausted} coefficient(s) hit the scan cap -- their degree is unknown, "
                  f"not zero")
        record(a.csv, dict(n=n, e0=e0, e1=e1, nvalues=q, values_done=q, deg_lclm=maxdeg,
                           D_interp=D, D_over_q=round(D / (q - 1), 4),
                           anf=max(ad) if ad else "", degree_drops=drops,
                           per_coefficient=";".join(map(str, ud)),
                           scan_exhausted=exhausted, lclm_seconds=round(t_lclm, 2),
                           ms_per_lclm=round(1000 * t_lclm / q, 3), tool="aim3.ore"),
               seconds=secs, censored=0, tier="B")


if __name__ == "__main__":
    main()
