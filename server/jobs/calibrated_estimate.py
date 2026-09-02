#!/usr/bin/env python3
"""Batch-2 task 1.5: hybrid-guessing cost with the solving degree corrected by a measured δ.

**This script prints a table.  It draws no security conclusion, and the report does not either** —
that judgement belongs to the paper, and it needs an extrapolation this batch is not entitled to
make (`COMPLEXITY_ACCOUNTING §4`: the measured δ comes from n ≤ 15).

The model, exactly as the batch prompt states it:

    d_reg(k) = dreg_boolean(N − k, m)              the Boolean semi-regular degree of regularity
    d(k, δ)  = d_reg(k) − δ                        δ = the measured shortfall of the real system
    cost(k)  = 2^k · ( Σ_{i ≤ d} C(N − k, i) )^ω   ω ∈ {2, 2.81}
    answer   = min over k

with N = (ℓ+1)n Boolean variables and m = 3n(ℓ+1) quadratics.

**A deliberate difference from B2.**  `E02/table4_reproduce.py` costs a step as
`C(N−k+d, d)^2` — the *commutative* monomial count, which is what reproduces the official Table 4.
Over a Boolean ring the monomials of degree ≤ d number `Σ_{i≤d} C(N−k, i)`, which is what the prompt
asks for and what this script uses.  The two differ, so the δ = 0, ω = 2 column here is **not**
literally B2; it is "B2's degree with the Boolean monomial count", and both are printed side by
side so the gap is visible rather than hidden.

Usage:
  python3 server/jobs/calibrated_estimate.py [--deltas 0,1] [--omegas 2,2.81]
                                             [--csv E04_server_calibrated.csv] [--kmax-frac 0.9]
"""
import argparse
import csv as _csv
import os
import sys
from math import comb, log2

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
for p in (os.path.join(ROOT, "experiments"),
          os.path.join(ROOT, "experiments", "E02_official_baseline")):
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record                                          # noqa: E402
from table4_reproduce import dreg_boolean                           # noqa: E402

# (name, n, ell) -- AIM3-V has four exponents, so ell = 3
PARAMS = [("AIM3-I", 128, 2), ("AIM3-III", 192, 2), ("AIM3-V", 256, 3)]


def log2_boolean_monomials(N, d):
    """log2 of Σ_{i≤d} C(N, i), the number of Boolean monomials of degree ≤ d in N variables."""
    s = 0
    for i in range(min(d, N) + 1):
        s += comb(N, i)
    return log2(s) if s else 0.0


def baselines():
    """B1 / B2 / B3 from the committed Tier-A tables, so the comparison is not re-derived here."""
    out = {}
    p = os.path.join(ROOT, "results", "E02_baselines.csv")
    if os.path.exists(p):
        for r in _csv.DictReader(open(p)):
            out.setdefault(r["scheme"], {}).update(
                B1=float(r["B1_log2T"]), B2=float(r["B2_log2T"]),
                B1_k=int(r["B1_k"]), B1_dreg=int(r["B1_dreg"]),
                B2_k=int(r["B2_k"]), B2_dreg=int(r["B2_dreg"]))
    p = os.path.join(ROOT, "results", "E11_theory_eval_cost.csv")
    if os.path.exists(p):
        for r in _csv.DictReader(open(p)):
            # the karatsuba row is the cheaper (i.e. the stronger) B3
            d = out.setdefault(r["scheme"], {})
            b3 = float(r["log2_B3"])
            if "B3" not in d or b3 < d["B3"]:
                d["B3"] = b3
                d["B3_mode"] = r["mode"]
    return out


def best_k(N, m, delta, omega, kmax):
    out = None
    for k in range(0, kmax + 1):
        dr = dreg_boolean(N - k, m)
        if not dr:
            continue
        d = dr - delta
        if d < 1:
            continue
        t = k + omega * log2_boolean_monomials(N - k, d)
        if out is None or t < out[0]:
            out = (t, k, dr, d)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deltas", default="0,1")
    ap.add_argument("--omegas", default="2,2.81")
    ap.add_argument("--kmax-frac", type=float, default=0.9,
                    help="search k up to this fraction of N (k = N would be brute force)")
    ap.add_argument("--csv", default="E04_server_calibrated.csv")
    a = ap.parse_args()
    deltas = [int(v) for v in a.deltas.split(",")]
    omegas = [float(v) for v in a.omegas.split(",")]
    base = baselines()

    print(f"# params: deltas={deltas} omegas={omegas} kmax_frac={a.kmax_frac}")
    print("# cost = 2^k * (sum_{i<=d} C(N-k,i))^omega,  d = dreg_boolean(N-k, m) - delta,  min over k")
    print("# NOTE: B1/B2 below use the COMMUTATIVE monomial count C(N-k+d,d)^2 of E02, not this one.\n")
    hdr = (f"{'scheme':<9} {'n':>4} {'l':>2} {'N':>5} {'m':>5} {'d':>3} {'w':>5} "
           f"{'k*':>5} {'d_reg':>6} {'d':>4} {'log2 T':>9}   {'B1':>7} {'B2':>7} {'B3':>7}")
    print(hdr)
    print("-" * len(hdr))
    for name, n, ell in PARAMS:
        N, m = (ell + 1) * n, 3 * n * (ell + 1)
        b = base.get(name, {})
        for delta in deltas:
            for omega in omegas:
                r = best_k(N, m, delta, omega, int(N * a.kmax_frac))
                if r is None:
                    continue
                t, k, dr, d = r
                print(f"{name:<9} {n:>4} {ell:>2} {N:>5} {m:>5} {delta:>3} {omega:>5} "
                      f"{k:>5} {dr:>6} {d:>4} {t:>9.1f}   "
                      f"{b.get('B1', float('nan')):>7.1f} {b.get('B2', float('nan')):>7.1f} "
                      f"{b.get('B3', float('nan')):>7.1f}")
                record(a.csv, dict(
                    scheme=name, n=n, ell=ell, nvars=N, neqs=m, delta=delta, omega=omega,
                    best_k=k, dreg_boolean=dr, d_used=d, log2_T=round(t, 2),
                    monomial_count="boolean sum_{i<=d} C(N-k,i)",
                    B1_log2T=b.get("B1", ""), B2_log2T=b.get("B2", ""), B3_log2T=b.get("B3", ""),
                    B3_mode=b.get("B3_mode", ""),
                    beats_B1=int(t < b["B1"]) if "B1" in b else "",
                    beats_B2=int(t < b["B2"]) if "B2" in b else "",
                    beats_B3=int(t < b["B3"]) if "B3" in b else ""),
                    seconds=0, tier="A")
    print("\n# B1/B2 from results/E02_baselines.csv; B3 (cheapest mode) from "
          "results/E11_theory_eval_cost.csv.")
    print("# No security conclusion is drawn here: delta is measured at n <= 15 and "
          "COMPLEXITY_ACCOUNTING section 4 does not license extrapolating it to n = 128..256.")


if __name__ == "__main__":
    main()
