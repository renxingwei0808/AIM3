#!/usr/bin/env python3
"""Reproduce AIMer v3.0 Table 4 exactly, and compute the alternative baselines every new attack must be compared to.

FINDING (verified 2026-08-31): Table 4 = generic semi-regular series (1-z^2)^m/(1-z)^(N-k) with cost
2^k * C(N-k+d, d)^2 (omega=2). The Boolean series (1+z)^N/((1-z)(1+z^2)^m) that the spec writes as eq.(2)
gives *lower* numbers (215.7/308.8/523.6). Both are recorded as baselines B1 (official) and B2 (Boolean).
Usage: python3 table4_reproduce.py [--csv]
"""
import sys, os, csv
from math import comb, log2
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS

PARAMS = [("AIM3-I", 128, 2, 57), ("AIM3-III", 192, 2, 102), ("AIM3-V", 256, 3, 150)]

def dreg_generic(N, m, maxd=140):
    a = [0] * (maxd + 1)
    for k in range(0, min(m, maxd // 2) + 1):
        a[2 * k] += (-1) ** k * comb(m, k)
    res = [0] * (maxd + 1)
    for i in range(maxd + 1):
        ci = comb(N + i - 1, i)
        for j in range(maxd + 1 - i):
            res[i + j] += a[j] * ci
    for d in range(1, maxd + 1):
        if res[d] <= 0:
            return d

def dreg_boolean(N, m, maxd=None):
    maxd = maxd or min(N, 140)
    num = [comb(N, i) for i in range(maxd + 1)]
    cum, s = [], 0
    for i in range(maxd + 1):
        s += num[i]; cum.append(s)
    res = [0] * (maxd + 1)
    for k in range(0, maxd // 2 + 1):
        c = (-1) ** k * comb(m + k - 1, k)
        for i in range(maxd + 1 - 2 * k):
            res[i + 2 * k] += c * cum[i]
    for d in range(1, maxd + 1):
        if res[d] <= 0:
            return d

def cost(N, k, d, mono="N+d"):
    return k + 2 * log2(comb(N - k + d, d) if mono == "N+d" else comb(N - k, d))

def best(N, m, series, mono, kmax):
    out = None
    for k in range(0, kmax):
        d = series(N - k, m)
        if not d:
            continue
        t = cost(N, k, d, mono)
        if out is None or t < out[0]:
            out = (round(t, 1), k, d)
    return out

if __name__ == "__main__":
    rows = []
    for name, n, ell, koff in PARAMS:
        N, m = (ell + 1) * n, 3 * (ell + 1) * n
        d = dreg_generic(N - koff, m)
        official = round(cost(N, koff, d), 1)
        b1 = best(N, m, dreg_generic, "N+d", N - 60)
        b2 = best(N, m, dreg_boolean, "N+d", N - 60)
        p1 = best(N, (ell + 1) * n, dreg_boolean, "N+d", N - 60)
        p13 = best(N, 2 * (ell + 1) * n, dreg_boolean, "N+d", N - 60)
        print(f"{name}: official k={koff} d={d} -> {official}  | B1 generic-series optimum {b1} | "
              f"B2 Boolean-series optimum {b2} | Boolean P1-only {p1} | Boolean P1+P3 {p13}")
        rows.append(dict(scheme=name, n=n, ell=ell, official_k=koff, official_dreg=d, official_log2T=official,
                         B1=b1, B2=b2, P1only=p1, P1P3=p13))
    if "--csv" in sys.argv:
        os.makedirs(RESULTS, exist_ok=True)
        with open(os.path.join(RESULTS, "E02_baselines.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
        print("written results/E02_baselines.csv")
