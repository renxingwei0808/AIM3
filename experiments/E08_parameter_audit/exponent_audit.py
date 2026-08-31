#!/usr/bin/env python3
"""Exponent audit for all admissible e (gcd(e,n)=1, 2<=e<=n-2) at n in {128,192,256} (or a given n):
   wt(r)=e^{-1} mod n, wt(-r)=n-wt(r), degenerate Frobenius-orbit tags, and a 'goodness' score.
Usage: python3 exponent_audit.py [n]   -> prints table; writes results/E08_audit_<n>.csv
"""
import sys, os, time
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS, record_table
from aim3.structures import weight_data, degenerate_class
REC = {128: (3, 7, 11), 192: (11, 23, 47), 256: (3, 7, 11, 19)}
def audit(n):
    rows = []
    for e in range(2, n - 1):
        if gcd(e, n) != 1:
            continue
        w = weight_data(n, e)
        tags = degenerate_class(n, e)
        score = min(w["wt_r"], w["wt_minus_r"]) / (n / 2)   # 1.0 = balanced, small = suspicious
        rows.append(dict(n=n, e=e, wt_r=w["wt_r"], wt_minus_r=w["wt_minus_r"], balance=round(score, 3),
                         degenerate=";".join(tags), recommended=e in REC.get(n, ())))
    return rows
if __name__ == "__main__":
    ns = [int(sys.argv[1])] if len(sys.argv) > 1 else [128, 192, 256]
    print(f"# params: n in {ns}")
    for n in ns:
        t0 = time.time()
        rows = audit(n)
        record_table(f"E08_audit_{n}.csv", rows, seconds=time.time() - t0)
        worst = sorted(rows, key=lambda r: r["balance"])[:8]
        print(f"n={n}: {len(rows)} admissible e; recommended:", [(r['e'], r['wt_r'], r['wt_minus_r'], r['balance']) for r in rows if r['recommended']])
        print("   least balanced:", [(r['e'], r['wt_r'], r['wt_minus_r']) for r in worst])
        print("   degenerate-tagged:", [(r['e'], r['degenerate']) for r in rows if r['degenerate']][:12], "...")
        print(f"   written results/E08_audit_{n}.csv  ({time.time()-t0:.1f}s)")
