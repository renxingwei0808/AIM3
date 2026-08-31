#!/usr/bin/env python3
"""Single enhanced-Mersenne S-box: (i) dimension of Boolean relations of degree<=3 vanishing on the graph
V_e = {(x,S_e(x)) : x != 0}; (ii) degree-3 Macaulay rank / syzygies of {P1,P2,P3}; (iii) degree-4 Macaulay
degree falls to degree 3.  Establishes: first-fall degree 4 with ~2n new cubics per S-box (generic e).

Usage: python3 sbox_degree_fall.py n e [--no-mac4] [--censored]
  --no-mac4   skip the (expensive) D=4 pass
  --censored  mark the CSV row as a censored run (its D=4 numbers are then missing, never a value)
Valid only when 2^n - 1 >= C(2n,<=3) for (i) to be meaningful as an ideal statement (n>=14), but the
Macaulay counts (ii),(iii) are meaningful for all n. Pure Python: n<=13 for Mac4 in minutes; use magma/ for larger.
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS, record
from aim3_lab import F, sbox, rank_gf2
from aim3_sbox_syz import sbox_system
from aim3.boolmac import macaulay_stats
from aim3.structures import degenerate_class

def relation_dim(points, nv, deg):
    ms = []
    for d in range(deg + 1):
        ms.extend(itertools.combinations(range(nv), d))
    masks = [sum(1 << i for i in m) for m in ms]
    rows = []
    for p in points:
        v = 0
        for k, mask in enumerate(masks):
            if p & mask == mask:
                v |= 1 << k
        rows.append(v)
    rk, _ = rank_gf2(rows)
    return len(ms) - rk, len(ms), rk

def macaulay(eqs, nv, D):
    """(rows, rank, dim(span ∩ deg<=D-1)) — src/aim3/boolmac.py, identical to legacy/aim3_macaulay.macaulay
    (checked by tests/test_core.py::test_boolmac_matches_legacy)."""
    st = macaulay_stats(eqs, nv, D)
    return st.rows, st.rank, st.dim_low


if __name__ == "__main__":
    n, e = int(sys.argv[1]), int(sys.argv[2])
    print(f"# params: n={n} e={e} flags={[a for a in sys.argv[3:]]}")
    fld, P1, P2, P3 = sbox_system(n, e)
    allq = P1 + P2 + P3
    nv = 2 * n
    t0 = time.time()
    pts = [x | (sbox(fld, x, e) << n) for x in range(1, fld.q)]
    dim3, tot3, _ = relation_dim(pts, nv, 3)
    r3, rk3, _ = macaulay(allq, nv, 3)
    rec = dict(n=n, e=e, IV_deg3=dim3, mac3_rows=r3, mac3_rank=rk3, syz3=r3 - rk3, excess_over_mac3=dim3 - rk3,
               mac4_rows="", mac4_rank="", mac4_dim_deg_le3="", falls_D4="", outside_ideal="",
               degenerate=";".join(degenerate_class(n, e)))
    print(f"n={n} e={e}: dim I(V)_<=3 = {dim3} (monomials {tot3}); Mac3 rank {rk3} rows {r3} syz3 {r3-rk3} "
          f"(=8n? {r3-rk3==8*n}); excess over Mac3 = {dim3-rk3} (2n={2*n})")
    if "--no-mac4" not in sys.argv:
        r4, rk4, low4 = macaulay(allq, nv, 4)
        rec.update(mac4_rows=r4, mac4_rank=rk4, mac4_dim_deg_le3=low4, falls_D4=low4 - rk3, outside_ideal=dim3 - low4)
        print(f"   Mac4: rows {r4} rank {rk4} dim(span\u2229deg<=3) {low4}; degree falls D4->3 = {low4-rk3}; "
              f"cubics outside <P1,P2,P3> up to D=4: {dim3-low4}")
    print(f"   {time.time()-t0:.1f}s")
    record("E03_sbox_degree_fall.csv", rec, seconds=time.time() - t0,
           censored=1 if "--censored" in sys.argv else 0)
