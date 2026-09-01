#!/usr/bin/env python3
"""Single enhanced-Mersenne S-box: (i) dimension of Boolean relations of degree<=3 vanishing on the graph
V_e = {(x,S_e(x)) : x != 0}; (ii) degree-3 Macaulay rank / syzygies of {P1,P2,P3}; (iii) degree-4 Macaulay
degree falls to degree 3.  Establishes: first-fall degree 4 with ~2n new cubics per S-box (generic e).

Usage: python3 sbox_degree_fall.py n e [--no-mac4] [--no-points] [--censored] [--backend python|sage]
  --no-mac4   skip the (expensive) D=4 pass
  --no-points skip dim I(V)_<=3 (the O(2^n) point matrix); the D=4 fall count is unaffected
  --pointset-deg D
              measure dim I(V_e)_<=D instead of <=3 and write the row to
              results/E03_server_pointset.csv (the D=3 table keeps its own header).  Only
              meaningful when 2^n - 1 >> C(2n,<=D): at D=4 that needs n >= 16 (at n=14 the graph
              has 16383 points against C(28,<=4) = 24158 monomials, so ANY instance has >= 7775
              relations from the point deficit alone).
  --censored  mark the CSV row as a censored run (its D=4 numbers are then missing, never a value)
  --backend   'python' (default, src/aim3/boolmac.py) or 'sage' (M4RI via
              experiments/E04_full_system_hilbert/sage/rank_backend.sage; scratch dir $AIM3_SCRATCH).
              Both engines are verified to agree exactly by server/jobs/verify_backend.py.
Valid only when 2^n - 1 >= C(2n,<=3) for (i) to be meaningful as an ideal statement (n>=14), but the
Macaulay counts (ii),(iii) are meaningful for all n. Pure Python: n<=13 for Mac4 in minutes; use magma/ for larger.
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS, record
import aim3.legacy_fields                     # noqa: F401  (adds n=17..24 to legacy CONWAY)
from aim3_lab import F, sbox, rank_gf2
from aim3_sbox_syz import sbox_system
from aim3.boolmac import macaulay_stats
from aim3.boolmac_export import macaulay_stats_backend, relation_dim_sage
from aim3.structures import degenerate_class

BACKEND = "python"
SCRATCH = os.environ.get("AIM3_SCRATCH")

def relation_dim(points, nv, deg, backend=None):
    """(dim of the relations, #monomials, rank).  `backend='sage'` uses the M4RI engine."""
    if (backend or BACKEND) == "sage":
        return relation_dim_sage(points, nv, deg, workdir=SCRATCH, tag=f"pts_nv{nv}")
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

def macaulay(eqs, nv, D, backend=None):
    """(rows, rank, dim(span ∩ deg<=D-1)) — src/aim3/boolmac.py, identical to legacy/aim3_macaulay.macaulay
    (checked by tests/test_core.py::test_boolmac_matches_legacy).  `backend='sage'` routes the same
    matrix through M4RI (server/jobs/verify_backend.py checks the two agree exactly)."""
    st = macaulay_stats_backend(eqs, nv, D, backend=backend or BACKEND, workdir=SCRATCH,
                                tag=f"e03_nv{nv}")
    return st.rows, st.rank, st.dim_low


if __name__ == "__main__":
    n, e = int(sys.argv[1]), int(sys.argv[2])
    PDEG = int(sys.argv[sys.argv.index("--pointset-deg") + 1]) if "--pointset-deg" in sys.argv else 3
    if "--backend" in sys.argv:
        BACKEND = sys.argv[sys.argv.index("--backend") + 1]
        if BACKEND not in ("python", "sage"):
            raise SystemExit(f"--backend must be python or sage, not {BACKEND!r}")
    print(f"# params: n={n} e={e} backend={BACKEND} pointset_deg={PDEG} "
          f"flags={[a for a in sys.argv[3:]]}")
    fld, P1, P2, P3 = sbox_system(n, e)
    allq = P1 + P2 + P3
    nv = 2 * n
    t0 = time.time()
    if PDEG != 3:
        # dim I(V_e)_<=PDEG on its own, written to its own table; the Macaulay side of the
        # comparison is D=PDEG+1 (results/E03_server_D5.csv) and is produced by server/jobs/e03_d5.py
        from math import comb as _comb
        pts = [x | (sbox(fld, x, e) << n) for x in range(1, fld.q)]
        cols = sum(_comb(nv, d) for d in range(PDEG + 1))
        meaningful = (fld.q - 1) >= cols
        print(f"   point set: {len(pts)} points vs C({nv},<={PDEG}) = {cols} monomials -- "
              f"{'meaningful' if meaningful else 'DOMINATED BY THE POINT DEFICIT, not an ideal statement'}")
        dimP, totP, rkP = relation_dim(pts, nv, PDEG)
        print(f"n={n} e={e}: dim I(V)_<={PDEG} = {dimP} of {totP} monomials "
              f"(point-matrix rank {rkP}); deficit floor = {max(0, cols - (fld.q - 1))}")
        record("E03_server_pointset.csv",
               dict(n=n, e=e, deg=PDEG, npoints=len(pts), monomials=totP, point_rank=rkP,
                    dim_IV=dimP, deficit_floor=max(0, cols - (fld.q - 1)),
                    meaningful=int(meaningful),
                    degenerate=";".join(degenerate_class(n, e))),
               seconds=time.time() - t0, tier="B")
        raise SystemExit(0)
    if "--no-points" in sys.argv:
        dim3, tot3 = "", ""
    else:
        pts = [x | (sbox(fld, x, e) << n) for x in range(1, fld.q)]
        dim3, tot3, _ = relation_dim(pts, nv, 3)
    r3, rk3, _ = macaulay(allq, nv, 3)
    rec = dict(n=n, e=e, IV_deg3=dim3, mac3_rows=r3, mac3_rank=rk3, syz3=r3 - rk3,
               excess_over_mac3=(dim3 - rk3) if dim3 != "" else "",
               mac4_rows="", mac4_rank="", mac4_dim_deg_le3="", falls_D4="", outside_ideal="",
               degenerate=";".join(degenerate_class(n, e)))
    print(f"n={n} e={e}: dim I(V)_<=3 = {dim3 if dim3 != '' else 'skipped'} (monomials {tot3}); "
          f"Mac3 rank {rk3} rows {r3} syz3 {r3-rk3} (=8n? {r3-rk3==8*n}); "
          f"excess over Mac3 = {dim3-rk3 if dim3 != '' else '-'} (2n={2*n})")
    if "--no-mac4" not in sys.argv:
        r4, rk4, low4 = macaulay(allq, nv, 4)
        rec.update(mac4_rows=r4, mac4_rank=rk4, mac4_dim_deg_le3=low4, falls_D4=low4 - rk3,
                   outside_ideal=(dim3 - low4) if dim3 != "" else "")
        print(f"   Mac4: rows {r4} rank {rk4} dim(span\u2229deg<=3) {low4}; degree falls D4->3 = {low4-rk3}; "
              f"cubics outside <P1,P2,P3> up to D=4: {dim3-low4 if dim3 != '' else 'n/a (--no-points)'}")
    print(f"   {time.time()-t0:.1f}s")
    record("E03_sbox_degree_fall.csv", rec, seconds=time.time() - t0,
           censored=1 if "--censored" in sys.argv else 0)
