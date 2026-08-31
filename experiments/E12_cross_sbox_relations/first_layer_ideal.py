#!/usr/bin/env python3
"""E12 — the FIRST-LAYER-only AIM3 system {P1,P2,P3} of two S-boxes in the 3n variables (pt,t0,t1),
against the point-set ideal measured by `cross_relations.py`.

Same instance as `cross_relations.py n deg seed`: identical unit-LU matrices and constants, because the
random stream (rand_unit_lu, rand_unit_lu, b0, b1) and the exponent list are reproduced verbatim.

What is measured
  * `mac3_rank`  — rank of the degree-3 Macaulay matrix of the 6n quadratics (3n per S-box) in 3n variables.
  * `IV_deg3`    — dim of ALL Boolean cubics vanishing on the graph {(pt,t0,t1)} (cross_relations.py).
  * `excess = IV_deg3 - mac3_rank` — cubic relations of the joint graph that the degree-3 Macaulay matrix
    of the two S-box ideals does NOT yet produce, i.e. that only become visible at D >= 4.
    The per-S-box first-fall count of E03 (2n each) predicts `excess = 2 * 2n = 4n` if the two S-box
    ideals contribute independently and there is no cross-S-box cubic (E12 verifies the latter directly).
  * `falls_D4 = dim(span_{D=4} n deg<=3) - rank_{D=3}` — the same quantity as E03/E04, to be compared with 4n.
  * PER-S-BOX CONTROL (the point of this script).  Comparing falls_D4 of the joint system with 2n per
    S-box is NOT a valid direct-sum test: the joint Macaulay matrix may also multiply by cross-block
    monomials (t0*t1, pt*t1, ...) that simply do not exist in the 2n-variable single-S-box ring of E03,
    so it can produce extra falls without any cross-S-box relation.  We therefore also run each S-box's
    own 3n quadratics ALONE in the same 3n variables, and report

        cross_excess = falls_D4(both S-boxes) - falls_D4(S-box 0 alone) - falls_D4(S-box 1 alone).

    cross_excess = 0  <=>  the D=4 falls are the direct sum (T1.3);  > 0  <=>  genuine cross-S-box falls.

Usage:
  python3 first_layer_ideal.py n [--D4] [--no-points] [--seed S] [--budget SEC]
                                 [--backend python|sage]
`--backend sage` routes every rank (Macaulay and the point-set matrix) through M4RI; the two
engines are verified to agree exactly by server/jobs/verify_backend.py.
Appends to results/E12_first_layer_ideal.csv.
"""
import sys, os, random, time, argparse
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3_lab import F, sbox, rand_unit_lu, apply_mat
from aim3_macaulay import AffExpr, prod_coords, lin_coords, add_polys
from aim3.boolmac import macaulay_stats, MacaulayTimeout
from aim3.boolmac_export import macaulay_stats_backend
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "E03_syzygy_and_degree_fall"))
from sbox_degree_fall import relation_dim

SCRATCH = os.environ.get("AIM3_SCRATCH")


def build_first_layer(n, exps2, seed=1):
    """The two first-layer S-boxes only: variables pt (bits 0..n-1), t0 (n..2n-1), t1 (2n..3n-1);
    equations P1,P2,P3 per S-box (3n each, 6n in total).  Returns (nv, eqsP1, eqsAll, sol, points).

    The instance is generated with exactly the random stream of cross_relations.py, so the two scripts
    describe the same object for the same (n, seed).
    """
    fld = F(n)
    rng = random.Random(seed)
    A = [rand_unit_lu(n, rng), rand_unit_lu(n, rng)]
    b = [rng.randrange(fld.q), rng.randrange(fld.q)]
    nv = 3 * n

    # a valid key for the sanity check
    while True:
        pt = rng.randrange(1, fld.q)
        xs = [apply_mat(A[i], pt, n) ^ b[i] for i in range(2)]
        if 0 in xs:
            continue
        ts = [sbox(fld, xs[i], exps2[i]) for i in range(2)]
        break
    sol = [0] * nv
    for k in range(n):
        sol[k] = (pt >> k) & 1
        for i in range(2):
            sol[(i + 1) * n + k] = (ts[i] >> k) & 1

    def colimgs(M):
        return [apply_mat(M, 1 << j, n) for j in range(n)]

    eqsP1, eqsAll = [], []
    for i in range(2):
        ci = colimgs(A[i])
        u = AffExpr(fld, {j: ci[j] for j in range(n)}, b[i])                  # x_i = A_i pt + b_i
        t = AffExpr(fld, {(i + 1) * n + k: 1 << k for k in range(n)}, 0)      # t_i
        t2 = AffExpr(fld, {(i + 1) * n + k: fld.pw(1 << k, 2) for k in range(n)}, 0)   # t_i^2
        e = exps2[i]
        P1 = add_polys(prod_coords(u, t), lin_coords(u.frob(e)), lin_coords(AffExpr(fld, {}, 1)))
        P2 = add_polys(prod_coords(u.frob(1), t), prod_coords(u.frob(e), u), lin_coords(u))
        P3 = add_polys(prod_coords(u, t2), prod_coords(u.frob(e), t), lin_coords(t))
        eqsP1 += P1
        eqsAll += P1 + P2 + P3

    def ev(poly):
        s = 0
        for m in poly:
            v = 1
            for j in m:
                v &= sol[j]
            s ^= v
        return s
    assert all(ev(p) == 0 for p in eqsAll), "first-layer system does not vanish at the key"

    pts = []
    for p in range(fld.q):
        x0 = apply_mat(A[0], p, n) ^ b[0]
        x1 = apply_mat(A[1], p, n) ^ b[1]
        if x0 == 0 or x1 == 0:
            continue
        pts.append(p | (sbox(fld, x0, exps2[0]) << n) | (sbox(fld, x1, exps2[1]) << (2 * n)))
    return nv, eqsP1, eqsAll, sol, pts


def measure(name, eqs, nv, D, budget, backend="python"):
    t0 = time.time()
    try:
        return macaulay_stats_backend(eqs, nv, D, backend=backend, budget=budget, workdir=SCRATCH,
                                      tag=f"e12_{name.replace(' ', '')}_nv{nv}"), 0, time.time() - t0
    except MacaulayTimeout as exc:
        print(f"   !! {name} D={D} exceeded {budget:.0f}s -> lower bound only")
        return exc.partial, 1, time.time() - t0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--D4", action="store_true", help="also run the D=4 pass (falls_D4)")
    ap.add_argument("--no-points", action="store_true", help="skip the O(2^n) point-set ideal dimension")
    ap.add_argument("--budget", type=float, default=3600.0)
    ap.add_argument("--backend", choices=["python", "sage"], default="python",
                    help="rank engine: pure Python (default) or Sage/M4RI")
    a = ap.parse_args()
    es = [e for e in range(3, a.n - 1) if gcd(e, a.n) == 1][:2]     # identical to cross_relations.py
    print(f"# params: n={a.n} seed={a.seed} exps={es} D4={a.D4} points={not a.no_points} "
          f"budget={a.budget}s backend={a.backend}")
    t_all = time.time()
    nv, eqsP1, eqsAll, sol, pts = build_first_layer(a.n, es, a.seed)
    n = a.n
    print(f"first-layer system: n={n} exps={es}: {len(eqsAll)} quadratics (6n) in {nv} vars (3n), "
          f"{len(pts)} graph points; vanishes at the key: OK")
    if len(pts) < sum(__import__("math").comb(nv, d) for d in range(4)):
        print(f"   WARNING: only {len(pts)} graph points for "
              f"{sum(__import__('math').comb(nv, d) for d in range(4))} monomials of degree <=3 — the "
              f"point-set ideal dimension is dominated by the point deficit, not by the ideal "
              f"(meaningful only for 2^n >> C(3n,<=3), i.e. n >= 15).")
    st3, cens, sc3 = measure("S_quad(firstlayer)", eqsAll, nv, 3, a.budget, a.backend)
    print(f"   D=3: rows={st3.rows} cols={st3.cols} rank={st3.rank} syz={st3.syz} "
          f"(2*8n={16*n}) dim(span∩deg<=2)={st3.dim_low}")
    rec = dict(n=n, exps=";".join(map(str, es)), seed=a.seed, nvars=nv, neqs=len(eqsAll),
               mac3_rows=st3.rows, mac3_rank=st3.rank, mac3_syz=st3.syz, syz_target=16 * n,
               IV_deg3="", excess_over_mac3="", excess_over_4n="",
               mac4_rows="", mac4_rank="", mac4_dim_le3="", falls_D4="", falls_minus_4n="",
               falls_sbox0="", falls_sbox1="", cross_excess="")
    if not a.no_points:
        t0 = time.time()
        d3, tot3, _ = relation_dim(pts, nv, 3, backend=a.backend)
        exc = d3 - st3.rank
        rec.update(IV_deg3=d3, excess_over_mac3=exc, excess_over_4n=exc - 4 * n)
        print(f"   point-set ideal: dim I(V)_<=3 = {d3} of {tot3} monomials on {len(pts)} points "
              f"({time.time()-t0:.1f}s)")
        print(f"   => excess over the D=3 Macaulay rank = {d3} - {st3.rank} = {exc}   "
              f"(direct-sum prediction 2*2n = {4*n}; difference {exc - 4*n})")
    if a.D4:
        st4, c4, sc4 = measure("S_quad(firstlayer)", eqsAll, nv, 4, a.budget, a.backend)
        cens |= c4
        print(f"   D=4: rows={st4.rows} cols={st4.cols} rank={st4.rank} syz={st4.syz} "
              f"dim(span∩deg<=3)={st4.dim_low}"
              + ("   [CENSORED: lower bound]" if c4 else ""))
        if not c4:
            falls = st4.dim_low - st3.rank
            rec.update(mac4_rows=st4.rows, mac4_rank=st4.rank, mac4_dim_le3=st4.dim_low,
                       falls_D4=falls, falls_minus_4n=falls - 4 * n)
            print(f"   => falls_D4 = {st4.dim_low} - {st3.rank} = {falls}  (4n = {4*n}; "
                  f"difference {falls - 4*n}); low-degree columns C(3n,<=3) = "
                  f"{sum(__import__('math').comb(nv, d) for d in range(4))}")
            # per-S-box control: each S-box alone, in the SAME 3n variables
            per = []
            for i in range(2):
                sub = eqsAll[i * 3 * n:(i + 1) * 3 * n]
                s3, cc3, _ = measure(f"S-box{i}alone", sub, nv, 3, a.budget, a.backend)
                s4, cc4, _ = measure(f"S-box{i}alone", sub, nv, 4, a.budget, a.backend)
                cens |= cc3 | cc4
                f_i = s4.dim_low - s3.rank
                per.append(f_i)
                print(f"   [control] S-box {i} alone in {nv} vars: D3 rank={s3.rank} syz={s3.syz} (8n={8*n}); "
                      f"D4 dim(span∩deg<=3)={s4.dim_low} -> falls_D4 = {f_i} ({f_i/n:.2f} n)")
            cross = falls - sum(per)
            rec.update(falls_sbox0=per[0], falls_sbox1=per[1], cross_excess=cross)
            print(f"   => DIRECT-SUM TEST: falls(both) - falls(S0) - falls(S1) = {falls} - {per[0]} - "
                  f"{per[1]} = {cross}  -> {'DIRECT SUM' if cross == 0 else 'CROSS-S-BOX FALLS'}")
    print(f"   total {time.time()-t_all:.1f}s")
    record("E12_first_layer_ideal.csv", rec, seconds=time.time() - t_all, censored=cens)
