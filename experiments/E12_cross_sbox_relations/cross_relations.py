#!/usr/bin/env python3
"""Cross-S-box relations: dimension of Boolean relations of degree <= deg vanishing on the joint graph
{(t0,t1)} and {(pt,t0,t1)} of a toy AIM3 first layer (unit-LU matrices). The spec only excludes *quadratic*
cross relations; this measures cubic (and higher) ones. Meaningful when 2^n >> C(#vars, <=deg).
Usage: python3 cross_relations.py n deg [seed] [--backend python|sage]
       e.g.  n=14 deg=3 (2n=28 vars: 3683 monomials vs 16383 points)
`--backend sage` computes the point-matrix rank with M4RI (identical output, verified by
server/jobs/verify_backend.py); at n=16 deg=4 the matrix is 65534 x 41449 and pure Python is
out of reach.
"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS, record
from math import gcd, comb
from aim3_lab import F, sbox, rand_unit_lu, apply_mat, rank_gf2
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "E03_syzygy_and_degree_fall"))
from sbox_degree_fall import relation_dim
if __name__ == "__main__":
    backend = sys.argv[sys.argv.index("--backend") + 1] if "--backend" in sys.argv else "python"
    pos = [a for a in sys.argv[1:] if not a.startswith("--") and a not in ("python", "sage")]
    n, deg = int(pos[0]), int(pos[1]); seed = int(pos[2]) if len(pos) > 2 else 1
    print(f"# params: n={n} deg={deg} seed={seed} backend={backend}")
    t0 = time.time()
    rng = random.Random(seed); fld = F(n)
    es = [e for e in range(3, n - 1) if gcd(e, n) == 1][:2]
    A0, A1 = rand_unit_lu(n, rng), rand_unit_lu(n, rng); b0, b1 = rng.randrange(fld.q), rng.randrange(fld.q)
    pts_t, pts_all = [], []
    for pt in range(fld.q):
        x0 = apply_mat(A0, pt, n) ^ b0; x1 = apply_mat(A1, pt, n) ^ b1
        if x0 == 0 or x1 == 0: continue
        t0_, t1_ = sbox(fld, x0, es[0]), sbox(fld, x1, es[1])
        pts_t.append(t0_ | (t1_ << n)); pts_all.append(pt | (t0_ << n) | (t1_ << (2 * n)))
    d_t, tot_t, _ = relation_dim(pts_t, 2 * n, deg, backend=backend)
    print(f"n={n} e={es} deg<={deg}: (t0,t1): dim={d_t} of {tot_t} monomials on {len(pts_t)} points  "
          f"[expected 0 if no cross relations; note the trace-image constraints Tr(t_i^-r)=0 have degree wt(-r)]")
    rec = dict(n=n, exps=";".join(map(str, es)), deg=deg, seed=seed, npoints=len(pts_t),
               dim_t0t1=d_t, mono_t0t1=tot_t, dim_pt_t0_t1="", mono_pt_t0_t1="")
    if 3 * n <= 48 and comb(3 * n, deg) < 2 ** n:
        d_a, tot_a, _ = relation_dim(pts_all, 3 * n, deg, backend=backend)
        rec.update(dim_pt_t0_t1=d_a, mono_pt_t0_t1=tot_a)
        print(f"   (pt,t0,t1): dim={d_a} of {tot_a}   [degree-3 span of the two S-box ideals predicts 2*(6n^2-5n)+? ; excess = cross falls]")
    print(f"   {time.time()-t0:.1f}s")
    record("E12_cross_relations.csv", rec, seconds=time.time() - t0)
