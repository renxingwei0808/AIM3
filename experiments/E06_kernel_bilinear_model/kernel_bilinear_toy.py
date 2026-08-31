#!/usr/bin/env python3
"""Symbolic kernel-variable bilinear model of toy AIM3 (l=2), Focus Report §4.3:
  variables U = (pt, t0, t1) [3n bits],  W = (a0, a1, a2) [3n bits]
  B_{e_i}(a_i, A_i pt + b_i) = a_i ,  t_i a_i = a_i^{2^{e_i}}   (i=0,1)
  B_{e_2}(a_2, b_2 + B_0 t_0 + B_1 t_1) = a_2 ,  (pt+ct) a_2 = a_2^{2^{e_2}}
6n Boolean equations, all bilinear in (U,W). Verifies vanishing at the key and reports Macaulay ranks / degree falls,
side by side with S_quad and P1-only on the same instance.
Also builds a RANDOM BILINEAR CONTROL of exactly the same shape (`random_bilinear_system`): 6n equations,
variables split into two blocks U (3n bits) and W (3n bits), every equation a random U x W bilinear form
plus a random linear part plus the constant that plants a random root.  This isolates what in the kernel
model's degree-fall behaviour is due to bilinearity alone and what is due to the AIM3 structure.

Usage: python3 kernel_bilinear_toy.py n [Dmax] [--no-control] [--budget SEC]
                                        [--backend python|sage] [--export DIR]
       (n=5,7 in seconds; n=9 D<=3; larger -> server).  Appends to results/E06_kernel_model_falls.csv.
  --backend  rank engine; 'sage' routes every Macaulay matrix through M4RI (verified identical to
             the pure-Python engine by server/jobs/verify_backend.py).  Scratch: $AIM3_SCRATCH.
  --export   also write the kernel model (and the systems it is compared against) to DIR in msolve
             and Magma BooleanPolynomialRing syntax, for the F4 step-degree runs of E04 task 6.
"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS, record
from math import gcd
from aim3_lab import F, sbox, rand_unit_lu, apply_mat
from aim3_macaulay import AffExpr, prod_coords, lin_coords, add_polys, build_aim3, macaulay, random_system
from aim3.boolmac import macaulay_stats, MacaulayTimeout
from aim3.boolmac_export import macaulay_stats_backend
from aim3.polysys_export import write_system

SCRATCH = os.environ.get("AIM3_SCRATCH")


def random_bilinear_system(n, seed=3):
    """Control group for the kernel model: 6n Boolean equations in 6n variables split U = 0..3n-1,
    W = 3n..6n-1; each equation is a random U x W bilinear form + a random affine part, with the
    constant chosen so that a random planted point is a root.  Same shape as the kernel model
    (same #vars, #eqs, same bipartite support), no AIM3 structure."""
    rng = random.Random(seed)
    nv = 6 * n
    sol = [rng.getrandbits(1) for _ in range(nv)]
    eqs = []
    for _ in range(6 * n):
        p = set()
        for i in range(3 * n):
            for j in range(3 * n, 6 * n):
                if rng.getrandbits(1):
                    p ^= {frozenset((i, j))}
        for i in range(nv):
            if rng.getrandbits(1):
                p ^= {frozenset((i,))}
        s = 0
        for m in p:
            v = 1
            for j in m:
                v &= sol[j]
            s ^= v
        if s:
            p ^= {frozenset()}
        eqs.append(p)
    return nv, eqs, sol

def build_kernel_model(n, exps, seed=1):
    ell = len(exps) - 1; assert ell == 2
    fld = F(n); rng = random.Random(seed)
    mats = [rand_unit_lu(n, rng) for _ in range(2 * ell)]
    consts = [rng.randrange(fld.q) for _ in range(ell + 1)]
    while True:
        pt = rng.randrange(1, fld.q)
        xs = [apply_mat(mats[i], pt, n) ^ consts[i] for i in range(ell)]
        if 0 in xs: continue
        ts = [sbox(fld, xs[i], exps[i]) for i in range(ell)]
        xl = consts[ell]
        for i in range(ell): xl ^= apply_mat(mats[ell + i], ts[i], n)
        if xl == 0 or 0 in ts: continue
        tl = sbox(fld, xl, exps[ell])
        if tl == 0: continue
        ct = pt ^ tl; break
    rs = [fld.pw(t, pow((1 << e) - 1, -1, fld.q - 1)) for t, e in zip(ts + [tl], exps)]   # kernel values a_i
    # variable layout: pt 0..n-1 | t0 n..2n-1 | t1 2n..3n-1 | a0 3n.. | a1 4n.. | a2 5n..
    nv = 6 * n
    sol = [0] * nv
    for k in range(n):
        sol[k] = (pt >> k) & 1
        for i in range(2): sol[(1 + i) * n + k] = (ts[i] >> k) & 1
        for i in range(3): sol[(3 + i) * n + k] = (rs[i] >> k) & 1
    def var(block): return AffExpr(fld, {block * n + k: 1 << k for k in range(n)}, 0)
    def colimgs(A): return [apply_mat(A, 1 << j, n) for j in range(n)]
    PT, T0, T1, A0, A1, A2 = (var(b) for b in range(6))
    eqs = []
    for i, (Ai, Ti) in enumerate(((A0, T0), (A1, T1))):
        ci = colimgs(mats[i]); e = exps[i]
        xi = AffExpr(fld, {j: ci[j] for j in range(n)}, consts[i])
        Bi = add_polys(prod_coords(Ai, xi.frob(e)), prod_coords(Ai.frob(e), xi), lin_coords(Ai))   # B(a_i,x_i)+a_i
        Ki = add_polys(prod_coords(Ti, Ai), lin_coords(Ai.frob(e)))                                    # t_i a_i + a_i^{2^e}
        eqs += Bi + Ki
    xl = AffExpr(fld, {}, consts[2])
    for i, Ti in enumerate((T0, T1)):
        ci = colimgs(mats[2 + i])
        for k in range(n):
            xl.coefs[(1 + i) * n + k] = xl.coefs.get((1 + i) * n + k, 0) ^ ci[k]
    e = exps[2]
    tl = AffExpr(fld, {k: 1 << k for k in range(n)}, ct)
    eqs += add_polys(prod_coords(A2, xl.frob(e)), prod_coords(A2.frob(e), xl), lin_coords(A2))
    eqs += add_polys(prod_coords(tl, A2), lin_coords(A2.frob(e)))
    def ev(p):
        s = 0
        for m in p:
            v = 1
            for j in m: v &= sol[j]
            s ^= v
        return s
    assert all(ev(p) == 0 for p in eqs), "kernel model does not vanish at the key"
    return nv, eqs, sol

def measure(name, eqs, nv, D, budget, backend="python"):
    t0 = time.time()
    try:
        return macaulay_stats_backend(eqs, nv, D, backend=backend, budget=budget, workdir=SCRATCH,
                                      tag=f"e06_{name}_nv{nv}"), 0, time.time() - t0
    except MacaulayTimeout as exc:
        print(f"     !! {name} D={D} exceeded {budget:.0f}s -> lower bounds only")
        return exc.partial, 1, time.time() - t0


if __name__ == "__main__":
    n = int(sys.argv[1])
    Dmax = int(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else 3
    budget = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else 3600.0
    backend = sys.argv[sys.argv.index("--backend") + 1] if "--backend" in sys.argv else "python"
    export_dir = sys.argv[sys.argv.index("--export") + 1] if "--export" in sys.argv else None
    # e up to n-1 (the spec only asks gcd(e,n)=1); needed to get three exponents at n=5.
    # For n=7,9 this is the same list [2,3,4] / [2,4,5] the script used before.
    exps = [e for e in range(2, n) if gcd(e, n) == 1][:3]
    print(f"# params: n={n} Dmax={Dmax} exps={exps} control={'--no-control' not in sys.argv} "
          f"budget={budget}s backend={backend} export={export_dir}")
    nv, eqs, sol = build_kernel_model(n, exps)
    print(f"kernel bilinear model: n={n} exps={exps}: {len(eqs)} eqs in {nv} vars (6n) - vanishes at key: OK")
    nv2, eqsP1, eqsAll, _ = build_aim3(n, exps)
    systems = [("kernel", "kernel-bilinear (6n eq / 6n var)", eqs, nv),
               ("S_quad", "S_quad (9n / 3n)", eqsAll, nv2),
               ("P1only", "P1-only (3n / 3n)", eqsP1, nv2)]
    if "--no-control" not in sys.argv:
        nv3, reqs, _ = random_bilinear_system(n)
        systems.append(("random_bilinear", "random bilinear control (6n / 6n)", reqs, nv3))
        rq, _ = random_system(nv2, len(eqsAll), seed=7)
        systems.append(("random_quadratic", "random quadratic control (9n / 3n)", rq, nv2))
    if export_dir:
        for tag, name, E, V in systems:
            info = write_system(E, V, export_dir, f"E06_n{n}_{tag}", sysname=tag,
                                sol=sol if tag == "kernel" else None,
                                meta=dict(experiment="E06", n=n, ell=2, exps=exps, system=tag,
                                          seed=1))
            print(f"  [export] {tag}: {info['neqs']} eqs / {info['nvars']} vars, "
                  f"{info['terms']} terms -> {export_dir}/E06_n{n}_{tag}.{{ms,magma.m,json}}")
    for tag, name, E, V in systems:
        print(f"  {name}")
        prev_rank = None
        for D in range(2, Dmax + 1):
            st, cens, sc = measure(tag, E, V, D, budget, backend)
            falls = "" if prev_rank is None else st.dim_low - prev_rank
            print(f"     D={D}: rows={st.rows:6d} cols={st.cols:7d} rank={st.rank:6d} syz={st.syz:5d} "
                  f"dim(span\u2229deg<={D-1})={st.dim_low}"
                  + (f"  falls={falls} ({falls/n:.2f} n)" if falls != "" else "")
                  + ("   [CENSORED: lower bound]" if cens else ""))
            record("E06_kernel_model_falls.csv",
                   dict(n=n, exps=";".join(map(str, exps)), system=tag, D=D, nvars=V, neqs=len(E),
                        rows=st.rows, cols=st.cols, rank=st.rank, syz=st.syz, dim_low=st.dim_low,
                        falls_from_prev_D=falls, twice_neqs=2 * len(E),
                        dim_low_over_n=round(st.dim_low / n, 4)),
                   seconds=sc, censored=cens)
            if cens:
                break
            prev_rank = st.rank
