#!/usr/bin/env python3
"""Gate 1, localised: the D=4 degree falls of every SUBSET of the S-box blocks (task 3).

`full_system_falls.py --per-block` answers one question — is `falls(all)` equal to
`sum_i falls(S-box i alone)`?  When the answer is no, the batch asks *which pairs* are responsible
and whether the excess grows with n.  This script measures the whole subset lattice instead:

    for every non-empty subset B of the (ell+1) S-box blocks
        falls(B) = dim(span_{D=4}(eqs of B) n deg<=3) - rank_{D=3}(eqs of B)

from which

    cross_excess(B) = falls(B) - sum_{i in B} falls({i})          (0 <=> direct sum on B)
    pair_excess(i,j) = falls({i,j}) - falls({i}) - falls({j})

All of it is run in the SAME (ell+1)n variables, so the cross-block multipliers are present for
every subset and cancel out of the comparison exactly as in `full_system_falls.py`.

**The control matters more than the measurement.**  A same-shape `random` quadratic system, cut
into blocks of 3n equations at the same positions, is put through the identical lattice.  A
non-zero `cross_excess` is evidence about AIM3 only if the random control's is zero: otherwise it
only says that more equations reach further into the (huge) cubic part of the ideal of a single
point.  `saturation = dim_low / C(nv,<=3)` is reported for every subset; the batch rule is that
only runs below 0.9 enter a verdict.

Usage:
  python3 server/jobs/gate1_pairs.py N -l ELL [--seed 1] [--systems S_quad,random]
                                     [--max-size K] [--backend sage] [--budget SEC]
                                     [--csv E04_server_gate1_pairs.csv]
"""
import argparse
import itertools
import os
import sys
import time
from math import comb, gcd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP,
          os.path.join(_EXP, "E04_full_system_hilbert")]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record                                          # noqa: E402
from aim3_macaulay import build_aim3, random_system                 # noqa: E402
from aim3.boolmac import MacaulayTimeout                            # noqa: E402
from aim3.boolmac_export import macaulay_stats_backend              # noqa: E402
from full_system_falls import default_exps                          # noqa: E402

SCRATCH = os.environ.get("AIM3_SCRATCH")


def falls_of(eqs, nv, backend, budget, tag):
    """(falls_D4, dim_low4, rank3, syz3, censored, seconds) for one equation set."""
    t0 = time.time()
    try:
        s3 = macaulay_stats_backend(eqs, nv, 3, backend=backend, budget=budget, workdir=SCRATCH,
                                    tag=f"{tag}_D3", verbose=False)
        s4 = macaulay_stats_backend(eqs, nv, 4, backend=backend, budget=budget, workdir=SCRATCH,
                                    tag=f"{tag}_D4", verbose=False)
    except MacaulayTimeout:
        return None, None, None, None, 1, time.time() - t0
    return s4.dim_low - s3.rank, s4.dim_low, s3.rank, s3.syz, 0, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("-l", "--ell", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--systems", default="S_quad,random")
    ap.add_argument("--exps", default=None,
                    help="comma-separated exponents (default: full_system_falls.default_exps, "
                         "which always starts at e=2 -- a DEGENERATE class. Pass generic ones to "
                         "check that the verdict does not depend on that choice)")
    ap.add_argument("--max-size", type=int, default=0,
                    help="largest subset size to measure (0 = all, i.e. up to ell+1)")
    ap.add_argument("--backend", choices=["python", "sage"], default="sage")
    ap.add_argument("--budget", type=float, default=86400.0)
    ap.add_argument("--csv", default="E04_server_gate1_pairs.csv")
    a = ap.parse_args()

    exps = [int(v) for v in a.exps.split(",")] if a.exps else default_exps(a.n, a.ell)
    if len(exps) != a.ell + 1:
        raise SystemExit(f"need {a.ell+1} exponents for ell={a.ell}, got {exps}")
    nv, eqsP1, eqsAll, _ = build_aim3(a.n, exps, seed=a.seed)
    nblocks = a.ell + 1
    blocksize = 3 * a.n
    low_cols = sum(comb(nv, d) for d in range(4))
    maxsize = a.max_size or nblocks
    pool = {"S_quad": eqsAll}
    if "random" in a.systems:
        pool["random"] = random_system(nv, len(eqsAll), seed=1000 + a.seed)[0]
    print(f"# params: n={a.n} ell={a.ell} seed={a.seed} exps={exps} backend={a.backend} "
          f"max_size={maxsize} budget={a.budget}s")
    print(f"{nv} variables, {nblocks} blocks of {blocksize} equations, "
          f"degree-<=3 columns C({nv},<=3) = {low_cols}")

    for sysname in a.systems.split(","):
        if sysname not in pool:
            continue
        eqs = pool[sysname]
        print(f"\n== {sysname}: {len(eqs)} equations")
        single = {}
        for size in range(1, maxsize + 1):
            for sub in itertools.combinations(range(nblocks), size):
                part = [e for i in sub for e in eqs[i * blocksize:(i + 1) * blocksize]]
                tag = f"g1_{sysname}_{''.join(map(str, sub))}"
                f, low, r3, syz3, cens, secs = falls_of(part, nv, a.backend, a.budget, tag)
                if size == 1:
                    single[sub[0]] = f
                sum_single = (None if any(single.get(i) is None for i in sub)
                              else sum(single[i] for i in sub))
                cross = None if (f is None or sum_single is None) else f - sum_single
                sat = None if low is None else round(low / low_cols, 4)
                label = "+".join(map(str, sub))
                print(f"   blocks {{{label}}}: {len(part):4d} eqs  D3 rank={r3} syz={syz3} "
                      f"(8n*{size}={8*a.n*size})  falls_D4={f}  sum_singles={sum_single}  "
                      f"cross_excess={cross}  saturation={sat}"
                      + ("   [CENSORED: lower bound]" if cens else ""), flush=True)
                record(a.csv, dict(
                    n=a.n, ell=a.ell, seed=a.seed, exps=";".join(map(str, exps)), system=sysname,
                    blocks=label, nblocks=len(sub), nvars=nv, neqs=len(part),
                    D3_rank=r3 if r3 is not None else "", D3_syz=syz3 if syz3 is not None else "",
                    D4_dim_low=low if low is not None else "",
                    falls_D4=f if f is not None else "",
                    sum_singles=sum_single if sum_single is not None else "",
                    cross_excess=cross if cross is not None else "",
                    low_cols=low_cols, saturation=sat if sat is not None else ""),
                    seconds=secs, censored=cens, tier="B")


if __name__ == "__main__":
    main()
