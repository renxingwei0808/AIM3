#!/usr/bin/env python3
"""E03 at D = 5 for a single S-box (task 4, second half).

E03 established the first-fall degree 4 and the 2n new cubics it produces.  The batch asks whether
D = 5 produces *further* falls — new elements of degree ≤ 4 beyond what D = 4 already reaches —
which is the missing data point for T1.2.  Definitions match E03/E04 exactly:

    falls_D4 = dim(span_{D=4} ∩ deg≤3) − rank_{D=3}
    falls_D5 = dim(span_{D=5} ∩ deg≤4) − rank_{D=4}          <-- what this script adds

`saturation_D5 = dim(span_{D=5} ∩ deg≤4) / C(2n,≤4)` is reported for the same reason as in E04: at
saturation the number stops measuring structure.  The same-shape random quadratic control is run on
request (`--control`), because a fall count is only evidence when the control has none.

`server/SIZING.md` puts the dense D=5 matrix at ≈ 2.4 GB for n=14 and ≈ 8 GB for n=16; both are run
with the Sage/M4RI backend.

Usage:
  python3 server/jobs/e03_d5.py N [--exps 3,5] [--backend sage] [--budget SEC] [--control]
                                  [--csv E03_server_D5.csv]
"""
import argparse
import os
import sys
import time
from math import comb, gcd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record                                          # noqa: E402
from aim3_sbox_syz import sbox_system                               # noqa: E402
from aim3_macaulay import random_system                             # noqa: E402
from aim3.boolmac import MacaulayTimeout                            # noqa: E402
from aim3.boolmac_export import macaulay_stats_backend              # noqa: E402
from aim3.structures import predicted_falls_D4, degenerate_class    # noqa: E402

SCRATCH = os.environ.get("AIM3_SCRATCH")


def measure(eqs, nv, D, backend, budget, tag):
    t0 = time.time()
    try:
        st = macaulay_stats_backend(eqs, nv, D, backend=backend, budget=budget, workdir=SCRATCH,
                                    tag=tag, verbose=True)
        return st, 0, time.time() - t0
    except MacaulayTimeout as exc:
        print(f"   !! D={D} exceeded the {budget:.0f}s budget -> lower bound only", flush=True)
        return exc.partial, 1, time.time() - t0


def run(n, e, backend, budget, control, csvname):
    nv = 2 * n
    fld, P1, P2, P3 = sbox_system(n, e)
    eqs = P1 + P2 + P3
    cols4 = sum(comb(nv, d) for d in range(5))
    print(f"\n== n={n} e={e} ({3*n} quadratics in {nv} variables); "
          f"degree-<=4 columns C({nv},<=4) = {cols4}", flush=True)
    systems = [("S_box", eqs)]
    if control:
        systems.append(("random", random_system(nv, len(eqs), seed=1000 + e)[0]))
    for name, E in systems:
        print(f"  {name}: {len(E)} eqs", flush=True)
        cens = 0
        st = {}
        secs = 0.0
        for D in (3, 4, 5):
            st[D], c, s = measure(E, nv, D, backend, budget, f"e03d5_{name}_n{n}_e{e}")
            cens |= c
            secs += s
            print(f"     D={D}: rows={st[D].rows:7d} cols={st[D].cols:8d} rank={st[D].rank:8d} "
                  f"syz={st[D].syz:7d} dim(span n deg<={D-1})={st[D].dim_low:7d}  ({s:.1f}s)"
                  + ("   [CENSORED: lower bound]" if c else ""), flush=True)
            if c:
                break
        rec = dict(n=n, e=e, system=name, nvars=nv,
                   degenerate=";".join(degenerate_class(n, e)),
                   D3_rank=st[3].rank if 3 in st else "", D3_syz=st[3].syz if 3 in st else "",
                   D4_rank=st[4].rank if 4 in st else "", D4_syz=st[4].syz if 4 in st else "",
                   D4_dim_low=st[4].dim_low if 4 in st else "",
                   D5_rows="", D5_rank="", D5_syz="", D5_dim_low="",
                   falls_D4="", predicted_falls_D4=predicted_falls_D4(n, e) if name == "S_box" else "",
                   falls_D5="", falls_D5_over_n="", cols_le4=cols4, saturation_D5="")
        if 4 in st and not cens:
            falls4 = st[4].dim_low - st[3].rank
            rec["falls_D4"] = falls4
            print(f"     falls_D4 = {st[4].dim_low} - {st[3].rank} = {falls4} ({falls4/n:.2f} n)"
                  + (f"; rule predicts {predicted_falls_D4(n, e)}" if name == "S_box" else ""),
                  flush=True)
        if 5 in st and not cens:
            falls5 = st[5].dim_low - st[4].rank
            sat = st[5].dim_low / cols4
            rec.update(D5_rows=st[5].rows, D5_rank=st[5].rank, D5_syz=st[5].syz,
                       D5_dim_low=st[5].dim_low, falls_D5=falls5,
                       falls_D5_over_n=round(falls5 / n, 4), saturation_D5=round(sat, 4))
            print(f"     falls_D5 = {st[5].dim_low} - {st[4].rank} = {falls5} ({falls5/n:.2f} n); "
                  f"dim(span n deg<=4)/C({nv},<=4) = {sat:.4f}"
                  + ("   [SATURATED: not a structural number]" if sat >= 0.9 else ""), flush=True)
        record(csvname, rec, seconds=secs, censored=cens, tier="B")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--exps", default=None,
                    help="comma-separated exponents (default: the first two admissible, "
                         "one generic and, if present, one degenerate)")
    ap.add_argument("--backend", choices=["python", "sage"], default="sage")
    ap.add_argument("--budget", type=float, default=86400.0)
    ap.add_argument("--control", action="store_true", default=True)
    ap.add_argument("--no-control", dest="control", action="store_false")
    ap.add_argument("--csv", default="E03_server_D5.csv")
    a = ap.parse_args()
    if a.exps:
        exps = [int(v) for v in a.exps.split(",")]
    else:
        adm = [e for e in range(2, a.n - 1) if gcd(e, a.n) == 1]
        generic = [e for e in adm if not degenerate_class(a.n, e)]
        degen = [e for e in adm if degenerate_class(a.n, e)]
        exps = generic[:1] + degen[:1]
    print(f"# params: n={a.n} exps={exps} backend={a.backend} budget={a.budget}s "
          f"control={a.control}")
    for e in exps:
        run(a.n, e, a.backend, a.budget, a.control, a.csv)


if __name__ == "__main__":
    main()
