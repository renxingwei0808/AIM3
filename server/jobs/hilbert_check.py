#!/usr/bin/env python3
"""Batch-2 task 1.4: measured Macaulay ranks against the Boolean semi-regular Hilbert function.

B2 assumes `S_quad` behaves like a semi-regular system of m = 3n(ℓ+1) quadratics in N = (ℓ+1)n
Boolean variables.  That assumption has a Hilbert series,

    HS(z) = [ (1+z)^N / (1+z²)^m ]_+          (truncated at its first non-positive coefficient)

whose partial sums give the dimension of the degree-≤D part of the quotient, hence the dimension of
the degree-≤D part of the ideal, hence the rank a degree-D Macaulay matrix should have if the
assumption held:

    rank_pred(D) = C(N,≤D) − Σ_{d≤D} [z^d] HS(z),      syz_pred(D) = rows(D) − rank_pred(D)

This script prints that prediction next to the measured rank, for whatever the CSVs already hold and
for any (n, ℓ, D) given on the command line (which it then measures).  The gap is the quantitative
form of the ROADMAP's C3/C4 objection: `syz₃ = 8n(ℓ+1)` is already a departure at D = 3, and the
question is whether it grows.

Usage:
  python3 server/jobs/hilbert_check.py --n 11 --ell 2 --D 5 [--system S_quad] [--backend sage]
                                       [--exps 3,4,7] [--budget SEC]
  python3 server/jobs/hilbert_check.py --from-csv          (no new computation: use what is stored)
"""
import argparse
import csv
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
from full_system_falls import default_exps, generic_exps            # noqa: E402

SCRATCH = os.environ.get("AIM3_SCRATCH")


def boolean_semiregular_series(N, m, maxd):
    """[z^d] of (1+z)^N / (1+z^2)^m for d = 0..maxd, truncated at the first non-positive term."""
    out = []
    for d in range(maxd + 1):
        c = 0
        for k in range(d // 2 + 1):
            c += (-1) ** k * comb(m + k - 1, k) * comb(N, d - 2 * k)
        if c <= 0:
            break
        out.append(c)
    return out


def predicted_rank(N, m, D):
    """(rank_pred, quotient_dim_<=D, series) for a semi-regular system, at Macaulay degree D."""
    hs = boolean_semiregular_series(N, m, D)
    quot = sum(hs[: D + 1])
    cols = sum(comb(N, i) for i in range(D + 1))
    return cols - quot, quot, hs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int)
    ap.add_argument("--ell", type=int, default=2)
    ap.add_argument("--D", type=int, default=5)
    ap.add_argument("--system", default="S_quad", choices=["S_quad", "P1only", "random"])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exps", default=None)
    ap.add_argument("--generic-exps", action="store_true")
    ap.add_argument("--backend", choices=["python", "sage"], default="sage")
    ap.add_argument("--budget", type=float, default=86400.0)
    ap.add_argument("--from-csv", action="store_true",
                    help="only tabulate what results/E04_server_gate1.csv already holds")
    ap.add_argument("--csv", default="E04_server_hilbert.csv")
    a = ap.parse_args()

    if a.from_csv:
        p = os.path.join(ROOT, "results", "E04_server_gate1.csv")
        print(f"{'l':>2} {'n':>3} {'sys':<8} {'D':>2} {'rows':>8} {'rank':>8} {'pred':>8} "
              f"{'rank-pred':>10} {'syz':>8} {'syz_pred':>9}")
        for r in csv.DictReader(open(p)):
            if not r["D3_rank"] or r["censored"] == "1":
                continue
            N, m = int(r["nvars"]), int(r["neqs"])
            for D, rk, rows, syz in ((3, r["D3_rank"], r["D3_rows"], r["D3_syz"]),
                                     (4, r["D4_rank"], r["D4_rows"], r["D4_syz"])):
                if not rk:
                    continue
                pr, _, _ = predicted_rank(N, m, D)
                print(f"{r['ell']:>2} {r['n']:>3} {r['system']:<8} {D:>2} {int(rows):>8} "
                      f"{int(rk):>8} {pr:>8} {int(rk)-pr:>10} {int(syz):>8} "
                      f"{int(rows)-pr:>9}")
        return

    if a.n is None:
        raise SystemExit("--n is required unless --from-csv is given")
    n, ell, D = a.n, a.ell, a.D
    if a.exps:
        exps = [int(v) for v in a.exps.split(",")]
    else:
        exps = generic_exps(n, ell) if a.generic_exps else default_exps(n, ell)
    nv, eqsP1, eqsAll, _ = build_aim3(n, exps, seed=a.seed)
    eqs = {"S_quad": eqsAll, "P1only": eqsP1,
           "random": random_system(nv, len(eqsAll), seed=1000 + a.seed)[0]}[a.system]
    N, m = nv, len(eqs)
    cols = sum(comb(N, i) for i in range(D + 1))
    pr, quot, hs = predicted_rank(N, m, D)
    print(f"# params: n={n} ell={ell} exps={exps} system={a.system} D={D} backend={a.backend}")
    print(f"N = {N} variables, m = {m} quadratics; C(N,<={D}) = {cols}")
    print(f"semi-regular series coefficients [z^0..]: {hs}")
    print(f"predicted quotient dim <= {D}: {quot};  predicted rank: {pr}")

    t0 = time.time()
    cens = 0
    try:
        st = macaulay_stats_backend(eqs, nv, D, backend=a.backend, budget=a.budget,
                                    workdir=SCRATCH, tag=f"hil_n{n}_l{ell}_{a.system}_D{D}")
    except MacaulayTimeout as exc:
        st, cens = exc.partial, 1
    secs = time.time() - t0
    print(f"measured: rows={st.rows} cols={st.cols} rank={st.rank} syz={st.syz} "
          f"dim_low={st.dim_low}   ({secs:.1f}s)" + ("   [CENSORED]" if cens else ""))
    if not cens:
        print(f"   rank − predicted = {st.rank - pr};  "
              f"syz − predicted = {st.syz - (st.rows - pr)}")
    record(a.csv, dict(n=n, ell=ell, exps=";".join(map(str, exps)), seed=a.seed, system=a.system,
                       D=D, nvars=N, neqs=m, cols=cols,
                       rows=st.rows, rank=st.rank, syz=st.syz, dim_low=st.dim_low,
                       predicted_rank=pr, predicted_syz=st.rows - pr,
                       rank_minus_predicted=st.rank - pr if not cens else "",
                       syz_minus_predicted=st.syz - (st.rows - pr) if not cens else "",
                       semiregular_series=";".join(map(str, hs))),
           seconds=secs, censored=cens, tier="B")


if __name__ == "__main__":
    main()
