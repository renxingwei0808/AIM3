#!/usr/bin/env python3
"""Gate for the Sage/M4RI rank backend (CLAUDE.md ground rule / batch task 2 + task 0.3).

Nothing in this batch may use `--backend sage` at a size beyond Tier-A reach until this script
prints `BACKEND ALIGNED`.  It checks two independent things:

  A. python vs sage on the SAME instance, for the cases the batch prompt names
     (E03 n=11 e=3 and n=13 e=7, E12 n=11 --D4, E06 n=7 D=4) — every reported number must agree
     exactly: rows, cols, rank, syz, dim(span n deg<=D-1), and for E03 also the point-set ideal
     dimension `dim I(V)_{<=3}`.
  B. sage vs the committed Tier-A CSV rows for the same (n, e) / (n) — this catches an instance
     that is reproducible but different (a changed random stream) as well as a wrong rank.

Usage:
  python3 server/jobs/verify_backend.py [--only A|B] [--json OUT.json] [--workdir DIR]

Exit status 0 iff every comparison agreed.
"""
import argparse
import csv
import json
import os
import sys
import time
from math import gcd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP] + \
         [os.path.join(_EXP, d) for d in sorted(os.listdir(_EXP)) if d.startswith("E")]:
    if os.path.isdir(p) and p not in sys.path:
        sys.path.insert(0, p)

from aim3.boolmac import macaulay_stats                                    # noqa: E402
from aim3.boolmac_export import macaulay_stats_backend, relation_dim_sage  # noqa: E402

CHECKS = []


def record(name, quantity, python, sage, extra=""):
    ok = python == sage
    CHECKS.append(dict(check=name, quantity=quantity, python=python, sage=sage, ok=ok, extra=extra))
    print(f"  {'OK  ' if ok else 'FAIL'} {name:28s} {quantity:22s} python={python!s:>10} "
          f"sage={sage!s:>10} {extra}", flush=True)
    return ok


def both(name, eqs, nv, D, workdir):
    """Run one Macaulay matrix through both engines and compare every field."""
    t0 = time.time()
    a = macaulay_stats(eqs, nv, D)
    ta = time.time() - t0
    t0 = time.time()
    b = macaulay_stats_backend(eqs, nv, D, backend="sage", workdir=workdir, verbose=False)
    tb = time.time() - t0
    ok = True
    for field in ("rows", "cols", "rank", "syz", "dim_low"):
        ok &= record(name, f"D={D} {field}", getattr(a, field), getattr(b, field),
                     extra=f"({ta:.1f}s vs {tb:.1f}s)" if field == "rows" else "")
    return a, b, ok


# --------------------------------------------------------------------------------- A: same instance
def check_E03(n, e, workdir, do_points=True):
    from aim3_lab import F, sbox
    from aim3_sbox_syz import sbox_system
    from sbox_degree_fall import relation_dim
    print(f"\n[A] E03 single S-box n={n} e={e}")
    fld, P1, P2, P3 = sbox_system(n, e)
    allq = P1 + P2 + P3
    nv = 2 * n
    ok = True
    dpy = None
    a3, b3, o = both(f"E03 n={n} e={e}", allq, nv, 3, workdir)
    ok &= o
    a4, b4, o = both(f"E03 n={n} e={e}", allq, nv, 4, workdir)
    ok &= o
    ok &= record(f"E03 n={n} e={e}", "falls_D4", a4.dim_low - a3.rank, b4.dim_low - b3.rank)
    if do_points:
        pts = [x | (sbox(fld, x, e) << n) for x in range(1, fld.q)]
        t0 = time.time()
        dpy, tot, rkpy = relation_dim(pts, nv, 3)
        tpy = time.time() - t0
        t0 = time.time()
        dsg, tots, rksg = relation_dim_sage(pts, nv, 3, workdir=workdir, verbose=False)
        tsg = time.time() - t0
        ok &= record(f"E03 n={n} e={e}", "dim I(V)_<=3", dpy, dsg, extra=f"({tpy:.1f}s vs {tsg:.1f}s)")
        ok &= record(f"E03 n={n} e={e}", "point-matrix rank", rkpy, rksg)
    return ok, dict(n=n, e=e, syz3=a3.syz, falls_D4=a4.dim_low - a3.rank,
                    mac3_rank=a3.rank, mac4_dim_le3=a4.dim_low,
                    IV_deg3=dpy)


def check_E12(n, workdir):
    from first_layer_ideal import build_first_layer
    print(f"\n[A] E12 first-layer system n={n} --D4")
    es = [e for e in range(3, n - 1) if gcd(e, n) == 1][:2]
    nv, eqsP1, eqsAll, sol, pts = build_first_layer(n, es, 1)
    ok = True
    a3, b3, o = both(f"E12 n={n}", eqsAll, nv, 3, workdir)
    ok &= o
    a4, b4, o = both(f"E12 n={n}", eqsAll, nv, 4, workdir)
    ok &= o
    ok &= record(f"E12 n={n}", "falls_D4", a4.dim_low - a3.rank, b4.dim_low - b3.rank)
    per_py, per_sg = [], []
    for i in range(2):
        sub = eqsAll[i * 3 * n:(i + 1) * 3 * n]
        s3 = macaulay_stats(sub, nv, 3)
        s4 = macaulay_stats(sub, nv, 4)
        g3 = macaulay_stats_backend(sub, nv, 3, backend="sage", workdir=workdir, verbose=False)
        g4 = macaulay_stats_backend(sub, nv, 4, backend="sage", workdir=workdir, verbose=False)
        per_py.append(s4.dim_low - s3.rank)
        per_sg.append(g4.dim_low - g3.rank)
        ok &= record(f"E12 n={n}", f"falls_D4 S-box{i} alone", per_py[-1], per_sg[-1])
    ok &= record(f"E12 n={n}", "cross_excess",
                 (a4.dim_low - a3.rank) - sum(per_py), (b4.dim_low - b3.rank) - sum(per_sg))
    return ok, dict(n=n, exps=es, mac3_rank=a3.rank, mac3_syz=a3.syz,
                    falls_D4=a4.dim_low - a3.rank, falls_sbox0=per_py[0], falls_sbox1=per_py[1],
                    cross_excess=(a4.dim_low - a3.rank) - sum(per_py))


def check_E06(n, D, workdir):
    from kernel_bilinear_toy import build_kernel_model
    print(f"\n[A] E06 kernel bilinear model n={n} D={D}")
    exps = [e for e in range(2, n) if gcd(e, n) == 1][:3]
    nv, eqs, sol = build_kernel_model(n, exps)
    ok = True
    rows = {}
    for d in range(2, D + 1):
        a, b, o = both(f"E06 n={n}", eqs, nv, d, workdir)
        ok &= o
        rows[d] = a
    return ok, dict(n=n, exps=exps, per_D={d: dict(rank=r.rank, syz=r.syz, dim_low=r.dim_low)
                                           for d, r in rows.items()})


# ------------------------------------------------------------------------ B: against committed CSVs
def csv_rows(name):
    path = os.path.join(ROOT, "results", name)
    if not os.path.exists(path):
        print(f"  (missing {path})")
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def check_against_csv(measured):
    print("\n[B] measured values vs the committed Tier-A CSV rows")
    ok = True
    for m in measured.get("E03", []):
        for r in csv_rows("E03_sbox_degree_fall.csv"):
            if int(r["n"]) == m["n"] and int(r["e"]) == m["e"] and r["censored"] != "1":
                for key, col in (("syz3", "syz3"), ("falls_D4", "falls_D4"),
                                 ("mac3_rank", "mac3_rank"), ("IV_deg3", "IV_deg3")):
                    if m.get(key) is None or not r[col]:
                        continue
                    ok &= record(f"E03 csv n={m['n']} e={m['e']}", col, int(r[col]), m[key])
                break
    for m in measured.get("E12", []):
        for r in csv_rows("E12_first_layer_ideal.csv"):
            if int(r["n"]) == m["n"] and r["censored"] != "1" and r["falls_D4"]:
                for key, col in (("mac3_rank", "mac3_rank"), ("mac3_syz", "mac3_syz"),
                                 ("falls_D4", "falls_D4"), ("cross_excess", "cross_excess")):
                    if not r[col]:
                        continue
                    ok &= record(f"E12 csv n={m['n']}", col, int(r[col]), m[key])
                break
    for m in measured.get("E06", []):
        rows = [r for r in csv_rows("E06_kernel_model_falls.csv")
                if int(r["n"]) == m["n"] and r["system"] == "kernel" and r["censored"] != "1"]
        for r in rows:
            d = int(r["D"])
            if d in m["per_D"]:
                for key, col in (("rank", "rank"), ("syz", "syz"), ("dim_low", "dim_low")):
                    ok &= record(f"E06 csv n={m['n']} D={d}", col, int(r[col]), m["per_D"][d][key])
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["A", "B"], default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--workdir", default=os.environ.get("AIM3_SCRATCH"))
    ap.add_argument("--skip-e03-points", action="store_true",
                    help="skip dim I(V)_<=3 (the O(2^n) point matrix) in the E03 checks")
    a = ap.parse_args()
    if a.workdir:
        os.makedirs(a.workdir, exist_ok=True)
    print(f"# verify_backend: root={ROOT} sage={os.environ.get('AIM3_SAGE', 'sage')} "
          f"workdir={a.workdir}")
    ok = True
    measured = {"E03": [], "E12": [], "E06": []}
    if a.only != "B":
        for n, e in ((11, 3), (13, 7)):
            o, m = check_E03(n, e, a.workdir, do_points=not a.skip_e03_points)
            ok &= o
            measured["E03"].append(m)
        o, m = check_E12(11, a.workdir)
        ok &= o
        measured["E12"].append(m)
        o, m = check_E06(7, 4, a.workdir)
        ok &= o
        measured["E06"].append(m)
    if a.only != "A":
        ok &= check_against_csv(measured)
    nfail = sum(1 for c in CHECKS if not c["ok"])
    print(f"\n{len(CHECKS)-nfail}/{len(CHECKS)} comparisons agreed")
    print("BACKEND ALIGNED" if ok else "BACKEND NOT ALIGNED — do not use --backend sage")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(dict(ok=ok, checks=CHECKS, measured=measured), fh, indent=1, default=str)
        print(f"written {a.json}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
