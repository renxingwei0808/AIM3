#!/usr/bin/env python3
"""Batch-2 task 4: cross-validate the Tier-A T4 measurement in Sage, then scale it past n = 14.

The Tier-A branch (now on `main`) measured the Ore-ring elimination with a hand-written skew-polynomial
class and reached n = 14 (`results/E07_theory_skew_degrees.csv`).  This driver

  1. exports the SAME instance the Tier-A probe builds (same script, same default seed), so the two
     implementations are compared on identical data rather than on two random instances;
  2. runs `experiments/E07_ore_skew_elimination/sage/ore_scale.sage` on it -- Sage's
     `OrePolynomialRing`, an independent implementation of the same algebra;
  3. compares every quantity the Tier-A CSV holds for that n and records agreement/disagreement;
  4. records the new large-n rows (ring degrees to n = 24; the degree in t0 symbolically to n = 24 and
     by evaluation-interpolation as far as the budget reaches).

Usage:
  python3 server/jobs/t4_ore_scale.py --ns 8,9,...,24 [--interp-max 16] [--interp-budget 3600]
                                      [--symbolic-max 24] [--trials 32]
                                      [--csv E07_server_ore_scale.csv]
"""
import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from _common import record                                          # noqa: E402

PROBE = os.path.join(ROOT, "experiments", "E07_ore_skew_elimination", "skew_elimination_probe.py")
SAGE_SCRIPT = os.path.join(ROOT, "experiments", "E07_ore_skew_elimination", "sage", "ore_scale.sage")
TIER_A_CSV = os.path.join(ROOT, "results", "E07_theory_skew_degrees.csv")


def tier_a_rows(n):
    """{measurement: (value, detail)} from the Tier-A CSV for this n, if it has this n at all."""
    import csv
    out = {}
    if not os.path.exists(TIER_A_CSV):
        return out
    for r in csv.DictReader(open(TIER_A_CSV)):
        if int(r["n"]) == n:
            out[r["measurement"]] = (r["value"], r.get("detail", ""), r.get("e0"), r.get("e1"))
    return out


def compare(n, res, ta):
    """Return (checks_done, disagreements[]) against the Tier-A numbers for this n."""
    bad, done = [], 0
    if not ta:
        return 0, bad
    e0, e1 = res["e0"], res["e1"]
    ta_e = None
    for v in ta.values():
        ta_e = (int(v[2]), int(v[3]))
        break
    if ta_e and ta_e != (e0, e1):
        bad.append(f"exponents differ: sage {(e0, e1)} vs Tier-A {ta_e}")
        return 0, bad
    if "deg_lclm_true" in ta:
        done += 1
        if int(ta["deg_lclm_true"][0]) != res["deg_lclm_true"]:
            bad.append(f"deg_lclm_true: sage {res['deg_lclm_true']} vs Tier-A {ta['deg_lclm_true'][0]}")
    if "kerint_true" in ta:
        done += 1
        if int(ta["kerint_true"][0]) != res["kerint_true"]:
            bad.append(f"kerint_true: sage {res['kerint_true']} vs Tier-A {ta['kerint_true'][0]}")
    if "affine_solvable" in ta:
        done += 1
        if bool(int(ta["affine_solvable"][0])) != res["affine_solvable_true"]:
            bad.append(f"affine_solvable: sage {res['affine_solvable_true']} vs "
                       f"Tier-A {ta['affine_solvable'][0]}")
    if "deg_lclm_rand" in ta:
        done += 1
        lo_hi = ta["deg_lclm_rand"][0]
        # the Tier-A row stores the gcrd range over random pairs as "lo..hi"
        try:
            lo, hi = [int(v) for v in lo_hi.split("..")]
        except ValueError:
            lo = hi = int(lo_hi)
        if not (res["deg_gcrd_rand_min"] >= 0 and res["deg_gcrd_rand_max"] <= max(hi, 0)):
            bad.append(f"random gcrd range: sage {res['deg_gcrd_rand_min']}..{res['deg_gcrd_rand_max']}"
                       f" vs Tier-A {lo_hi}")
    return done, bad


def run(cmd, log, budget):
    t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL, start_new_session=True)
        try:
            rc = p.wait(timeout=budget)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(p.pid), signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                p.kill()
            p.wait()
            rc = -9
    return rc, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default=",".join(str(n) for n in range(8, 25)))
    ap.add_argument("--interp-max", type=int, default=16,
                    help="largest n for which the 2^n evaluation-interpolation is attempted")
    ap.add_argument("--interp-budget", type=float, default=3600.0)
    ap.add_argument("--symbolic-max", type=int, default=24)
    ap.add_argument("--symbolic-budget", type=float, default=3600.0)
    ap.add_argument("--trials", type=int, default=32)
    ap.add_argument("--workdir", default=os.environ.get("AIM3_SCRATCH") or "/tmp/t4_ore")
    ap.add_argument("--csv", default="E07_server_ore_scale.csv")
    a = ap.parse_args()
    sage = os.environ.get("AIM3_SAGE") or shutil.which("sage")
    if not sage:
        raise SystemExit("sage not found (set AIM3_SAGE)")
    ns = [int(v) for v in a.ns.split(",")]
    wd = os.path.abspath(a.workdir)
    os.makedirs(wd, exist_ok=True)
    print(f"# params: ns={ns} interp_max={a.interp_max} symbolic_max={a.symbolic_max} "
          f"trials={a.trials} budgets=({a.interp_budget},{a.symbolic_budget})s sage={sage}")
    print("# instances come from the Tier-A probe itself (same seed), so the two implementations "
          "are compared on identical data")

    agree_total, disagree_total = 0, 0
    for n in ns:
        inst = os.path.join(wd, f"inst_n{n}.json")
        rc, _ = run([sys.executable, PROBE, "--export", inst, "--n", str(n), "--no-csv"],
                    os.path.join(wd, f"export_n{n}.log"), 600)
        if rc != 0 or not os.path.exists(inst):
            print(f"n={n}: export failed (rc={rc}) -- skipped")
            continue
        outj = os.path.join(wd, f"ore_n{n}.json")
        cmd = [sage, SAGE_SCRIPT, inst, "--out", outj, "--trials", str(a.trials),
               "--interp-budget", str(a.interp_budget), "--symbolic-budget", str(a.symbolic_budget)]
        if n > a.interp_max:
            cmd.append("--no-interp")
        if n > a.symbolic_max:
            cmd.append("--no-symbolic")
        log = os.path.join(wd, f"ore_n{n}.log")
        budget = a.interp_budget + a.symbolic_budget + 1800
        print(f"\n== n={n}  (interp={'yes' if n <= a.interp_max else 'no'}, "
              f"symbolic={'yes' if n <= a.symbolic_max else 'no'})", flush=True)
        rc, secs = run(cmd, log, budget)
        if rc != 0 or not os.path.exists(outj):
            print(f"   sage failed or timed out (rc={rc}, {secs:.0f}s) -- recorded censored")
            record(a.csv, dict(n=n, measurement="ore_scale", value="", detail=f"rc={rc}",
                               tool="sage-ore"), seconds=secs, censored=1, tier="B")
            continue
        res = json.load(open(outj))
        ta = tier_a_rows(n)
        done, bad = compare(n, res, ta)
        agree_total += done - len(bad)
        disagree_total += len(bad)
        status = "no Tier-A row" if not ta else (f"AGREES on {done - len(bad)}/{done}"
                                                if not bad else "DISAGREES: " + "; ".join(bad))
        print(f"   deg lclm {res['deg_lclm_true']} (sum {res['deg_sum']}), gcrd {res['deg_gcrd_true']}, "
              f"kerint {res['kerint_true']}, solvable {res['affine_solvable_true']}  --  {status}",
              flush=True)
        sym = res.get("symbolic") or {}
        itp = res.get("interp") or {}
        if sym:
            print(f"   symbolic: D_sym(lclm)={sym['D_sym_lclm']} den={sym['D_sym_lclm_den']} "
                  f"D_sym(u)={sym['D_sym_u']} log2 formal deg E={sym['log2_deg_E_formal']:.1f}")
        if itp and not itp.get("censored"):
            print(f"   interpolation: D_interp={itp['D_interp_lclm']} (max {2**n - 1}), "
                  f"ANF={itp['anf_deg_lclm']} (max {n}), drops={itp['degree_drops']}")
        record(a.csv, dict(
            n=n, e0=res["e0"], e1=res["e1"], measurement="ore_scale", tool="sage-ore",
            deg_Lambda=";".join(map(str, res["deg_Lambda"])),
            deg_gcrd_true=res["deg_gcrd_true"], deg_lclm_true=res["deg_lclm_true"],
            deg_sum=res["deg_sum"], kerint_true=res["kerint_true"],
            affine_solvable_true=int(res["affine_solvable_true"]),
            deg_gcrd_rand=f"{res['deg_gcrd_rand_min']}..{res['deg_gcrd_rand_max']}",
            affine_solvable_rand=f"{res['affine_solvable_rand']}/{res['trials']}",
            deg_gcrd_sym=sym.get("deg_gcrd_sym", ""), deg_lclm_sym=sym.get("deg_lclm_sym", ""),
            D_sym_lclm=sym.get("D_sym_lclm", ""), D_sym_lclm_den=sym.get("D_sym_lclm_den", ""),
            D_sym_u=sym.get("D_sym_u", ""), deg_u=sym.get("deg_u", ""), deg_v=sym.get("deg_v", ""),
            cofactors_check=int(sym["cofactors_check"]) if sym else "",
            log2_deg_E_formal=round(sym["log2_deg_E_formal"], 2) if sym else "",
            D_interp_lclm=itp.get("D_interp_lclm", ""), anf_deg_lclm=itp.get("anf_deg_lclm", ""),
            interp_degree_drops=itp.get("degree_drops", ""),
            interp_censored=int(bool(itp.get("censored", 0))) if itp else "",
            tierA_checks=done, tierA_disagreements=len(bad),
            detail=("; ".join(bad) if bad else ("agrees with Tier-A" if ta else "no Tier-A row"))),
            seconds=secs, censored=0, tier="B")

    print(f"\n# cross-validation total: {agree_total} agreements, {disagree_total} disagreements")
    if disagree_total:
        print("# a disagreement is a REAL finding: two independent implementations of the same "
              "algebra differ, and the report must say which one was checked further")


if __name__ == "__main__":
    main()
