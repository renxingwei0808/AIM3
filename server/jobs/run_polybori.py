#!/usr/bin/env python3
"""Batch-2 task 3.3: drive `sage/polybori_block_orders.sage` and record its step degrees.

PolyBoRi's `prot=True` prints one `Current Degree: d` line per degree it works at, computing
**directly in the declared order** — which is the block-order solving degree Magma 2.24 could not
supply (server batch 1, E10).  This driver runs one system under several orders, parses those
lines, and writes one CSV row per (system, order) with

    d_solv_polybori = max over the printed degrees
    gb_size, gb_max_degree, seconds, peak RSS

so grevlex, U > W and W > U can be compared on the same instance with the same tool.

Usage:
  python3 server/jobs/run_polybori.py --systems DIR --tags SUBSTR [--orders grevlex,blockU,blockW]
                                      [--timeout 3600] [--csv E10_server_polybori.csv]

Symbolic orders, resolved per system from its `.json` descriptor:
  `blockU`  = the first half of the variables first  (kernel model: U = pt,t0,t1)
  `blockW`  = the last  half first                    (kernel model: W = a0,a1,a2)
  `blockPT` = the first n variables first             (S_quad-shaped: the pt block)
  `blockT`  = the last nvars−n first                  (S_quad-shaped: the t blocks)
"""
import argparse
import glob
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from _common import record                                          # noqa: E402

SAGE_SCRIPT = os.path.join(ROOT, "experiments", "E10_block_order_groebner", "sage",
                           "polybori_block_orders.sage")
STEP = re.compile(r"Current Degree:\s*(\d+)")


def resolve(order, info):
    nv, n, exp = info["nvars"], info.get("n"), info.get("experiment", "")
    if order in ("blockU", "blockW"):
        if nv % 2:
            return None
        return f"block:{nv//2}" if order == "blockU" else f"blocklast:{nv//2}"
    if order in ("blockPT", "blockT"):
        if not n or nv <= n:
            return None
        return f"block:{n}" if order == "blockPT" else f"blocklast:{nv-n}"
    return order


def run_one(sysfile, order, logpath, jsonpath, sage, timeout):
    cmd = [sage, SAGE_SCRIPT, sysfile, "--order", order, "--out", jsonpath]
    tf = logpath + ".time"
    if shutil.which("/usr/bin/time"):
        cmd = ["/usr/bin/time", "-v", "-o", tf] + cmd
    t0 = time.time()
    timed_out = False
    with open(logpath, "w") as fh:
        proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, start_new_session=True)
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                proc.kill()
            proc.wait()
            rc = -9
    secs = time.time() - t0
    text = open(logpath, errors="replace").read()
    degs = [int(m) for m in STEP.findall(text)]
    peak = None
    if os.path.exists(tf):
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(tf).read())
        if m:
            peak = int(m.group(1))
    res = {}
    if os.path.exists(jsonpath):
        try:
            res = json.load(open(jsonpath))
        except Exception:
            res = {}
    return dict(rc=rc, seconds=secs, peak_rss_kb=peak, censored=int(timed_out or not res),
                d_solv=max(degs) if degs else None, nsteps=len(degs), steps=degs, **res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", required=True)
    ap.add_argument("--tags", default="")
    ap.add_argument("--orders", default="grevlex,blockU,blockW")
    ap.add_argument("--timeout", type=float, default=3600.0)
    ap.add_argument("--csv", default="E10_server_polybori.csv")
    ap.add_argument("--logdir", default=None)
    a = ap.parse_args()
    sage = os.environ.get("AIM3_SAGE") or shutil.which("sage")
    if not sage:
        raise SystemExit("sage not found (set AIM3_SAGE)")
    sysdir = os.path.abspath(a.systems)
    logdir = os.path.abspath(a.logdir or os.path.join(sysdir, "..", "polybori_logs"))
    os.makedirs(logdir, exist_ok=True)
    orders = a.orders.split(",")
    infos = []
    for p in sorted(glob.glob(os.path.join(sysdir, "*.json"))):
        if os.path.basename(p) == "index.json":
            continue
        info = json.load(open(p))
        if a.tags and a.tags not in info["tag"]:
            continue
        infos.append(info)
    print(f"# params: systems={sysdir} n={len(infos)} orders={orders} timeout={a.timeout}s "
          f"sage={sage}")

    for info in infos:
        print(f"\n== {info['tag']}  ({info['neqs']} eqs / {info['nvars']} vars)", flush=True)
        for label in orders:
            concrete = resolve(label, info)
            if concrete is None:
                continue
            sysfile = os.path.join(sysdir, info["tag"] + ".magma.m")
            base = os.path.join(logdir, f"{info['tag']}.{label}")
            r = run_one(sysfile, concrete, base + ".log", base + ".json", sage, a.timeout)
            print(f"   {label:<8} ({concrete:<14}) d_solv={r['d_solv']} steps={r['nsteps']} "
                  f"gb={r.get('gb_size')}/{r.get('gb_max_degree')} {r['seconds']:.1f}s "
                  f"peak={(r['peak_rss_kb'] or 0)/1024:.0f}MB"
                  + ("  [CENSORED: lower bound]" if r["censored"] else ""), flush=True)
            record(a.csv, dict(
                n=info.get("n", ""), ell=info.get("ell", ""), seed=info.get("seed", ""),
                exps=info.get("exps", ""), system=info.get("system", ""),
                experiment=info.get("experiment", ""), tag=info["tag"],
                nvars=info["nvars"], neqs=info["neqs"],
                order_label=label, order=concrete, ring_order=r.get("ring_order", ""),
                tool="polybori",
                d_solv=r["d_solv"] if r["d_solv"] is not None else "",
                nsteps=r["nsteps"], steps=";".join(map(str, r["steps"])),
                gb_size=r.get("gb_size", ""), gb_max_degree=r.get("gb_max_degree", ""),
                trivial_gb=r.get("trivial_gb", ""),
                peak_rss_kb=r["peak_rss_kb"] or "", rc=r["rc"]),
                seconds=r["seconds"], censored=r["censored"], tier="B")
    print("\n# done")


if __name__ == "__main__":
    main()
