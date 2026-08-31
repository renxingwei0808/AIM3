#!/usr/bin/env python3
"""E04 task 6 — run F4 on the exported systems and record the real solving degree `d_solv`.

Saturation destroys the D=4 fall count of the full system (REPORT_tierA_batch1 §5.3); the F4 step
degrees do not.  This driver runs each exported system through

  * **msolve** `-v 2 -t <threads>`  — its verbose log prints one line per matrix, "<degree> ... R x C"
  * **Magma** `SetVerbose("Faugere", 2)` via `magma/aim3_f4_degrees.m`

and records, per (system, order, tool): `d_solv` = the largest step degree actually reached,
`max_rows` = the largest F4 matrix, wall-clock seconds and peak RSS (from `/usr/bin/time -v` when
available, else from the tool's own report).  A run that hits `--timeout` is written with
`censored=1` and its `d_solv` is a LOWER bound, never a value (CLAUDE.md ground rule 4).

At n = 7 the two tools are cross-validated against each other (`--cross-check`): a disagreement in
`d_solv` is reported, never averaged away.

Usage:
  python3 run_f4_degrees.py --systems DIR [--tags SUBSTR] [--tools msolve,magma]
                            [--orders grevlex] [--threads 32] [--timeout 7200]
                            [--csv E04_f4_degrees.csv] [--logdir DIR]
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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record                                        # noqa: E402
from aim3.polysys_export import parse_msolve_log, parse_magma_log  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MAGMA_DRIVER = os.path.join(HERE, "magma", "aim3_f4_degrees.m")


def _time_v(cmd, log_path, timeout):
    """Run `cmd`, tee stdout+stderr to `log_path`, return (rc, seconds, peak_rss_kb, timed_out).

    The child runs in its own process group and the whole group is killed on timeout: killing only
    the direct child would leave `/usr/bin/time`'s grandchild (magma, msolve) running for hours
    after the driver has moved on.
    """
    tf = log_path + ".time"
    use_time = shutil.which("/usr/bin/time") is not None
    full = (["/usr/bin/time", "-v", "-o", tf] if use_time else []) + cmd
    t0 = time.time()
    timed_out = False
    with open(log_path, "w") as fh:
        proc = subprocess.Popen(full, stdout=fh, stderr=subprocess.STDOUT,
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
    peak = None
    if use_time and os.path.exists(tf):
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(tf).read())
        if m:
            peak = int(m.group(1))
    return rc, secs, peak, timed_out


def run_msolve(info, sysdir, logdir, threads, timeout, gb_only=True):
    ms = os.path.join(sysdir, info["tag"] + ".ms")
    log = os.path.join(logdir, info["tag"] + ".msolve.log")
    out = os.path.join(logdir, info["tag"] + ".msolve.out")
    exe = os.environ.get("AIM3_MSOLVE") or shutil.which("msolve")
    if not exe:
        return None
    # -g 2 stops after the reduced Groebner basis.  Without it msolve goes on to a rational
    # parametrization, which re-runs F4 ~20 times with different random choices: the step degrees
    # are identical each time, but the step count and the wall clock would not describe one F4 run.
    cmd = [exe, "-v", "2", "-t", str(threads), "-f", ms, "-o", out]
    if gb_only:
        cmd += ["-g", "2"]
    rc, secs, peak, to = _time_v(cmd, log, timeout)
    parsed = parse_msolve_log(open(log, errors="replace").read())
    return dict(tool="msolve", order="grevlex(drl)", rc=rc, seconds=secs, peak_rss_kb=peak,
                censored=int(to or rc != 0), log=log, **parsed)


def run_magma(info, sysdir, logdir, order, threads, timeout):
    sysfile = os.path.join(sysdir, info["tag"] + ".magma.m")
    safe = order.replace(":", "_")
    wrapper = os.path.join(logdir, f"{info['tag']}.{safe}.run.m")
    jsonout = os.path.join(logdir, f"{info['tag']}.{safe}.f4.json")
    log = os.path.join(logdir, f"{info['tag']}.{safe}.magma.log")
    exe = os.environ.get("AIM3_MAGMA") or shutil.which("magma")
    if not exe:
        return None
    # Magma's `load` takes a string literal only, hence the generated wrapper.
    with open(wrapper, "w") as fh:
        fh.write(f'load "{os.path.abspath(sysfile)}";\n')
        fh.write(f'ord := "{order}";\n')
        fh.write(f'out := "{os.path.abspath(jsonout)}";\n')
        fh.write(f'nthreads := "{threads}";\n')
        fh.write(f'load "{MAGMA_DRIVER}";\n')
    rc, secs, peak, to = _time_v([exe, "-b", wrapper], log, timeout)
    parsed = parse_magma_log(open(log, errors="replace").read())
    extra = {}
    if os.path.exists(jsonout):
        try:
            extra = json.load(open(jsonout))
        except Exception:
            extra = {}
    return dict(tool="magma", order=order, rc=rc, seconds=secs, peak_rss_kb=peak,
                censored=int(to or rc != 0), log=log, gb_size=extra.get("gb_size"),
                gb_max_degree=extra.get("gb_max_degree"), magma_seconds=extra.get("seconds"),
                ring=extra.get("ring"),
                **parsed)


def resolve_order(order, info):
    """Turn a symbolic block order into the concrete `elim:` string for this system, or None.

    Variable layout: `S_quad`/`P1only`/`random` use pt = 1..n then t_0..t_{l-1}; the E06 kernel
    model uses U = 1..3n then W = 3n+1..6n.  A symbolic order that does not fit the system's shape
    returns None and the run is skipped rather than silently mislabelled.
    """
    nv = info["nvars"]
    n = info.get("n")
    exp = info.get("experiment")
    if order == "genref":
        return "gen:grevlex"
    if order in ("blockPT", "blockT"):
        if exp != "E04" or not n or nv <= n:
            return None
        return f"elim:{n}" if order == "blockPT" else f"elimlast:{nv - n}"
    if order in ("blockU", "blockW"):
        if exp != "E06" or nv % 2:
            return None
        return f"elim:{nv // 2}" if order == "blockU" else f"elimlast:{nv // 2}"
    return order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", required=True, help="directory written by export_systems.py")
    ap.add_argument("--tags", default="", help="only run systems whose tag contains this substring")
    ap.add_argument("--tools", default="msolve,magma")
    ap.add_argument("--orders", default="grevlex",
                    help="Magma monomial orders, comma separated.  Besides the literal orders "
                         "(grevlex, lex, elim:K, elimlast:K) four symbolic block orders are "
                         "resolved per system (E10): genref is grevlex on the same general ring "
                         "the block orders need (the only fair baseline for them); blockPT / "
                         "blockT put the pt block resp. the "
                         "t block first in an S_quad-shaped system, blockU / blockW put U resp. W "
                         "first in the E06 kernel model.  An order that does not apply to a system "
                         "is skipped.")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--timeout", type=float, default=7200.0)
    ap.add_argument("--csv", default="E04_f4_degrees.csv")
    ap.add_argument("--logdir", default=None)
    ap.add_argument("--msolve-full-solve", action="store_true",
                    help="let msolve continue to the rational parametrization (default: -g 2, "
                         "reduced Groebner basis only, so the timing describes one F4 run)")
    a = ap.parse_args()

    sysdir = os.path.abspath(a.systems)
    logdir = os.path.abspath(a.logdir or os.path.join(sysdir, "..", "f4logs"))
    os.makedirs(logdir, exist_ok=True)
    tools = a.tools.split(",")
    orders = a.orders.split(",")
    infos = []
    for p in sorted(glob.glob(os.path.join(sysdir, "*.json"))):
        if os.path.basename(p) == "index.json":
            continue
        info = json.load(open(p))
        if a.tags and a.tags not in info["tag"]:
            continue
        infos.append(info)
    print(f"# params: systems={sysdir} n={len(infos)} tools={tools} orders={orders} "
          f"threads={a.threads} timeout={a.timeout}s logdir={logdir}")

    for info in infos:
        print(f"\n== {info['tag']}  ({info['neqs']} eqs / {info['nvars']} vars)", flush=True)
        runs = []
        if "msolve" in tools:
            r = run_msolve(info, sysdir, logdir, a.threads, a.timeout,
                           gb_only=not a.msolve_full_solve)
            if r:
                runs.append(r)
        if "magma" in tools:
            for order in orders:
                concrete = resolve_order(order, info)
                if concrete is None:
                    continue
                r = run_magma(info, sysdir, logdir, concrete, a.threads, a.timeout)
                if r:
                    r["order_label"] = order
                    runs.append(r)
        for r in runs:
            print(f"   {r['tool']:7s} order={r['order']:12s} d_solv={r['d_solv']} "
                  f"steps={r['nsteps']} max_rows={r['max_rows']} {r['seconds']:.1f}s "
                  f"peakRSS={(r['peak_rss_kb'] or 0)/1024:.0f}MB"
                  + ("  [CENSORED: lower bound]" if r["censored"] else ""), flush=True)
            record(a.csv, dict(
                n=info.get("n", ""), ell=info.get("ell", ""), exps=info.get("exps", ""),
                seed=info.get("seed", ""), system=info.get("system", info["tag"]),
                experiment=info.get("experiment", ""), tag=info["tag"],
                nvars=info["nvars"], neqs=info["neqs"], order=r["order"],
                order_label=r.get("order_label", r["order"]), tool=r["tool"],
                d_solv=r["d_solv"] if r["d_solv"] is not None else "",
                nsteps=r["nsteps"], max_rows=r["max_rows"] or "", max_cols=r["max_cols"] or "",
                gb_size=r.get("gb_size") or "", gb_max_degree=r.get("gb_max_degree") or "",
                ring=r.get("ring") or "",
                peak_rss_kb=r["peak_rss_kb"] or "", rc=r["rc"],
                log=os.path.relpath(r["log"], logdir)),
                seconds=r["seconds"], censored=r["censored"], tier="B")
        # cross-check when both tools ran on the same system
        ds = {r["tool"]: r["d_solv"] for r in runs if not r["censored"] and r["d_solv"] is not None}
        if len(ds) > 1 and len(set(ds.values())) > 1:
            print(f"   !! d_solv DISAGREEMENT between tools: {ds}", flush=True)
    print("\n# done")


if __name__ == "__main__":
    main()
