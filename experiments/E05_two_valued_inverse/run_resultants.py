#!/usr/bin/env python3
"""E05 task 7 (Tier C) — resultant of the reduced system for n = 8.. , with the T3 slope.

For each n: regenerate the instance with `reduced_system_toy.py --export`, run Magma's
`Resultant(Eq1, Eq2, t1)` through `magma/reduced_system_resultant.m` under a wall-clock and a
memory cap, and record time, memory, the degree and term count of the eliminated polynomial, the
number of its roots and whether the true key is one of them.

Escalation stops at the first n that times out or runs out of memory: that n is recorded with
`censored=1`, its numbers are LOWER bounds, and no larger n is attempted (batch rule 0.2 and
CLAUDE.md ground rule 4).

Finally `log2(time)` and `log2(degree)` are regressed against n by least squares.  ROADMAP §3.3
closes track T3 when the log2(time) slope is >= 2 and nothing factors.

Usage:
  python3 run_resultants.py [--ns 8,9,...,20] [--seed 1] [--timeout 7200] [--mem-gb 120]
                            [--workdir DIR] [--csv E05_server_resultants.csv] [--factor]
"""
import argparse
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from _common import record                                        # noqa: E402

TOY = os.path.join(HERE, "reduced_system_toy.py")
DRIVER = os.path.join(HERE, "magma", "reduced_system_resultant.m")


def least_squares(xs, ys):
    """slope, intercept of the least-squares line through (xs, ys); None if fewer than 2 points."""
    if len(xs) < 2:
        return None, None
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return None, None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
    return slope, my - slope * mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default=",".join(str(n) for n in range(8, 21)))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=7200.0)
    ap.add_argument("--mem-gb", type=int, default=120)
    ap.add_argument("--workdir", default=os.environ.get("AIM3_SCRATCH") or "/tmp/aim3_e05")
    ap.add_argument("--csv", default="E05_server_resultants.csv")
    ap.add_argument("--factor", action="store_true", help="also try Factorization (T3.3 lead)")
    a = ap.parse_args()
    ns = [int(v) for v in a.ns.split(",")]
    os.makedirs(a.workdir, exist_ok=True)
    magma = os.environ.get("AIM3_MAGMA") or shutil.which("magma")
    if not magma:
        raise SystemExit("magma not found (set AIM3_MAGMA)")
    print(f"# params: ns={ns} seed={a.seed} timeout={a.timeout}s mem={a.mem_gb}GB "
          f"workdir={a.workdir} factor={a.factor} magma={magma}")

    fitted = []
    for n in ns:
        print(f"\n== n={n}", flush=True)
        sysfile = os.path.join(a.workdir, f"E05_n{n}_s{a.seed}.magma.m")
        if not os.path.exists(sysfile):
            r = subprocess.run([sys.executable, TOY, str(n), str(a.seed), "--no-roots",
                                "--export", a.workdir],
                               capture_output=True, text=True)
            for line in r.stdout.splitlines():
                if line.startswith(("n=", "   Eq", "   Bez", "   [export]")):
                    print("  " + line.strip(), flush=True)
            if r.returncode != 0 or not os.path.exists(sysfile):
                print(f"  !! export failed: {r.stderr[-800:]}")
                break

        # Shape of the system, read from the exported file rather than from Magma's output: a
        # censored run writes no JSON, and a row that says only "censored" tells the next reader
        # nothing about what was attempted.
        shape = {}
        try:
            txt = open(sysfile).read()
            for name in ("EQ1", "EQ2"):
                body = txt.split(name + " := [", 1)[1].split("];", 1)[0]
                trips = [t for t in body.replace("\n", "").split("<") if t.strip()]
                degs = [(int(t.split(",")[0]), int(t.split(",")[1])) for t in trips]
                shape[name.lower() + "_terms"] = len(degs)
                shape[name.lower() + "_deg_t0"] = max((d[0] for d in degs), default=0)
                shape[name.lower() + "_deg_t1"] = max((d[1] for d in degs), default=0)
            for key in ("exps", "eps", "key"):
                for line in txt.splitlines():
                    if line.startswith(key + " :="):
                        shape[key] = line.split(":=", 1)[1].split(";")[0].strip()
        except Exception as exc:
            print(f"  (could not read the shape from {sysfile}: {exc})")

        wrapper = os.path.join(a.workdir, f"E05_n{n}_s{a.seed}.run.m")
        jsonout = os.path.join(a.workdir, f"E05_n{n}_s{a.seed}.res.json")
        log = os.path.join(a.workdir, f"E05_n{n}_s{a.seed}.magma.log")
        with open(wrapper, "w") as fh:
            fh.write(f'load "{sysfile}";\n')
            fh.write(f'out := "{jsonout}";\n')
            fh.write(f"dofactor := {1 if a.factor else 0};\n")
            fh.write("doroots := 1;\n")
            fh.write(f'load "{DRIVER}";\n')
        if os.path.exists(jsonout):
            os.remove(jsonout)

        # ulimit -v is applied to the Magma process only (the caps of server/README.md)
        cmd = ["bash", "-c",
               f"ulimit -v {a.mem_gb * 1024 * 1024}; exec /usr/bin/time -v -o {log}.time "
               f"{magma} -b {wrapper}"]
        t0 = time.time()
        censored = 0
        # own process group, killed as a group: otherwise the timeout only reaps the wrapper and
        # Magma keeps running (and holding its memory) long after we have moved on
        with open(log, "w") as fh:
            proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, start_new_session=True)
            try:
                rc = proc.wait(timeout=a.timeout)
            except subprocess.TimeoutExpired:
                censored = 1
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    proc.kill()
                proc.wait()
                rc = -9
        wall = time.time() - t0
        peak_kb = None
        if os.path.exists(log + ".time"):
            m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(log + ".time").read())
            if m:
                peak_kb = int(m.group(1))
        res = {}
        if os.path.exists(jsonout):
            try:
                res = json.load(open(jsonout))
            except Exception:
                res = {}
        if not res:
            censored = 1
        for line in open(log, errors="replace").read().splitlines():
            if line.startswith("###"):
                print("  " + line, flush=True)
        print(f"  -> rc={rc} wall={wall:.1f}s peakRSS={(peak_kb or 0)/2**20:.1f}GB "
              f"{'[CENSORED: lower bound]' if censored else ''}", flush=True)

        record(a.csv, dict(
            n=n, seed=a.seed,
            exps=res.get("exps") or shape.get("exps", ""),
            eps=res.get("eps") or shape.get("eps", ""),
            eq1_terms=res.get("eq1_terms") or shape.get("eq1_terms", ""),
            eq2_terms=res.get("eq2_terms") or shape.get("eq2_terms", ""),
            eq1_deg_t0=shape.get("eq1_deg_t0", ""), eq2_deg_t0=shape.get("eq2_deg_t0", ""),
            eq1_deg_t1=res.get("eq1_deg_t1") or shape.get("eq1_deg_t1", ""),
            eq2_deg_t1=res.get("eq2_deg_t1") or shape.get("eq2_deg_t1", ""),
            key_is_root_of_system=res.get("key_is_root_of_system", ""),
            resultant_seconds=res.get("resultant_seconds", ""),
            resultant_degree=res.get("resultant_degree", ""),
            resultant_terms=res.get("resultant_terms", ""),
            resultant_is_zero=res.get("resultant_is_zero", ""),
            n_roots=res.get("n_roots", ""), root_seconds=res.get("root_seconds", ""),
            key_is_root_of_resultant=res.get("key_is_root_of_resultant", ""),
            n_factors=res.get("n_factors", ""), factor_degrees=res.get("factor_degrees", ""),
            peak_rss_kb=peak_kb or "", magma_peak_bytes=res.get("peak_bytes", ""),
            wall_seconds=round(wall, 3), rc=rc, tool="magma",
            tool_version=res.get("tool_version", "")),
            seconds=wall, censored=censored, tier="C")

        if censored:
            print(f"  !! n={n} hit the {a.timeout:.0f}s / {a.mem_gb} GB cap -- reporting a lower "
                  f"bound and stopping the escalation (batch rule 0.2)", flush=True)
            break
        if res.get("resultant_seconds") is not None and res.get("resultant_degree", -1) > 0:
            fitted.append((n, float(res["resultant_seconds"]), int(res["resultant_degree"])))

    print("\n# fit (COMPLEXITY_ACCOUNTING section 4: >= 4 values of n before any extrapolation)")
    if len(fitted) < 2:
        print("  too few completed sizes to fit anything")
        return
    xs = [n for n, _, _ in fitted]
    st, _ = least_squares(xs, [math.log2(max(t, 1e-6)) for _, t, _ in fitted])
    sd, _ = least_squares(xs, [math.log2(d) for _, _, d in fitted])
    print(f"  completed n = {xs}")
    print(f"  slope of log2(resultant seconds) vs n = {st:.3f}")
    print(f"  slope of log2(resultant degree)  vs n = {sd:.3f}")
    print(f"  ROADMAP section 3.3 termination criterion (time slope >= 2): "
          f"{'MET' if st is not None and st >= 2 else 'NOT met'}"
          + ("" if len(xs) >= 4 else "   [fewer than 4 sizes -- not an extrapolation]"))


if __name__ == "__main__":
    main()
