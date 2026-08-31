#!/usr/bin/env python3
"""Pull the server batch back into the repository (run on the LOCAL machine).

Two kinds of artefact come back, and only these two (batch rule 0.4 — matrices and intermediate
bases stay on the server):

  * `results/server/<date>_<EXP>/*.{log,time,meta}` — the job logs and their accounting stamps;
    `.log` is git-ignored by the repository's own rules, `.meta` and `.time` are not.
  * the per-job CSVs written under `$AIM3_BASE/out/<job>/`, merged into the batch's new tables in
    `results/` (one output file per experiment, rows in job order, a single header).

Every server job writes its CSVs into its own `$AIM3_RESULTS` directory *outside* the checkout, so
a run can never rewrite a committed Tier-A table; this script is the only thing that puts server
numbers into `results/`.

Usage:
  python3 server/jobs/collect.py [--dry-run] [--host server] [--base /home/luck/rxw/AIM3_remote]
"""
import argparse
import csv
import io
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# per-job CSV name -> (glob of job directories, output file in results/, optional dedupe key).
# The dedupe key keeps ONE row per key, preferring a complete row (no empty cells) over a partial
# one: the E03 sweeps were run twice at n = 19, 20 -- once with the O(2^n) point matrix and once
# with `--no-points` to finish the fall-count table -- and only the (n, e) pair identifies a run.
MERGE = [
    ("E03_sbox_degree_fall.csv", "E03_sweep_n*", "E03_server_n17_20.csv", ("n", "e")),
    ("E03_server_D5.csv", "E03_D5_n*", "E03_server_D5.csv"),
    ("E04_full_system_falls.csv", "E04_gate1_l*", "E04_server_gate1.csv"),
    ("E04_server_gate1_pairs.csv", "E04_gate1pairs_*", "E04_server_gate1_pairs.csv"),  # noqa
    # The F4 jobs were launched with the block orders in their order list, before the Magma driver
    # could build a block-order ring at all (BooleanPolynomialRing refuses one -- see the E10
    # README).  Those rows are duplicates of E10_server_block_orders.csv, where they were measured
    # with the corrected driver, so only the grevlex rows are kept here and E04_f4_degrees.csv
    # stays what task 6 asks for: the d_solv table.
    ("E04_f4_degrees.csv", "E04_f4_degrees_*", "E04_f4_degrees.csv", None,
     lambda row, idx: row[idx["order"]].startswith("grevlex")),
    ("E10_block_orders.csv", "E10_blockorder_*", "E10_server_block_orders.csv"),
    ("E12_first_layer_ideal.csv", "E12_first_layer_n*", "E12_server_first_layer.csv"),
    ("E12_cross_relations.csv", "E12_cross_rel_*", "E12_server_cross_relations_deg4.csv"),
    ("E12_cross_relations_explained.csv", "E12_explain_*",
     "E12_server_cross_relations_explained_deg4.csv"),
    ("E12_server_criteria_subtract.csv", "E12_criteria_*", "E12_server_criteria_subtract.csv"),
    ("E06_kernel_model_falls.csv", "E06_kernel_n*", "E06_server_kernel.csv"),
    ("E05_server_resultants.csv", "E05_resultants*", "E05_server_resultants.csv", ("n",)),
]


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def remote_cat(host, ssh_config, path):
    r = run(["ssh", "-F", ssh_config, host, f"cat {path} 2>/dev/null"])
    return r.stdout


def remote_ls(host, ssh_config, pattern):
    r = run(["ssh", "-F", ssh_config, host, f"ls -d {pattern} 2>/dev/null"])
    return [p for p in r.stdout.split() if p]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="server")
    ap.add_argument("--ssh-config", default=os.path.expanduser("~/.ssh/config_server"))
    ap.add_argument("--base", default="/home/luck/rxw/AIM3_remote")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    # ---- 1. logs and accounting stamps
    src = f"{a.host}:{a.base}/AIM3/results/server/"
    dst = os.path.join(ROOT, "results", "server") + "/"
    os.makedirs(dst, exist_ok=True)
    cmd = ["rsync", "-az", "-e", f"ssh -F {a.ssh_config}", src, dst]
    print("+ " + " ".join(cmd))
    if not a.dry_run:
        r = run(cmd)
        if r.returncode:
            print(r.stderr[-2000:])
            return 1

    # ---- 1b. the backend-alignment record (the gate for --backend sage)
    for name in ("backend_align.json", "backend_align.log"):
        cmd = ["rsync", "-az", "-e", f"ssh -F {a.ssh_config}",
               f"{a.host}:{a.base}/out/{name}", os.path.join(ROOT, "results", "server", name)]
        print("+ " + " ".join(cmd))
        if not a.dry_run:
            run(cmd)

    # ---- 2. per-job CSVs, merged
    for entry in MERGE:
        csvname, jobglob, outname = entry[:3]
        dedupe = entry[3] if len(entry) > 3 else None
        keep_row = entry[4] if len(entry) > 4 else None
        dirs = remote_ls(a.host, a.ssh_config, f"{a.base}/out/{jobglob}/")
        rows, header = [], None
        used = []
        for d in sorted(dirs):
            text = remote_cat(a.host, a.ssh_config, os.path.join(d, csvname))
            if not text.strip():
                continue
            rdr = csv.reader(io.StringIO(text))
            h = next(rdr, None)
            if h is None:
                continue
            if header is None:
                header = h
            elif h != header:
                print(f"  !! {outname}: header mismatch in {d}\n     {h}\n     {header}")
                continue
            n0 = len(rows)
            rows.extend(r for r in rdr if r)
            used.append(f"{os.path.basename(d.rstrip('/'))}(+{len(rows)-n0})")
        if keep_row and rows:
            idx = {c: k for k, c in enumerate(header)}
            before = len(rows)
            rows = [r for r in rows if keep_row(r, idx)]
            if before != len(rows):
                print(f"     ({before - len(rows)} rows filtered out, see MERGE in collect.py)")
        if dedupe and rows:
            idx = {c: k for k, c in enumerate(header)}
            best = {}
            for r in rows:
                key = tuple(r[idx[c]] for c in dedupe)
                score = sum(1 for v in r if v != "")
                if key not in best or score > best[key][0]:
                    best[key] = (score, r)
            dropped = len(rows) - len(best)
            rows = [v[1] for v in best.values()]
            if dropped:
                print(f"     ({dropped} duplicate rows dropped on {dedupe}, "
                      f"the most complete kept)")
        if not rows:
            print(f"  -- {outname}: nothing yet")
            continue
        out = os.path.join(ROOT, "results", outname)
        print(f"  -> {outname}: {len(rows)} rows from {len(used)} jobs [{', '.join(used)}]")
        if not a.dry_run:
            with open(out, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(header)
                w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
