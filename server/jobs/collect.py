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

# per-job CSV name  ->  (glob of job directories, output file in results/)
MERGE = [
    ("E03_sbox_degree_fall.csv", "E03_sweep_n*", "E03_server_n17_20.csv"),
    ("E03_server_D5.csv", "E03_D5_n*", "E03_server_D5.csv"),
    ("E04_full_system_falls.csv", "E04_gate1_l*", "E04_server_gate1.csv"),
    ("E04_server_gate1_pairs.csv", "E04_gate1pairs_*", "E04_server_gate1_pairs.csv"),
    ("E04_f4_degrees.csv", "E04_f4_degrees_*", "E04_f4_degrees.csv"),
    ("E10_block_orders.csv", "E10_blockorder_*", "E10_server_block_orders.csv"),
    ("E12_first_layer_ideal.csv", "E12_first_layer_n*", "E12_server_first_layer.csv"),
    ("E12_cross_relations.csv", "E12_cross_rel_*", "E12_server_cross_relations_deg4.csv"),
    ("E12_cross_relations_explained.csv", "E12_explain_*",
     "E12_server_cross_relations_explained_deg4.csv"),
    ("E12_server_criteria_subtract.csv", "E12_criteria_*", "E12_server_criteria_subtract.csv"),
    ("E06_kernel_model_falls.csv", "E06_kernel_n*", "E06_server_kernel.csv"),
    ("E05_server_resultants.csv", "E05_resultants*", "E05_server_resultants.csv"),
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

    # ---- 2. per-job CSVs, merged
    for csvname, jobglob, outname in MERGE:
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
