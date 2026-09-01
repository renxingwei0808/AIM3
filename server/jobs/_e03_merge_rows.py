#!/usr/bin/env python3
"""Append the server E03 rows that `summarise.py` can read to a merged CSV; print how many were not.

Helper for `server/jobs/e03_summary_merged.sh`.  `summarise.py::load` parses `IV_deg3` and
`outside_ideal` as integers, so a row produced with `--no-points` (empty point-set columns) would
crash it.  Those rows are dropped here — never filled in with a placeholder, which would fabricate
data — and counted on stdout.

    python3 server/jobs/_e03_merge_rows.py SERVER_CSV MERGED_CSV
"""
import csv
import sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    with open(src) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        print(0)
        return
    keep = [r for r in rows if r.get("IV_deg3") and r.get("outside_ideal")]
    with open(dst, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        for r in keep:
            w.writerow(r)
    print(len(rows) - len(keep))


if __name__ == "__main__":
    main()
