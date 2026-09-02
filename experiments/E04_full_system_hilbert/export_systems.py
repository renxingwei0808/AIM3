#!/usr/bin/env python3
"""E04 task 6 — export toy AIM3 systems for the F4 step-degree runs (msolve and Magma).

Generates, for every (n, ell, seed) in the grid, the three systems that `full_system_falls.py`
compares — `S_quad` (the official 3n(ell+1) quadratics), `P1only` (the n(ell+1) bilinear P1
coordinates) and a same-shape `random` control with a planted root — from exactly the same
`legacy/aim3_macaulay.build_aim3` instance, and writes each of them as

    <tag>.ms        msolve input (characteristic 2, WITH the field equations x_i^2 + x_i)
    <tag>.magma.m   nvars / sysname / MONS for experiments/.../magma/aim3_f4_degrees.m
    <tag>.json      shape, provenance and the planted solution

with `<tag> = E04_n<n>_l<ell>_s<seed>_<system>`.  The kernel model (E06) and the first-layer system
(E12) are exported by their own scripts' `--export` flag, into the same directory, so that
`run_f4_degrees.py` can sweep all of them together.

Usage:
  python3 export_systems.py --out DIR [--ns 7,9,11,13,15] [--ells 1,2] [--seeds 1,2]
                            [--systems S_quad,P1only,random]
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import COMMIT, HOST                                  # noqa: E402
from aim3_macaulay import build_aim3, random_system               # noqa: E402
from aim3.polysys_export import write_system                      # noqa: E402
from full_system_falls import default_exps, generic_exps          # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="output directory (results/server/<date>_E04/systems)")
    ap.add_argument("--ns", default="7,9,11,13,15")
    ap.add_argument("--ells", default="1,2")
    ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--systems", default="S_quad,P1only,random")
    ap.add_argument("--exps", default=None,
                    help="explicit comma-separated exponents (one list, used for every n); needed "
                         "at the n where fewer than ell+1 non-degenerate exponents exist")
    ap.add_argument("--tag-prefix", default=None, help="override the tag prefix (default E04/E04g)")
    ap.add_argument("--generic-exps", action="store_true",
                    help="use the smallest ell+1 NON-DEGENERATE exponents (batch 2 default)")
    a = ap.parse_args()
    ns = [int(v) for v in a.ns.split(",")]
    ells = [int(v) for v in a.ells.split(",")]
    seeds = [int(v) for v in a.seeds.split(",")]
    systems = a.systems.split(",")
    os.makedirs(a.out, exist_ok=True)
    print(f"# params: ns={ns} ells={ells} seeds={seeds} systems={systems} out={a.out} "
          f"generic_exps={a.generic_exps}")
    index = []
    for n in ns:
        for ell in ells:
            try:
                if a.exps:
                    exps = [int(v) for v in a.exps.split(",")]
                    if len(exps) != ell + 1:
                        raise SystemExit(f"--exps needs {ell+1} values for ell={ell}, got {exps}")
                else:
                    exps = generic_exps(n, ell) if a.generic_exps else default_exps(n, ell)
            except SystemExit as exc:
                print(f"  skip n={n} ell={ell}: {exc}")
                continue
            for seed in seeds:
                nv, eqsP1, eqsAll, sol = build_aim3(n, exps, seed=seed)
                pool = {"S_quad": eqsAll, "P1only": eqsP1}
                if "random" in systems:
                    req, rsol = random_system(nv, len(eqsAll), seed=1000 + seed)
                    pool["random"] = req
                for name in systems:
                    pfx = a.tag_prefix or ("E04g" if (a.generic_exps or a.exps) else "E04")
                    tag = f"{pfx}_n{n}_l{ell}_s{seed}_{name}"
                    info = write_system(pool[name], nv, a.out, tag, sysname=name,
                                        sol=sol if name != "random" else None,
                                        meta=dict(experiment="E04", n=n, ell=ell, seed=seed,
                                                  exps=";".join(map(str, exps)), system=name,
                                                  builder="legacy/aim3_macaulay.build_aim3",
                                                  host=HOST, commit=COMMIT))
                    index.append(info)
                    print(f"  {tag}: {info['neqs']} eqs / {info['nvars']} vars, "
                          f"{info['terms']} terms, max deg {info['max_degree']}")
    with open(os.path.join(a.out, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1, sort_keys=True)
    print(f"# {len(index)} systems written; index.json lists them")


if __name__ == "__main__":
    main()
