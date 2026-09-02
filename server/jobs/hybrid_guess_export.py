#!/usr/bin/env python3
"""Batch-2 task 1.3: export the systems obtained by guessing the first k bits of pt correctly.

The hybrid-guessing cost of `COMPLEXITY_ACCOUNTING §6` is `min_k 2^k · (cost of solving after k bits
are fixed)`, and every published estimate takes the second factor from a semi-regular model.  This
script produces the systems that measure it directly: it builds a toy AIM3 instance, substitutes the
**true** values of pt bits 0..k−1 (so the planted root survives and the run is not a search), renumbers
the remaining 3n−k variables, and writes the result in msolve and Magma syntax.

`run_f4_degrees.py` then measures `d_solv` on them, and

    δ(n, k) = dreg_boolean(3n − k, 9n) − d_solv

is how far the B2 model over-estimates the solving degree once k bits are guessed.  δ is what
`calibrated_estimate.py` turns into a corrected complexity.

Guessing the *true* bits is deliberate: it isolates the solving degree of the correct branch.  The
wrong-branch cost is a separate question (the system becomes inconsistent and F4 finds 1 quickly),
and `COMPLEXITY_ACCOUNTING §6` counts it separately.

Usage:
  python3 server/jobs/hybrid_guess_export.py --out DIR --n 13 --ell 2 --ks 0,2,4,6
                                             [--seeds 1,2] [--exps 3,5,8] [--generic-exps]
                                             [--systems S_quad,P1only,random]
"""
import argparse
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP,
          os.path.join(_EXP, "E04_full_system_hilbert"),
          os.path.join(_EXP, "E02_official_baseline")]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import COMMIT, HOST                                    # noqa: E402
from aim3_macaulay import build_aim3, random_system                 # noqa: E402
from aim3.polysys_export import write_system                        # noqa: E402
from full_system_falls import default_exps, generic_exps            # noqa: E402
from table4_reproduce import dreg_boolean                           # noqa: E402


def substitute(poly, values):
    """Boolean polynomial with `values` (var -> 0/1) substituted in; XOR-reduced."""
    out = set()
    for m in poly:
        keep = []
        alive = True
        for v in m:
            if v in values:
                if values[v] == 0:
                    alive = False
                    break
            else:
                keep.append(v)
        if alive:
            out ^= {frozenset(keep)}
    return out


def fix_first_bits(eqs, nv, sol, k):
    """Fix variables 0..k-1 to their true values and renumber the rest to 0..nv-k-1."""
    values = {j: sol[j] for j in range(k)}
    ren = {v: v - k for v in range(k, nv)}
    out = []
    for f in eqs:
        g = substitute(f, values)
        out.append({frozenset(ren[v] for v in m) for m in g})
    return out, nv - k, sol[k:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--ell", type=int, default=2)
    ap.add_argument("--ks", default="0,2,4,6")
    ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--exps", default=None)
    ap.add_argument("--generic-exps", action="store_true")
    ap.add_argument("--systems", default="S_quad,P1only,random")
    ap.add_argument("--tag-prefix", default="E04h")
    a = ap.parse_args()

    n, ell = a.n, a.ell
    if a.exps:
        exps = [int(v) for v in a.exps.split(",")]
    else:
        exps = generic_exps(n, ell) if a.generic_exps else default_exps(n, ell)
    ks = [int(v) for v in a.ks.split(",")]
    seeds = [int(v) for v in a.seeds.split(",")]
    systems = a.systems.split(",")
    os.makedirs(a.out, exist_ok=True)
    print(f"# params: n={n} ell={ell} exps={exps} ks={ks} seeds={seeds} systems={systems} "
          f"out={a.out}")

    index = []
    for seed in seeds:
        nv, eqsP1, eqsAll, sol = build_aim3(n, exps, seed=seed)
        pool = {"S_quad": eqsAll, "P1only": eqsP1}
        if "random" in systems:
            pool["random"] = random_system(nv, len(eqsAll), seed=1000 + seed)[0]
        for k in ks:
            if k >= n:
                print(f"  skip k={k}: only n={n} pt bits exist")
                continue
            for name in systems:
                eqs = pool[name]
                # the random control has its own planted root; fixing OUR sol's bits would make it
                # inconsistent, so each system is cut with the root it actually has
                if name == "random":
                    req, rsol = random_system(nv, len(eqsAll), seed=1000 + seed)
                    sub, nv2, sub_sol = fix_first_bits(req, nv, rsol, k)
                else:
                    sub, nv2, sub_sol = fix_first_bits(eqs, nv, sol, k)
                # the planted root must survive the substitution
                for f in sub:
                    v = 0
                    for m in f:
                        t = 1
                        for j in m:
                            t &= sub_sol[j]
                        v ^= t
                    assert v == 0, f"{name}: the planted root did not survive k={k}"
                tag = f"{a.tag_prefix}_n{n}_l{ell}_s{seed}_k{k}_{name}"
                info = write_system(sub, nv2, a.out, tag, sysname=name, sol=sub_sol,
                                    meta=dict(experiment="E04h", n=n, ell=ell, seed=seed, k=k,
                                              exps=";".join(map(str, exps)), system=name,
                                              nvars_full=nv,
                                              b2_dreg=dreg_boolean(nv2, 3 * n * (ell + 1)),
                                              builder="legacy/aim3_macaulay.build_aim3 + fix pt bits",
                                              host=HOST, commit=COMMIT))
                index.append(info)
                print(f"  {tag}: {info['neqs']} eqs / {info['nvars']} vars "
                      f"(3n−k = {3*n}−{k}), B2 d_reg = {info['b2_dreg']}")
    with open(os.path.join(a.out, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1, sort_keys=True)
    print(f"# {len(index)} systems written")


if __name__ == "__main__":
    main()
