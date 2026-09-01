#!/usr/bin/env python3
"""Batch-2 task 0: the control that batch 1 listed as missing — S-boxes on DISJOINT variables.

`gate1_pairs.py` reports `cross_excess(B) = falls_D4(B) − Σ_{i∈B} falls_D4({i})`.  At ℓ = 3, n = 11
that came out 1639, and batch 1 argued from n = 13 and 15 that it is a small-parameter effect rather
than a coupling.  What that argument could not do is rule out a fault in the *measurement*: the
same-shape random control has no degree falls at all, so its `cross_excess` is trivially 0 and
proves nothing about the arithmetic.

This script supplies the missing control.  It puts k independent S-boxes on **pairwise disjoint
variable blocks** — S-box i on variables [2ni, 2n(i+1)) — inside one ring of 2nk variables, and runs
the identical subset-lattice measurement.  Here `cross_excess = 0` is forced by construction: the
ideals live on disjoint variable sets, so their degree-4 row spaces intersect only in the constants
and the fall counts must add.  **A non-zero answer would mean the instrument is wrong**, and every
`cross_excess` in batches 1 and 2 would have to be discarded.

Two configurations are worth running, and they test different things:

* `--k 2 --n 11`  — 44 variables, the same ring size and the same per-subset matrix size as the
  ℓ = 3, n = 11 pairs whose `cross_excess` was measured as 0 or n;
* `--k 4 --n 7`   — 56 variables, the same *number of blocks* as the ℓ = 3 instance where 1639
  appeared, which is the configuration the anomaly actually lives in.

Usage:
  python3 server/jobs/decoupled_control.py --n 11 --k 2 [--exps 3,4] [--backend sage]
                                           [--csv E04_server_decoupled.csv] [--budget SEC]
"""
import argparse
import itertools
import os
import sys
import time
from math import comb, gcd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record                                          # noqa: E402
from aim3_sbox_syz import sbox_system                               # noqa: E402
from aim3.boolmac import MacaulayTimeout                            # noqa: E402
from aim3.boolmac_export import macaulay_stats_backend              # noqa: E402
from aim3.structures import degenerate_class                        # noqa: E402

SCRATCH = os.environ.get("AIM3_SCRATCH")


def shift_poly(poly, offset):
    """Move a Boolean polynomial (set of frozensets of variable indices) up by `offset`."""
    return {frozenset(v + offset for v in m) for m in poly}


def build_decoupled(n, exps):
    """k S-boxes on disjoint variable blocks: S-box i uses variables [2ni, 2n(i+1)).

    Returns (nv, blocks) where blocks[i] is the list of 3n Boolean quadratics of S-box i.
    """
    blocks = []
    for i, e in enumerate(exps):
        _, P1, P2, P3 = sbox_system(n, e)
        blocks.append([shift_poly(p, 2 * n * i) for p in (P1 + P2 + P3)])
    return 2 * n * len(exps), blocks


def falls_of(eqs, nv, backend, budget, tag):
    t0 = time.time()
    try:
        s3 = macaulay_stats_backend(eqs, nv, 3, backend=backend, budget=budget, workdir=SCRATCH,
                                    tag=f"{tag}_D3", verbose=False)
        s4 = macaulay_stats_backend(eqs, nv, 4, backend=backend, budget=budget, workdir=SCRATCH,
                                    tag=f"{tag}_D4", verbose=False)
    except MacaulayTimeout:
        return None, None, None, None, 1, time.time() - t0
    return s4.dim_low - s3.rank, s4.dim_low, s3.rank, s3.syz, 0, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--k", type=int, default=2, help="number of independent S-boxes")
    ap.add_argument("--exps", default=None,
                    help="comma-separated exponents (default: the first k non-degenerate ones)")
    ap.add_argument("--backend", choices=["python", "sage"], default="sage")
    ap.add_argument("--budget", type=float, default=86400.0)
    ap.add_argument("--csv", default="E04_server_decoupled.csv")
    a = ap.parse_args()

    n, k = a.n, a.k
    if a.exps:
        exps = [int(v) for v in a.exps.split(",")]
    else:
        exps = [e for e in range(2, n - 1) if gcd(e, n) == 1 and not degenerate_class(n, e)][:k]
    if len(exps) != k:
        raise SystemExit(f"n={n}: need {k} non-degenerate exponents, found {exps}")

    nv, blocks = build_decoupled(n, exps)
    low_cols = sum(comb(nv, d) for d in range(4))
    print(f"# params: n={n} k={k} exps={exps} backend={a.backend} budget={a.budget}s")
    print(f"{k} S-boxes on DISJOINT variable blocks: {nv} variables ({2*n} each), "
          f"{3*n} quadratics each; degree-<=3 columns C({nv},<=3) = {low_cols}")
    print("cross_excess must be 0 for every subset -- it is forced by the disjoint supports.\n")

    single = {}
    bad = 0
    for size in range(1, k + 1):
        for sub in itertools.combinations(range(k), size):
            part = [e for i in sub for e in blocks[i]]
            tag = f"dec_n{n}_k{k}_{''.join(map(str, sub))}"
            f, low, r3, syz3, cens, secs = falls_of(part, nv, a.backend, a.budget, tag)
            if size == 1:
                single[sub[0]] = f
            ssum = None if any(single.get(i) is None for i in sub) else sum(single[i] for i in sub)
            cross = None if (f is None or ssum is None) else f - ssum
            sat = None if low is None else round(low / low_cols, 4)
            label = "+".join(map(str, sub))
            verdict = "" if cross is None else ("OK" if cross == 0 else "!! NON-ZERO -- INSTRUMENT FAULT")
            print(f"   blocks {{{label}}}: {len(part):4d} eqs  D3 rank={r3} syz={syz3} "
                  f"(8n*{size}={8*n*size})  falls_D4={f}  sum_singles={ssum}  "
                  f"cross_excess={cross}  saturation={sat}  {verdict}"
                  + ("   [CENSORED: lower bound]" if cens else ""), flush=True)
            if cross not in (None, 0):
                bad += 1
            record(a.csv, dict(
                n=n, k=k, exps=";".join(map(str, exps)), blocks=label, nblocks=len(sub),
                nvars=nv, neqs=len(part),
                D3_rank=r3 if r3 is not None else "", D3_syz=syz3 if syz3 is not None else "",
                syz_target=8 * n * size,
                D4_dim_low=low if low is not None else "", falls_D4=f if f is not None else "",
                sum_singles=ssum if ssum is not None else "",
                cross_excess=cross if cross is not None else "",
                low_cols=low_cols, saturation=sat if sat is not None else ""),
                seconds=secs, censored=cens, tier="B")
    print(f"\n# {bad} subset(s) with a non-zero cross_excess"
          + ("  -> the instrument is sound" if bad == 0 else
             "  -> THE INSTRUMENT IS NOT SOUND; every cross_excess in batches 1-2 is suspect"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
