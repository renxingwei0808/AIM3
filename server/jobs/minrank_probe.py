#!/usr/bin/env python3
"""Batch-2 task 3.2: does the first layer's rank condition carry any information?

ROADMAP T2.1 proposes reading the first layer as a MinRank/Kipnis–Shamir instance.  Concretely, the
kernel-model relation of first-layer S-box i,

    B_{e_i}(a_i, x_i) = a_i,        x_i = A_i·pt + b_i,

is, coordinate by coordinate, **affine in the pt bits with coefficients affine in the a_i bits**.
Stacking the two first-layer S-boxes gives

    N(W) · (pt, 1)^T = 0,     N of size 2n × (n+1), entries affine in W = (a_0, a_1),

so the true W makes N drop rank: `(pt, 1)` is a kernel vector, hence rank N ≤ n.  MinRank is only
an attack if that is **unusual** — if a random W gives rank n+1.  This script measures exactly that:

  1. rank N at the true W (expected n, i.e. corank 1) and whether `(pt,1)` really spans the kernel;
  2. the rank distribution over uniformly random W — **if a random W also gives rank ≤ n, the rank
     condition is free and MinRank has nothing to key on**;
  3. `--export DIR` writes the KS system — the 2n equations `N(W)·(u,1)^T = 0` in the 3n unknowns
     (u, a_0, a_1), i.e. the first layer with the t variables gone — in msolve and Magma syntax, so
     `run_f4_degrees.py` can put its `d_solv` next to `S_quad`'s.

Usage:
  python3 server/jobs/minrank_probe.py --n 9 [--samples 2000] [--seed 1] [--generic-exps]
                                       [--export DIR] [--csv E06_server_minrank.csv]
"""
import argparse
import json
import os
import random
import sys
import time
from collections import Counter
from math import gcd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP,
          os.path.join(_EXP, "E06_kernel_bilinear_model")]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record, COMMIT, HOST                            # noqa: E402
from aim3.structures import degenerate_class                        # noqa: E402
from aim3.polysys_export import write_system                        # noqa: E402
from kernel_bilinear_toy import build_kernel_model                  # noqa: E402


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


def affine_rows(polys, n, wvals):
    """Substitute the a-bits, then read each polynomial as a row (coef of pt_0..pt_{n-1}, const).

    Raises if a substituted polynomial is not affine in pt -- that would mean the layout assumed
    here is wrong, and a silent wrong answer is worse than a crash.
    """
    rows = []
    for f in polys:
        g = substitute(f, wvals)
        row = 0
        for m in g:
            if len(m) == 0:
                row |= 1 << n                       # the constant column
            elif len(m) == 1:
                (v,) = tuple(m)
                assert v < n, f"variable {v} is not a pt bit after substitution"
                row ^= 1 << v
            else:
                raise AssertionError(f"non-affine monomial {sorted(m)} after substituting W")
        rows.append(row)
    return rows


def rank_and_kernel(rows, ncols):
    """(rank, basis of the right kernel as ints over ncols columns) of a GF(2) matrix."""
    piv = {}
    for r in rows:
        rr, tag = r, 0
        while rr:
            b = rr.bit_length() - 1
            if b in piv:
                rr ^= piv[b][0]
                tag ^= piv[b][1]
            else:
                piv[b] = (rr, tag)
                break
    rank = len(piv)
    # right kernel by back-substitution over the column space
    pivcols = sorted(piv)
    free = [c for c in range(ncols) if c not in piv]
    ker = []
    for fc in free:
        vec = 1 << fc
        for c in sorted(pivcols):
            row = piv[c][0]
            if bin(row & vec).count("1") % 2:
                vec ^= 1 << c
        ker.append(vec)
    return rank, ker


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--samples", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exps", default=None)
    ap.add_argument("--generic-exps", action="store_true")
    ap.add_argument("--export", default=None)
    ap.add_argument("--csv", default="E06_server_minrank.csv")
    a = ap.parse_args()
    n = a.n
    t0 = time.time()

    if a.exps:
        exps = [int(v) for v in a.exps.split(",")]
    elif a.generic_exps:
        cand = [e for e in range(2, n) if gcd(e, n) == 1 and not degenerate_class(n, e)]
        if len(cand) < 3:
            raise SystemExit(f"n={n}: only {len(cand)} non-degenerate exponents {cand}, need 3")
        exps = cand[:3]
    else:
        exps = [e for e in range(2, n) if gcd(e, n) == 1][:3]

    nv, eqs, sol = build_kernel_model(n, exps)
    # layout: pt 0..n-1 | t0 n..2n-1 | t1 2n..3n-1 | a0 3n..4n-1 | a1 4n..5n-1 | a2 5n..6n-1
    # equation order: B_0 (n), K_0 (n), B_1 (n), K_1 (n), B_2 (n), K_2 (n)
    B = eqs[0:n] + eqs[2 * n:3 * n]                 # the two first-layer B blocks, 2n equations
    wvars = list(range(3 * n, 5 * n))               # a_0, a_1
    classes = ";".join(f"{e}:{';'.join(degenerate_class(n, e)) or 'generic'}" for e in exps)
    print(f"# params: n={n} exps={exps} ({classes}) samples={a.samples} seed={a.seed}")
    print(f"N(W) is {len(B)} x {n+1}: {2*n} first-layer B coordinates, affine in pt, "
          f"coefficients affine in W = (a_0, a_1) ({2*n} bits)")

    # ---- 1. the true W
    wtrue = {v: sol[v] for v in wvars}
    rows = affine_rows(B, n, wtrue)
    rk, ker = rank_and_kernel(rows, n + 1)
    ptvec = sum(sol[j] << j for j in range(n)) | (1 << n)     # (pt, 1)
    in_ker = all(bin(r & ptvec).count("1") % 2 == 0 for r in rows)
    print(f"   true W: rank N = {rk} (n = {n}, corank {n + 1 - rk}); "
          f"(pt,1) is in the kernel: {in_ker}; kernel dimension {len(ker)}")
    assert in_ker, "(pt,1) is not a kernel vector -- the layout assumed here is wrong"

    # ---- 2. random W
    rng = random.Random(a.seed)
    hist = Counter()
    t1 = time.time()
    for _ in range(a.samples):
        w = {v: rng.getrandbits(1) for v in wvars}
        hist[rank_and_kernel(affine_rows(B, n, w), n + 1)[0]] += 1
    print(f"   random W ({a.samples} samples, {time.time()-t1:.1f}s): rank distribution "
          + ", ".join(f"{r}:{c} ({100*c/a.samples:.1f}%)" for r, c in sorted(hist.items())))
    le_n = sum(c for r, c in hist.items() if r <= n)
    print(f"   => rank <= n for {le_n}/{a.samples} random W ({100*le_n/a.samples:.2f}%)"
          + ("  -- the rank condition is FREE, MinRank has nothing to key on"
             if le_n == a.samples else
             "  -- the true W is distinguished by its rank"))

    # ---- 3. the KS system
    exported = ""
    if a.export:
        # the 2n B equations live on (pt, a_0, a_1); renumber those 3n variables to 0..3n-1
        keep = list(range(0, n)) + list(range(3 * n, 5 * n))
        ren = {v: i for i, v in enumerate(keep)}
        ks = [{frozenset(ren[v] for v in m) for m in f} for f in B]
        kssol = [sol[v] for v in keep]
        for f in ks:
            v = 0
            for m in f:
                t = 1
                for j in m:
                    t &= kssol[j]
                v ^= t
            assert v == 0, "the KS system does not vanish at the true key"
        tag = f"E06ks_n{n}_s{a.seed}_KS"
        info = write_system(ks, 3 * n, a.export, tag, sysname="KS_first_layer", sol=kssol,
                            meta=dict(experiment="E06ks", n=n, ell=2, seed=a.seed, k=0,
                                      exps=";".join(map(str, exps)), system="KS_first_layer",
                                      host=HOST, commit=COMMIT))
        exported = info["path"] if "path" in info else tag
        print(f"   [export] {tag}: {info['neqs']} eqs / {info['nvars']} vars -> {a.export}")

    record(a.csv, dict(n=n, exps=";".join(map(str, exps)), exp_classes=classes, seed=a.seed,
                       N_rows=len(B), N_cols=n + 1, rank_true_W=rk, corank_true_W=n + 1 - rk,
                       pt_in_kernel=int(in_ker), samples=a.samples,
                       rank_hist=";".join(f"{r}:{c}" for r, c in sorted(hist.items())),
                       random_rank_le_n=le_n,
                       random_rank_le_n_frac=round(le_n / a.samples, 6),
                       ks_exported=exported),
           seconds=time.time() - t0, tier="B")


if __name__ == "__main__":
    main()
