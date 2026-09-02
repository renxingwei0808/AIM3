#!/usr/bin/env python3
"""Batch-2 task 2.2: is the D=5 fall count just the multiples of the D=4 cubics?

Batch 1 measured, on four runs and exactly,

    falls_D5 = n(4n + 7)   for generic e,      falls_D5 = n(4n − 6)   for the {−3, (−3)⁻¹} class,

where `falls_D5 = dim(span_{D=5} ∩ deg≤4) − rank_{D=4}`.  The obvious candidate explanation is that
the D=5 pass mostly re-uses what D=4 already found: the 2n new cubics `c_j` times the 2n variables
`z_k` are 4n² degree-≤4 elements of the D=5 row space, and 4n² is the leading term of both formulas.
This script measures how much of `falls_D5` those products actually account for.

    X        = dim( span{ z_k · c_j }  modulo  rowspace(Mac_4) )
             = rank(Mac₄ rows ∪ products) − rank(Mac₄ rows)
    net new  = falls_D5 − X

If the batch-1 arithmetic `n(4n+7) = (2n)(2n) + 7n` is the right reading, X = 4n² and net = 7n; if
the products are not independent, X < 4n² and the shortfall is a count of linear relations among
them.  Both ranks are taken over the degree-≤4 column space (a product `z_k·c_j` has degree ≤ 4, so
nothing leaves it), with the Sage/M4RI backend.

The 2n cubics themselves are obtained the way Tier-A batch 2 does it: `fall_space(eqs, nv, 4)`
gives a basis of rowspace(Mac₄) ∩ deg≤3, and `reduce_modulo` against the Mac₃ rows leaves exactly
the new ones.

Usage:
  python3 server/jobs/d5_decomposition.py --n 14 --e 3 [--backend sage] [--falls-d5 882]
                                          [--csv E03_server_d5_decomposition.csv]
"""
import argparse
import csv
import itertools
import os
import sys
import time
from math import comb

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record                                          # noqa: E402
from aim3_sbox_syz import sbox_system                               # noqa: E402
from aim3.boolmac import (monomial_columns, poly_rows, fall_space,  # noqa: E402
                          reduce_modulo, macaulay_stats)
from aim3.boolmac_export import rank_sparse_sage                    # noqa: E402
from aim3.structures import degenerate_class                        # noqa: E402

SCRATCH = os.environ.get("AIM3_SCRATCH")


def int_to_monomials(x, inv_col):
    """Bit-packed row over the degree-<=D columns -> the set of monomial bitmasks it contains."""
    out = []
    while x:
        b = (x & -x).bit_length() - 1
        out.append(inv_col[b])
        x &= x - 1
    return out


def lookup_falls_d5(n, e):
    """falls_D5 from results/E03_server_D5.csv if batch 1 already measured this (n, e)."""
    p = os.path.join(ROOT, "results", "E03_server_D5.csv")
    if not os.path.exists(p):
        return None
    for r in csv.DictReader(open(p)):
        if (int(r["n"]), int(r["e"])) == (n, e) and r["system"] == "S_box" and r["falls_D5"]:
            return int(r["falls_D5"])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--e", type=int, required=True)
    ap.add_argument("--backend", choices=["python", "sage"], default="sage")
    ap.add_argument("--falls-d5", type=int, default=None,
                    help="override; by default read from results/E03_server_D5.csv")
    ap.add_argument("--csv", default="E03_server_d5_decomposition.csv")
    a = ap.parse_args()
    n, e = a.n, a.e
    nv = 2 * n
    t_all = time.time()

    cls = ";".join(degenerate_class(n, e)) or "generic"
    print(f"# params: n={n} e={e} class={cls} backend={a.backend}")

    fld, P1, P2, P3 = sbox_system(n, e)
    eqs = P1 + P2 + P3

    # ---- the 2n new cubics
    t0 = time.time()
    col4, ncol4, top4 = monomial_columns(nv, 4)
    inv4 = {v: k for k, v in col4.items()}
    mac3 = list(poly_rows(eqs, nv, 3))                 # D=3 columns are a prefix of the D=4 ones
    low4 = fall_space(eqs, nv, 4)["low"]
    cubics, ncub = reduce_modulo(low4, mac3)
    print(f"   fall space at D=4: {len(low4)} elements; new modulo Mac3: {ncub} "
          f"(2n = {2*n}, 2n+1 = {2*n+1})   [{time.time()-t0:.1f}s]")

    # ---- the products z_k * c_j, still of degree <= 4
    t0 = time.time()
    prods = []
    for c in cubics:
        mons = int_to_monomials(c, inv4)
        for k in range(nv):
            bit = 1 << k
            acc = set()
            for m in mons:
                acc ^= {m | bit}
            if not acc or any(x.bit_count() > 4 for x in acc):
                continue
            prods.append(sorted(col4[x] for x in acc))
    print(f"   products z_k*c_j: {len(prods)} of {ncub*nv} kept "
          f"(a product is dropped only if it vanishes or leaves degree 4)   [{time.time()-t0:.1f}s]")

    # ---- ranks over the degree-<=4 column space
    def mac4_rows():
        for r in poly_rows(eqs, nv, 4, col=col4):
            out = []
            x = r
            while x:
                b = (x & -x).bit_length() - 1
                out.append(b)
                x &= x - 1
            yield out

    t0 = time.time()
    if a.backend == "sage":
        r_mac4 = rank_sparse_sage(mac4_rows(), ncol4, workdir=SCRATCH, tag=f"d5dec_n{n}_e{e}_A")
        r_both = rank_sparse_sage(itertools.chain(mac4_rows(), iter(prods)), ncol4,
                                  workdir=SCRATCH, tag=f"d5dec_n{n}_e{e}_AP")
    else:
        from aim3.boolmac import span_rank
        A = [sum(1 << c for c in row) for row in mac4_rows()]
        P = [sum(1 << c for c in row) for row in prods]
        r_mac4, r_both = span_rank(A), span_rank(A + P)
    X = r_both - r_mac4
    print(f"   rank(Mac4) = {r_mac4}, rank(Mac4 u products) = {r_both}  ->  X = {X}   "
          f"[{time.time()-t0:.1f}s]")

    falls5 = a.falls_d5 if a.falls_d5 is not None else lookup_falls_d5(n, e)
    net = None if falls5 is None else falls5 - X
    # the natural ceiling is the NUMBER of products, ncub*nv; 4n^2 is that ceiling only when the
    # cubic count is the generic 2n, so both are reported
    ceil = len(prods)
    pred = 4 * n * n
    print(f"   ceiling (number of products) = {ceil};  X = {X}  ->  "
          f"{'all independent' if X == ceil else f'{ceil - X} linear relations among the products'}")
    print(f"   4n^2 = {pred} (the batch-1 reading, valid only when #cubics = 2n = {2*n}; "
          f"here #cubics = {ncub});  X - 4n^2 = {X - pred}")
    if falls5 is not None:
        print(f"   falls_D5 = {falls5} = {falls5/n:g}n;  net new quartics = {falls5} - {X} = {net}"
              f"  ({net/n:g}n)" if net is not None else "")
    else:
        print("   falls_D5 not available -- run server/jobs/e03_d5.py for this (n, e) first")

    record(a.csv, dict(n=n, e=e, deg_class=cls, nvars=nv,
                       cubics=ncub, products=len(prods),
                       rank_mac4=r_mac4, rank_mac4_plus_products=r_both, X=X,
                       products_ceiling=len(prods), relations_among_products=len(prods) - X,
                       four_n_squared=pred, X_minus_4n2=X - pred,
                       falls_D5=falls5 if falls5 is not None else "",
                       net_new_quartics=net if net is not None else "",
                       net_over_n=round(net / n, 4) if net is not None else "",
                       cols_le4=ncol4),
           seconds=time.time() - t_all, tier="B")


if __name__ == "__main__":
    main()
