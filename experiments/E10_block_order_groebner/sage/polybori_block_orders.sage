#!/usr/bin/env sage
# -*- mode: python -*-
r"""Block-order Gröbner bases in PolyBoRi — the measurement Magma 2.24 cannot make.

    sage polybori_block_orders.sage SYS.magma.m --order ORDER --out OUT.json [--quiet]

Server batch 1 established that Magma's GF(2) F4 runs in grevlex whatever monomial order the ring
declares and converts afterwards, so a "block-order step degree" read from its Faugère log is a
grevlex step degree (`experiments/E10_block_order_groebner/README.md`).  PolyBoRi computes
**directly** in the declared order and, with `prot=True`, prints one `Current Degree: d` line per
degree it works at — so `max d` is a genuine solving degree *in that order*.

`SYS.magma.m` is a file written by `src/aim3/polysys_export.py`: it defines `nvars` and `MONS`, the
list of monomial index lists (1-based).  Only those two are read, by text, so no Magma is involved.

ORDER is one of

    grevlex          degrevlex on all nvars variables
    block:K          block order, the FIRST K variables in the leading block  (K > rest)
    blocklast:K      block order, the LAST  K variables in the leading block  (rest < K);
                     implemented by rotating the variables so those K come first, which is an
                     algebra-preserving renaming — the step degrees are unaffected by the names

For the E06 kernel model (U = variables 1..3n, W = 3n+1..6n) `block:3n` is U > W and
`blocklast:3n` is W > U.  For an S_quad-shaped system (pt = 1..n) `block:n` is "pt block first".

The JSON records gb_size, gb_max_degree, the parsed step degrees, seconds and the order actually
built, so a run can be checked against what was asked for.
"""
import json
import os
import re
import sys
import time


def read_system(path):
    """(nvars, MONS) from a polysys_export Magma file, by text -- 1-based indices as written."""
    txt = open(path).read()
    m = re.search(r"nvars\s*:=\s*(\d+)\s*;", txt)
    if not m:
        raise SystemExit("%s: no nvars" % path)
    nvars = int(m.group(1))
    body = txt.split("MONS := [", 1)[1].rsplit("];", 1)[0]
    import ast
    mons = ast.literal_eval("[" + body + "]")
    return nvars, mons


def build_ring(nvars, order):
    """(ring, permutation) -- perm[i] is the NEW position of original variable i (0-based)."""
    perm = list(range(nvars))
    if order == "grevlex":
        R = BooleanPolynomialRing(nvars, 'x', order='degrevlex')
        return R, perm, "degrevlex"
    kind, k = order.split(":")
    k = int(k)
    if not 0 < k < nvars:
        raise SystemExit("block size %d out of range for %d variables" % (k, nvars))
    if kind == "blocklast":
        # rotate so the last k variables come first; a renaming, not a change of algebra
        perm = [(i - (nvars - k)) % nvars for i in range(nvars)]
    elif kind != "block":
        raise SystemExit("unknown order %s" % order)
    T = TermOrder('degrevlex', k) + TermOrder('degrevlex', nvars - k)
    R = BooleanPolynomialRing(nvars, 'x', order=T)
    return R, perm, str(R.term_order())


def main(argv):
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        raise SystemExit(2)
    path = argv[1]
    order, out, quiet = "grevlex", None, False
    i = 2
    while i < len(argv):
        if argv[i] == "--order":
            order = argv[i + 1]; i += 2
        elif argv[i] == "--out":
            out = argv[i + 1]; i += 2
        elif argv[i] == "--quiet":
            quiet = True; i += 1
        else:
            raise SystemExit("unknown option %s" % argv[i])
    out = out or (path + ".polybori.json")

    nvars, mons = read_system(path)
    R, perm, ordname = build_ring(nvars, order)
    x = R.gens()
    F = []
    for poly in mons:
        p = R(int(0))
        for m in poly:
            t = R(int(1))
            for v in m:                      # MONS is 1-based
                t = t * x[perm[int(v) - 1]]
            p = p + t
        if not p.is_zero():
            F.append(p)
    print("### system %s: %d polynomials in %d Boolean variables" % (os.path.basename(path),
                                                                     len(F), nvars))
    print("### order requested %s, ring order %s" % (order, ordname))
    sys.stdout.flush()

    I = R.ideal(F)
    t0 = time.time()
    G = I.groebner_basis(prot=True)
    secs = time.time() - t0
    degs = [int(g.deg()) for g in G]
    print("### GB size %d, max degree %d, %.3f s" % (len(G), max(degs) if degs else 0, secs))

    res = dict(path=os.path.abspath(path), order=order, ring_order=ordname,
               nvars=int(nvars), neqs=len(F), gb_size=len(G),
               gb_max_degree=int(max(degs)) if degs else 0,
               trivial_gb=bool(len(G) == 1 and G[0] == R(int(1))),
               seconds=float(round(secs, 3)), tool="polybori",
               sage_version=str(sage.version.version))
    with open(out, "w") as fh:
        json.dump(res, fh, indent=int(1), sort_keys=True)
    print("### wrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(int(main(sys.argv)))
