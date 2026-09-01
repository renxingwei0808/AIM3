#!/usr/bin/env sage
# -*- mode: python -*-
r"""GF(2) rank backend for the Boolean Macaulay matrices exported by `src/aim3/boolmac_export.py`.

    sage rank_backend.sage MATRIX.bin [--out RESULT.json] [--quiet]

Reads the `AIM3MAC1` sparse-row file, builds a *dense* `Matrix(GF(2), ...)` — Sage's dense GF(2)
matrices are M4RI — echelonises it once, and reports

    rank      = number of pivots
    rank_top  = number of pivots among the degree-D columns
    dim_low   = rank - rank_top = dim( rowspace  n  polynomials of degree <= D-1 )
    syz       = nrows - rank

**Column reversal.**  The exported file uses `boolmac.monomial_columns` order: degree grows with
the column index, so the degree-D monomials are the last `ncols - top_start` columns.  Sage
echelonises left to right, and the number of pivots inside a *leading* block of columns is the rank
of that block.  We therefore load exported column `j` at position `ncols - 1 - j`, which puts the
degree-D columns first; then `#pivots < ncols - top_start` is exactly `rank_top`.  This reproduces
`boolmac._eliminate`, which pivots on the highest bit (= the degree-D columns) first.

Exit status is 0 only if RESULT.json was written.
"""
import json
import os
import struct
import sys
import time
from array import array

MAGIC = b"AIM3MAC1"
HDR = struct.Struct("<4Q")


def parse_args(argv):
    if len(argv) < 2 or argv[1].startswith("--"):
        sys.stderr.write(__doc__)
        raise SystemExit(2)
    path = argv[1]
    out = path + ".rank.json"
    quiet = False
    i = 2
    while i < len(argv):
        if argv[i] == "--out":
            out = argv[i + 1]
            i += 2
        elif argv[i] == "--quiet":
            quiet = True
            i += 1
        else:
            sys.stderr.write("unknown option %s\n" % argv[i])
            raise SystemExit(2)
    return path, out, quiet


def load_dense(path, log):
    """-> (M, nrows, ncols, top_start, nnz, seconds) with the columns reversed as documented above."""
    fh = open(path, "rb")
    if fh.read(len(MAGIC)) != MAGIC:
        raise SystemExit("%s: not an AIM3MAC1 file" % path)
    nrows, ncols, top_start, nnz = HDR.unpack(fh.read(HDR.size))
    nrows, ncols, top_start, nnz = int(nrows), int(ncols), int(top_start), int(nnz)
    log("matrix %d x %d, %d nonzeros (%.2f per row), degree-D columns: %d"
        % (nrows, ncols, nnz, nnz / max(nrows, 1), ncols - top_start))
    t0 = time.time()
    M = matrix(GF(2), nrows, ncols)          # dense GF(2) == M4RI
    one = GF(2)(1)
    last = ncols - 1
    unpack = struct.Struct("<I").unpack
    for i in range(nrows):
        (cnt,) = unpack(fh.read(4))
        a = array("I")
        a.fromfile(fh, cnt)
        for j in a:
            M[i, last - j] = one
    fh.close()
    return M, nrows, ncols, top_start, nnz, time.time() - t0


def main(argv):
    path, out, quiet = parse_args(argv)

    def log(msg):
        if not quiet:
            sys.stdout.write("[rank_backend] %s\n" % msg)
            sys.stdout.flush()

    log("sage %s, file %s" % (sage.version.version, path))
    M, nrows, ncols, top_start, nnz, load_s = load_dense(path, log)
    log("loaded in %.1f s" % load_s)

    t0 = time.time()
    M.echelonize()
    ech_s = time.time() - t0
    piv = M.pivots()
    rank = len(piv)
    ntop = ncols - top_start                 # the degree-D columns now occupy positions [0, ntop)
    rank_top = sum(1 for p in piv if int(p) < ntop)
    dim_low = rank - rank_top
    log("echelonised in %.1f s: rank=%d rank_top=%d dim_low=%d syz=%d"
        % (ech_s, rank, rank_top, dim_low, nrows - rank))

    res = dict(path=os.path.abspath(path), nrows=int(nrows), ncols=int(ncols),
               top_start=int(top_start), nnz=int(nnz),
               rank=int(rank), rank_top=int(rank_top), dim_low=int(dim_low), syz=int(nrows - rank),
               load_seconds=float(round(load_s, 3)), echelon_seconds=float(round(ech_s, 3)),
               engine="sage/M4RI", sage_version=str(sage.version.version))
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
    log("wrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(int(main(sys.argv)))
