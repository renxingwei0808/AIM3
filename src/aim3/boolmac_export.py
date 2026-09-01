"""Export a Boolean Macaulay matrix for an external rank engine (Sage / M4RI).

`boolmac.macaulay_stats` streams the rows through a pure-Python bit-integer elimination.  That is
exact and dependency-free, but O(rank * rows * cols / 64) in *Python* word operations, which is
days at n = 13, D = 4.  This module writes the *same* matrix — same rows, same multipliers, same
column order — to a file that `experiments/E04_full_system_hilbert/sage/rank_backend.sage` reads
into an M4RI dense matrix.

**The column order is `boolmac.monomial_columns`**: columns are numbered by increasing degree, so
the degree-D monomials occupy `[top_start, ncol)` — the *highest* indices.  The Sage backend
reverses that when it builds the matrix, so that after a left-to-right echelonisation

    rank        = #pivots
    rank_top    = #pivots in the first (ncol - top_start) positions   (= the degree-D columns)
    dim_low     = rank - rank_top                                     (= boolmac's `dim_low`)

which is exactly the quantity `boolmac._eliminate` counts when it pivots on the highest bit first.
`tests/test_core.py::test_export_matches_boolmac` checks row-for-row agreement with `poly_rows`.

File format (little endian, `AIM3MAC1`):

    magic    8 bytes  b"AIM3MAC1"
    header   4 x uint64   nrows, ncols, top_start, nnz
    rows     nrows x [ uint32 count, count x uint32 column-index (ascending) ]

A `<path>.meta.json` sidecar carries the same header plus nv, D and mult_deg for humans.
"""
from __future__ import annotations
import itertools, json, os, struct, subprocess, sys, tempfile, time
from array import array

from .boolmac import monomial_columns, MacaulayResult

__all__ = ["poly_row_indices", "export_macaulay", "read_macaulay", "macaulay_stats_sage",
           "SageBackendError", "sage_backend_script"]

MAGIC = b"AIM3MAC1"
_HDR = struct.Struct("<4Q")


class SageBackendError(RuntimeError):
    """The Sage rank backend failed; the caller must NOT fall back silently to another engine."""


def _multipliers(nv: int, mult_deg: int):
    """The multiplier monomials of `boolmac.poly_rows`, in the same order."""
    out = []
    for d in range(mult_deg + 1):
        for m in itertools.combinations(range(nv), d):
            mask = 0
            for i in m:
                mask |= 1 << i
            out.append(mask)
    return out


def poly_row_indices(eqs, nv: int, D: int, col=None, mult_deg=None):
    """Yield each Macaulay row as a sorted list of column indices.

    Mirrors `boolmac.poly_rows` exactly — same multiplier order, same "drop the row if it vanishes
    or leaves degree D" rules — but never materialises the big integer.
    """
    if col is None:
        col, _, _ = monomial_columns(nv, D)
    if mult_deg is None:
        mult_deg = D - 2
    mults = _multipliers(nv, mult_deg)
    for f in eqs:
        fm = []
        for m in f:
            mask = 0
            for i in m:
                mask |= 1 << i
            fm.append(mask)
        for mu in mults:
            acc = set()
            for vm in fm:
                acc ^= {vm | mu}
            if not acc:
                continue
            over = False
            idx = []
            for key in acc:
                if key.bit_count() > D:
                    over = True
                    break
                idx.append(col[key])
            if over or not idx:
                continue
            idx.sort()
            yield idx


def export_macaulay(eqs, nv: int, D: int, path: str, mult_deg=None):
    """Write the degree-D Macaulay matrix of `eqs` to `path`; return the meta dict."""
    col, ncol, top_start = monomial_columns(nv, D)
    if mult_deg is None:
        mult_deg = D - 2
    nrows = nnz = 0
    t0 = time.time()
    with open(path, "wb") as fh:
        fh.write(MAGIC)
        fh.write(_HDR.pack(0, 0, 0, 0))                     # placeholder, rewritten below
        for idx in poly_row_indices(eqs, nv, D, col=col, mult_deg=mult_deg):
            fh.write(struct.pack("<I", len(idx)))
            array("I", idx).tofile(fh)
            nrows += 1
            nnz += len(idx)
        fh.seek(len(MAGIC))
        fh.write(_HDR.pack(nrows, ncol, top_start, nnz))
    meta = dict(path=os.path.abspath(path), nv=nv, D=D, mult_deg=mult_deg, nrows=nrows, ncols=ncol,
                top_start=top_start, nnz=nnz, neqs=len(eqs), export_seconds=round(time.time() - t0, 3),
                bytes=os.path.getsize(path))
    with open(path + ".meta.json", "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    return meta


def read_macaulay(path):
    """(meta, iterator of rows-as-index-lists) — the reader used by the tests and by Sage."""
    fh = open(path, "rb")
    if fh.read(len(MAGIC)) != MAGIC:
        fh.close()
        raise ValueError(f"{path}: not an AIM3MAC1 file")
    nrows, ncols, top_start, nnz = _HDR.unpack(fh.read(_HDR.size))
    meta = dict(nrows=nrows, ncols=ncols, top_start=top_start, nnz=nnz)

    def rows():
        try:
            for _ in range(nrows):
                (cnt,) = struct.unpack("<I", fh.read(4))
                a = array("I")
                a.fromfile(fh, cnt)
                yield a
        finally:
            fh.close()
    return meta, rows()


def sage_backend_script():
    """Absolute path of experiments/E04_full_system_hilbert/sage/rank_backend.sage."""
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    return os.path.join(root, "experiments", "E04_full_system_hilbert", "sage", "rank_backend.sage")


def macaulay_stats_sage(eqs, nv: int, D: int, mult_deg=None, workdir=None, keep=False,
                        sage=None, timeout=None, tag="mac", verbose=True) -> MacaulayResult:
    """Same return value as `boolmac.macaulay_stats`, computed by Sage/M4RI.

    `workdir` (default: a temp dir, removed unless `keep`) holds the exported matrix; matrices are
    scratch data and never belong in the repository.  Raises `SageBackendError` on any failure —
    a wrong rank must never be produced by silently falling back.
    """
    sage = sage or os.environ.get("AIM3_SAGE") or "sage"
    script = sage_backend_script()
    if not os.path.exists(script):
        raise SageBackendError(f"missing Sage backend script {script}")
    tmp = workdir or tempfile.mkdtemp(prefix="aim3mac_")
    os.makedirs(tmp, exist_ok=True)
    mat = os.path.join(tmp, f"{tag}_nv{nv}_D{D}.bin")
    out = mat + ".rank.json"
    try:
        meta = export_macaulay(eqs, nv, D, mat, mult_deg=mult_deg)
        if verbose:
            print(f"     [sage] exported {meta['nrows']}x{meta['ncols']} ({meta['nnz']} nnz, "
                  f"{meta['bytes']/2**20:.1f} MiB) in {meta['export_seconds']:.1f}s", flush=True)
        cmd = [sage, script, mat, "--out", out]
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0 or not os.path.exists(out):
            raise SageBackendError(
                f"{' '.join(cmd)} failed (rc={proc.returncode})\nstdout:\n{proc.stdout[-4000:]}\n"
                f"stderr:\n{proc.stderr[-4000:]}")
        with open(out) as fh:
            res = json.load(fh)
        if verbose:
            print(f"     [sage] rank={res['rank']} rank_top={res['rank_top']} "
                  f"dim_low={res['dim_low']} ({time.time()-t0:.1f}s wall, "
                  f"load {res.get('load_seconds')}s, echelon {res.get('echelon_seconds')}s)", flush=True)
        if res["nrows"] != meta["nrows"] or res["ncols"] != meta["ncols"]:
            raise SageBackendError(f"shape mismatch: exported {meta}, Sage reported {res}")
        return MacaulayResult(res["nrows"], res["ncols"], res["rank"], res["dim_low"])
    finally:
        if not keep and workdir is None:
            for p in (mat, mat + ".meta.json", out):
                try:
                    os.remove(p)
                except OSError:
                    pass
            try:
                os.rmdir(tmp)
            except OSError:
                pass


if __name__ == "__main__":                                    # tiny CLI for debugging an export
    import argparse
    ap = argparse.ArgumentParser(description="dump an exported Macaulay matrix header")
    ap.add_argument("path")
    a = ap.parse_args()
    meta, rows = read_macaulay(a.path)
    print(json.dumps(meta, indent=1))
    for i, r in enumerate(rows):
        if i >= 3:
            break
        print(f"row {i}: {len(r)} nnz, first 10 columns {list(r[:10])}")


# --------------------------------------------------------------------------- point-evaluation
# `E03_syzygy_and_degree_fall/sbox_degree_fall.relation_dim` also needs a GF(2) rank: the matrix of
# all monomials of degree <= deg evaluated at every graph point.  Its column set and order are the
# same `monomial_columns(nv, deg)`, so the same file format and the same Sage backend serve it.
# A point contributes a 1 in column `m` iff supp(m) is a subset of supp(point), so the nonzeros of
# a row are exactly the subsets of the point's support of size <= deg -- enumerating those is
# C(w,<=deg) work per point instead of scanning all C(nv,<=deg) monomials.

def point_row_indices(points, nv: int, deg: int, col=None):
    """Yield each point's evaluation row as a sorted list of column indices."""
    if col is None:
        col, _, _ = monomial_columns(nv, deg)
    for p in points:
        supp = []
        q = p
        while q:
            b = (q & -q).bit_length() - 1
            supp.append(b)
            q &= q - 1
        idx = []
        for d in range(min(deg, len(supp)) + 1):
            for c in itertools.combinations(supp, d):
                mask = 0
                for i in c:
                    mask |= 1 << i
                idx.append(col[mask])
        idx.sort()
        yield idx


def export_point_matrix(points, nv: int, deg: int, path: str):
    """Write the (points x monomials of degree <= deg) evaluation matrix; return the meta dict."""
    col, ncol, top_start = monomial_columns(nv, deg)
    nrows = nnz = 0
    t0 = time.time()
    with open(path, "wb") as fh:
        fh.write(MAGIC)
        fh.write(_HDR.pack(0, 0, 0, 0))
        for idx in point_row_indices(points, nv, deg, col=col):
            fh.write(struct.pack("<I", len(idx)))
            array("I", idx).tofile(fh)
            nrows += 1
            nnz += len(idx)
        fh.seek(len(MAGIC))
        fh.write(_HDR.pack(nrows, ncol, top_start, nnz))
    meta = dict(path=os.path.abspath(path), nv=nv, D=deg, kind="point-evaluation", nrows=nrows,
                ncols=ncol, top_start=top_start, nnz=nnz, export_seconds=round(time.time() - t0, 3),
                bytes=os.path.getsize(path))
    with open(path + ".meta.json", "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    return meta


def _run_sage(mat, out, sage=None, timeout=None):
    sage = sage or os.environ.get("AIM3_SAGE") or "sage"
    script = sage_backend_script()
    if not os.path.exists(script):
        raise SageBackendError(f"missing Sage backend script {script}")
    cmd = [sage, script, mat, "--out", out]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0 or not os.path.exists(out):
        raise SageBackendError(f"{' '.join(cmd)} failed (rc={proc.returncode})\n"
                               f"stdout:\n{proc.stdout[-4000:]}\nstderr:\n{proc.stderr[-4000:]}")
    with open(out) as fh:
        return json.load(fh)


def relation_dim_sage(points, nv: int, deg: int, workdir=None, keep=False, sage=None, timeout=None,
                      tag="pts", verbose=True):
    """Sage/M4RI version of `relation_dim`: -> (dim, ncols, rank), identical semantics."""
    tmp = workdir or tempfile.mkdtemp(prefix="aim3pts_")
    os.makedirs(tmp, exist_ok=True)
    mat = os.path.join(tmp, f"{tag}_nv{nv}_deg{deg}.bin")
    out = mat + ".rank.json"
    try:
        meta = export_point_matrix(points, nv, deg, mat)
        if verbose:
            print(f"     [sage] point matrix {meta['nrows']}x{meta['ncols']} ({meta['nnz']} nnz, "
                  f"{meta['bytes']/2**20:.1f} MiB) exported in {meta['export_seconds']:.1f}s", flush=True)
        res = _run_sage(mat, out, sage=sage, timeout=timeout)
        if verbose:
            print(f"     [sage] rank={res['rank']} -> dim = {res['ncols']} - {res['rank']} = "
                  f"{res['ncols'] - res['rank']} ({res.get('echelon_seconds')}s echelon)", flush=True)
        return res["ncols"] - res["rank"], res["ncols"], res["rank"]
    finally:
        if not keep and workdir is None:
            for p in (mat, mat + ".meta.json", out):
                try:
                    os.remove(p)
                except OSError:
                    pass
            try:
                os.rmdir(tmp)
            except OSError:
                pass


def macaulay_stats_backend(eqs, nv: int, D: int, backend="python", budget=None, mult_deg=None,
                           workdir=None, keep=False, sage=None, tag="mac", verbose=True):
    """Dispatch `boolmac.macaulay_stats` (backend='python') or the Sage/M4RI engine ('sage').

    Both raise `boolmac.MacaulayTimeout` when `budget` is exceeded.  For the Sage engine there is no
    meaningful partial elimination state, so the attached partial result carries rank 0 -- a valid
    lower bound, and `censored=1` in the CSV marks it as such (CLAUDE.md ground rule 4).
    """
    from .boolmac import macaulay_stats, MacaulayTimeout
    if backend == "python":
        return macaulay_stats(eqs, nv, D, mult_deg=mult_deg, budget=budget)
    if backend != "sage":
        raise ValueError(f"unknown backend {backend!r} (expected 'python' or 'sage')")
    try:
        return macaulay_stats_sage(eqs, nv, D, mult_deg=mult_deg, workdir=workdir, keep=keep,
                                   sage=sage, timeout=budget, tag=tag, verbose=verbose)
    except subprocess.TimeoutExpired:
        _, ncol, _ = monomial_columns(nv, D)
        raise MacaulayTimeout(MacaulayResult(0, ncol, 0, 0), budget or 0.0)


def export_sparse(rows, ncols, path, top_start=None, kind="sparse"):
    """Write an arbitrary sparse GF(2) matrix (rows = iterables of column indices) in AIM3MAC1."""
    if top_start is None:
        top_start = ncols                       # no "top block": rank_top == 0, dim_low == rank
    nrows = nnz = 0
    t0 = time.time()
    with open(path, "wb") as fh:
        fh.write(MAGIC)
        fh.write(_HDR.pack(0, 0, 0, 0))
        for idx in rows:
            idx = sorted(idx)
            fh.write(struct.pack("<I", len(idx)))
            array("I", idx).tofile(fh)
            nrows += 1
            nnz += len(idx)
        fh.seek(len(MAGIC))
        fh.write(_HDR.pack(nrows, ncols, top_start, nnz))
    meta = dict(path=os.path.abspath(path), kind=kind, nrows=nrows, ncols=ncols,
                top_start=top_start, nnz=nnz, export_seconds=round(time.time() - t0, 3),
                bytes=os.path.getsize(path))
    with open(path + ".meta.json", "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    return meta


def rank_sparse_sage(rows, ncols, workdir=None, keep=False, sage=None, timeout=None, tag="rk",
                     verbose=True):
    """Rank over GF(2) of an arbitrary sparse row set, via the Sage/M4RI backend."""
    tmp = workdir or tempfile.mkdtemp(prefix="aim3rk_")
    os.makedirs(tmp, exist_ok=True)
    mat = os.path.join(tmp, f"{tag}_{ncols}.bin")
    out = mat + ".rank.json"
    try:
        meta = export_sparse(rows, ncols, mat)
        if verbose:
            print(f"     [sage] {meta['nrows']}x{meta['ncols']} ({meta['nnz']} nnz, "
                  f"{meta['bytes']/2**20:.1f} MiB) exported in {meta['export_seconds']:.1f}s",
                  flush=True)
        res = _run_sage(mat, out, sage=sage, timeout=timeout)
        if verbose:
            print(f"     [sage] rank={res['rank']} ({res.get('echelon_seconds')}s echelon)", flush=True)
        return res["rank"]
    finally:
        if not keep and workdir is None:
            for p in (mat, mat + ".meta.json", out):
                try:
                    os.remove(p)
                except OSError:
                    pass
            try:
                os.rmdir(tmp)
            except OSError:
                pass
