"""Fast Boolean Macaulay-matrix statistics (pure standard library).

Same semantics as `legacy/aim3_macaulay.macaulay`, but

* monomials are integer bitmasks over the variables instead of `frozenset`s;
* rows are streamed straight into the elimination, so the matrix is never stored;
* the full rank and the rank of the top-degree projection come out of ONE pass.

The one-pass trick: `monomial_columns` numbers the monomials by increasing degree, so the
degree-D columns occupy the *highest* bit positions and `bit_length()`-driven elimination
pivots on them first.  In the resulting echelon form the pivots lying in the low region are
exactly a basis of rowspace ∩ (degree ≤ D−1), hence

    dim(span ∩ deg≤D−1) = rank − rank(top-degree projection) = #pivots below `top_start`.

`tests/test_core.py::test_boolmac_matches_legacy` checks this against the legacy two-pass
implementation on small systems; do not change the column ordering without re-running it.
"""
from __future__ import annotations
import itertools, time

__all__ = ["monomial_columns", "poly_rows", "multiplier_masks", "labelled_rows", "macaulay_stats", "fall_space", "span_rank", "reduce_modulo",
           "MacaulayResult", "MacaulayTimeout"]


class MacaulayTimeout(Exception):
    """Raised when a run exceeds its deadline. `.partial` is a MacaulayResult whose rank/dim_low are
    LOWER BOUNDS (rows counted so far), never values -- CLAUDE.md ground rule 4."""

    def __init__(self, partial, seconds):
        super().__init__(f"Macaulay elimination exceeded {seconds:.0f}s; partial (lower bound): {partial}")
        self.partial = partial
        self.seconds = seconds


class MacaulayResult(tuple):
    """(rows, cols, rank, syz, dim_low) with attribute access."""
    __slots__ = ()

    def __new__(cls, rows, cols, rank, dim_low):
        return super().__new__(cls, (rows, cols, rank, rows - rank, dim_low))

    rows = property(lambda s: s[0])
    cols = property(lambda s: s[1])
    rank = property(lambda s: s[2])
    syz = property(lambda s: s[3])
    dim_low = property(lambda s: s[4])

    def __repr__(self):
        return (f"MacaulayResult(rows={self.rows}, cols={self.cols}, rank={self.rank}, "
                f"syz={self.syz}, dim_low={self.dim_low})")


def monomial_columns(nv: int, D: int):
    """{variable-bitmask: column index} ordered by increasing degree, plus (ncol, top_start).

    Columns [top_start, ncol) are the degree-D monomials.
    """
    col = {}
    c = 0
    top_start = 0
    for d in range(D + 1):
        if d == D:
            top_start = c
        for m in itertools.combinations(range(nv), d):
            mask = 0
            for i in m:
                mask |= 1 << i
            col[mask] = c
            c += 1
    return col, c, top_start


def _as_masks(poly):
    """legacy polynomial (set of frozensets of variable indices) -> list of variable bitmasks."""
    out = []
    for m in poly:
        mask = 0
        for i in m:
            mask |= 1 << i
        out.append(mask)
    return out


def multiplier_masks(nv: int, mult_deg: int, mult_vars=None):
    """The multiplier monomials, as variable bitmasks, in the canonical order.

    `mult_vars` restricts the multipliers to squarefree monomials in that subset of the variables
    (the identity multiplier is always present).  It is what a *block-restricted* XL uses: for the
    E06 kernel model, multiplying only by the W block asks whether the bilinear structure lets the
    system linearise at a lower degree than the full multiplier set does.  `None` = every variable,
    which is the historical behaviour.
    """
    pool = range(nv) if mult_vars is None else sorted(set(mult_vars))
    out = []
    for d in range(mult_deg + 1):
        for m in itertools.combinations(pool, d):
            mask = 0
            for i in m:
                mask |= 1 << i
            out.append(mask)
    return out


def poly_rows(eqs, nv: int, D: int, col=None, mult_deg=None, mult_vars=None):
    """Yield the Macaulay rows of `eqs` at degree D as integers (bit i = column i).

    Multipliers are all squarefree monomials of degree <= `mult_deg` (default D-2, the
    convention used by legacy/aim3_macaulay.py: an equation of degree 2 times a multiplier of
    degree <= D-2 lands in degree <= D), optionally restricted to the variables `mult_vars`.
    Rows that vanish or exceed degree D are dropped.
    """
    if col is None:
        col, _, _ = monomial_columns(nv, D)
    if mult_deg is None:
        mult_deg = D - 2
    mults = multiplier_masks(nv, mult_deg, mult_vars)
    for f in eqs:
        fm = _as_masks(f)
        for mu in mults:
            acc = set()
            for vm in fm:
                acc ^= {vm | mu}
            if not acc:
                continue
            row = 0
            over = False
            for key in acc:
                if key.bit_count() > D:
                    over = True
                    break
                row |= 1 << col[key]
            if over or row == 0:
                continue
            yield row


def _eliminate(row_iter, top_start, ncol=None, deadline=None):
    """Streamed GF(2) echelon: returns (nrows_seen, rank, #pivots below top_start).

    Two shortcuts that do not change the result:
    * once every column is a pivot the matrix has full rank, so the remaining rows only need
      counting (they are all dependent);
    * `deadline` (an absolute time.monotonic() value) aborts with MacaulayTimeout carrying the
      partial counts as lower bounds.
    """
    piv = {}
    nrows = 0
    low = 0
    full = False
    check = 0
    for r in row_iter:
        nrows += 1
        if full:
            continue
        if deadline is not None:
            check += 1
            if check >= 256:
                check = 0
                if time.monotonic() > deadline:
                    raise MacaulayTimeout(MacaulayResult(nrows, ncol if ncol else 0, len(piv), low),
                                          time.monotonic() - deadline)
        while r:
            b = r.bit_length() - 1
            p = piv.get(b)
            if p is None:
                piv[b] = r
                if b < top_start:
                    low += 1
                if ncol is not None and len(piv) == ncol:
                    full = True
                break
            r ^= p
    return nrows, len(piv), low


def macaulay_stats(eqs, nv: int, D: int, mult_deg=None, budget=None,
                   mult_vars=None) -> MacaulayResult:
    """Rows / columns / rank / syzygies / dim(span ∩ deg<=D-1) of the degree-D Macaulay matrix.

    `budget` (seconds) raises MacaulayTimeout instead of running forever; the exception carries
    the partial counts as lower bounds.
    """
    col, ncol, top_start = monomial_columns(nv, D)
    deadline = None if budget is None else time.monotonic() + budget
    nrows, rank, low = _eliminate(poly_rows(eqs, nv, D, col=col, mult_deg=mult_deg,
                                            mult_vars=mult_vars),
                                  top_start, ncol=ncol, deadline=deadline)
    return MacaulayResult(nrows, ncol, rank, low)


def labelled_rows(eqs, nv: int, D: int, col=None, mult_deg=None):
    """Like `poly_rows`, but yields (row, eq_index, multiplier_mask) so a degree fall can be traced
    back to the multiplier x equation combinations that produced it."""
    if col is None:
        col, _, _ = monomial_columns(nv, D)
    if mult_deg is None:
        mult_deg = D - 2
    mults = []
    for d in range(mult_deg + 1):
        for m in itertools.combinations(range(nv), d):
            mask = 0
            for i in m:
                mask |= 1 << i
            mults.append(mask)
    for ei, f in enumerate(eqs):
        fm = _as_masks(f)
        for mu in mults:
            acc = set()
            for vm in fm:
                acc ^= {vm | mu}
            if not acc:
                continue
            row = 0
            over = False
            for key in acc:
                if key.bit_count() > D:
                    over = True
                    break
                row |= 1 << col[key]
            if over or row == 0:
                continue
            yield row, ei, mu


def fall_space(eqs, nv: int, D: int, mult_deg=None, track=False):
    """Basis of rowspace(Mac_D) n {degree <= D-1}, i.e. the elements whose top-degree part cancels.

    Eliminates using ONLY the degree-D columns as pivots (they are the highest bit positions, see
    `monomial_columns`); a row that reduces to zero there has fallen, and what is left is its
    degree-<=D-1 part.  Returns a dict with

        low        list of ints (bit i = column i), spanning the fall space
        combos     with track=True, the matching list of combination masks over the labelled rows,
                   so combos[k] bit r means "labelled row r took part in fall k"
        labels     with track=True, the (eq_index, multiplier_mask) of labelled row r
        rank_top   rank of the degree-D projection
        cols       total number of columns
        top_start  first degree-D column index

    dim(span(low)) equals `macaulay_stats(...).dim_low`; the test in tests/test_core.py checks it.
    """
    col, ncol, top_start = monomial_columns(nv, D)
    topmask = ((1 << ncol) - 1) ^ ((1 << top_start) - 1)
    piv = {}                      # top-region leading bit -> (row, combo)
    low, combos, labels = [], [], []
    src = (labelled_rows(eqs, nv, D, col=col, mult_deg=mult_deg) if track
           else ((r, None, None) for r in poly_rows(eqs, nv, D, col=col, mult_deg=mult_deg)))
    for idx, (r, ei, mu) in enumerate(src):
        if track:
            labels.append((ei, mu))
        c = 1 << idx if track else 0
        top = r & topmask
        while top:
            b = top.bit_length() - 1
            p = piv.get(b)
            if p is None:
                piv[b] = (r, c)
                break
            r ^= p[0]
            c ^= p[1]
            top = r & topmask
        else:
            if r:
                low.append(r)
                if track:
                    combos.append(c)
    return dict(low=low, combos=combos if track else None, labels=labels if track else None,
                rank_top=len(piv), cols=ncol, top_start=top_start)


def span_rank(rows):
    """rank of a list of ints over GF(2)."""
    piv = {}
    for r in rows:
        while r:
            b = r.bit_length() - 1
            if b in piv:
                r ^= piv[b]
            else:
                piv[b] = r
                break
    return len(piv)


def reduce_modulo(rows, basis):
    """Return (independent representatives of `rows` modulo span(basis), quotient dimension)."""
    piv = {}
    for r in basis:
        while r:
            b = r.bit_length() - 1
            if b in piv:
                r ^= piv[b]
            else:
                piv[b] = r
                break
    reps = []
    for r in rows:
        while r:
            b = r.bit_length() - 1
            if b in piv:
                r ^= piv[b]
            else:
                piv[b] = r
                reps.append(r)
                break
    return reps, len(reps)
