"""Field-valued Boolean polynomials: F = sum_m c_m * m, c_m in GF(2^n), m a squarefree monomial in
the Boolean variables.

This is the representation that makes the *field* identities behind the AIM3 degree falls visible.
`legacy/aim3_macaulay.py` computes the same objects but immediately splits them into n Boolean
coordinate polynomials; here the field coefficient is kept, so that

    Tr_gamma(F) = sum_m Tr(gamma * c_m) * m

gives every Boolean coordinate of F at once, and the n-dimensional F_2-space
{Tr(gamma F) : gamma in GF(2^n)} is exactly the span of F's coordinate polynomials.

Why the Boolean-monomial keys are the right thing: a field power x^{2^k} is F_2-linear in the bits
of x, so x^{2^k} = sum_j beta_j^{2^k} x_j is a FieldPoly with n degree-1 terms; multiplying two such
FieldPolys unions the monomials, and m u m = m implements the Boolean relation x_j^2 = x_j for free.
So `X.frob(i) * X` really is the field element x^{2^i+1}, with the correct Boolean degree
wt(2^i + 1) rather than 2.

Everything here is pure standard library.
"""
from __future__ import annotations
from .gf2n import GF2n

__all__ = ["FieldPoly", "var_field", "trace_poly", "boolean_coords"]


class FieldPoly:
    """sum_m c_m * m  with m a frozenset of Boolean variable indices and c_m in GF(2^n)."""
    __slots__ = ("F", "t")

    def __init__(self, F: GF2n, terms=None):
        self.F = F
        self.t = {}
        if terms:
            for m, c in terms.items():
                self.add(m, c)

    # ---- construction -------------------------------------------------------
    def add(self, mono, c):
        if not c:
            return
        mono = frozenset(mono)
        v = self.t.get(mono, 0) ^ c
        if v:
            self.t[mono] = v
        else:
            self.t.pop(mono, None)

    @classmethod
    def constant(cls, F, c):
        return cls(F, {frozenset(): c})

    # ---- arithmetic ---------------------------------------------------------
    def __add__(self, other):
        r = FieldPoly(self.F, self.t)
        for m, c in other.t.items():
            r.add(m, c)
        return r

    def __mul__(self, other):
        r = FieldPoly(self.F)
        mul = self.F.mul
        for m1, c1 in self.t.items():
            for m2, c2 in other.t.items():
                r.add(m1 | m2, mul(c1, c2))
        return r

    def scale(self, c):
        return FieldPoly(self.F, {m: self.F.mul(c, cm) for m, cm in self.t.items()})

    def frob(self, k: int):
        """F^{2^k}: squaring is F_2-linear and monomials are idempotent, so only the coefficients move."""
        k %= self.F.n
        return FieldPoly(self.F, {m: self.F.frob(c, k) for m, c in self.t.items()})

    # ---- inspection ---------------------------------------------------------
    def degree(self):
        return max((len(m) for m in self.t), default=-1)

    def part(self, d):
        """the homogeneous Boolean-degree-d part."""
        return FieldPoly(self.F, {m: c for m, c in self.t.items() if len(m) == d})

    def is_zero(self):
        return not self.t

    def __repr__(self):
        return f"FieldPoly(n={self.F.n}, {len(self.t)} terms, degree {self.degree()})"


def var_field(F: GF2n, offset: int, basis=None) -> FieldPoly:
    """the field variable whose bit j is Boolean variable `offset + j`, in the polynomial basis."""
    basis = basis or [1 << j for j in range(F.n)]
    return FieldPoly(F, {frozenset((offset + j,)): basis[j] for j in range(F.n)})


def trace_poly(P: FieldPoly, gamma: int):
    """Tr(gamma * P) as a Boolean polynomial: a set of monomials (frozensets), XOR semantics."""
    F = P.F
    out = set()
    for m, c in P.t.items():
        if F.trace(F.mul(gamma, c)):
            out ^= {m}
    return out


def boolean_coords(P: FieldPoly, basis=None):
    """the n Boolean coordinate polynomials of P (coordinate i = coefficient of X^i), in the same
    `set of frozensets` format as legacy/aim3_macaulay.py."""
    n = P.F.n
    out = [set() for _ in range(n)]
    for m, c in P.t.items():
        for i in range(n):
            if (c >> i) & 1:
                out[i] ^= {m}
    return out
