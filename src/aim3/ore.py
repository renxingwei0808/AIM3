"""Ore (skew) polynomials over GF(2^n): R = F_{2^n}<phi> with phi*a = a^2*phi.

An element f = sum_i c_i phi^i acts on F_{2^n} as the linearised polynomial
    f(v) = sum_i c_i v^{2^i},
and the ring multiplication is exactly composition of those maps:
    (a phi^i)(b phi^j) = a * b^{2^i} * phi^{i+j}.
`(f*g).apply(v) == f.apply(g.apply(v))` is asserted by tests/test_ore.py.

R is a left and right Euclidean domain, which is what makes track T4 (ROADMAP section 3.4) possible
at all: `divmod_right` divides on the right, `gcrd` is the greatest common RIGHT divisor and `lclm`
the least common LEFT multiple, with deg lclm = deg f + deg g - deg gcrd.

Two things to keep straight, because they are different objects:

* **The ring.**  phi^n is NOT 1 in R -- R is a domain and degrees add.  All of divmod_right/gcrd/lclm
  work here.
* **The action.**  On F_{2^n} we have v^{2^n} = v, so R only acts through R/(phi^n - 1), which is not
  a domain.  `apply` and `root_space` live here.  A skew polynomial of degree >= n therefore acts
  like a reduced one, and the root space of a degree-d element has dimension at most min(d, n)
  -- not d.  Any statement that mixes the two must say which one it means.

Pure standard library.
"""
from __future__ import annotations
from .gf2n import GF2n

__all__ = ["SkewPoly", "solve_gf2n", "linearise_values", "from_bit_matrix"]


# ---------------------------------------------------------------- GF(2^n) linear algebra
def solve_gf2n(F: GF2n, M, v):
    """Solve M h = v over GF(2^n) by Gaussian elimination (M a list of rows, v a list)."""
    n = len(M)
    A = [list(row) + [v[i]] for i, row in enumerate(M)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c]), None)
        if p is None:
            raise ValueError("singular matrix")
        A[c], A[p] = A[p], A[c]
        ic = F.inv(A[c][c])
        A[c] = [F.mul(ic, x) for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [x ^ F.mul(f, y) for x, y in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]


def linearise_values(F: GF2n, basis, values):
    """Coefficients h_j of sum_j h_j u^{2^j} agreeing with basis[k] |-> values[k].

    The matrix M[k][j] = basis[k]^{2^j} is the Moore matrix of the basis; it is invertible exactly
    when the basis is an F_2-basis of F_{2^n}.
    """
    n = F.n
    M = [[F.frob(b, j) for j in range(n)] for b in basis]
    return solve_gf2n(F, M, list(values))


def from_bit_matrix(F: GF2n, mat):
    """SkewPoly of the F_2-linear map given by a bit matrix in the `aim3_lab.apply_mat` convention
    (y_i = XOR_j mat[i][j] x_j)."""
    n = F.n

    def apply_mat(x):
        y = 0
        for i in range(n):
            s = 0
            row = mat[i]
            for j in range(n):
                if row[j] and (x >> j) & 1:
                    s ^= 1
            y |= s << i
        return y

    basis = [1 << k for k in range(n)]
    return SkewPoly(F, linearise_values(F, basis, [apply_mat(b) for b in basis]))


# ---------------------------------------------------------------- the ring
class SkewPoly:
    """sum_i coeffs[i] * phi^i over GF(2^n); coeffs[i] is an int in the gf2n representation."""
    __slots__ = ("F", "c")

    def __init__(self, F: GF2n, coeffs=()):
        self.F = F
        c = list(coeffs)
        while c and c[-1] == 0:
            c.pop()
        self.c = c

    # ---- basics ---------------------------------------------------------
    @classmethod
    def monomial(cls, F, coeff, k):
        return cls(F, [0] * k + [coeff])

    @classmethod
    def one(cls, F):
        return cls(F, [1])

    @classmethod
    def zero(cls, F):
        return cls(F, [])

    def is_zero(self):
        return not self.c

    def degree(self):
        """-1 for the zero polynomial."""
        return len(self.c) - 1

    def lead(self):
        return self.c[-1] if self.c else 0

    def __eq__(self, other):
        return isinstance(other, SkewPoly) and self.c == other.c and self.F.n == other.F.n

    def __repr__(self):
        if not self.c:
            return "SkewPoly(0)"
        return "SkewPoly(" + " + ".join(f"{c:x}*phi^{i}" for i, c in enumerate(self.c) if c) + ")"

    # ---- arithmetic -----------------------------------------------------
    def __add__(self, other):
        a, b = self.c, other.c
        m = max(len(a), len(b))
        return SkewPoly(self.F, [(a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0)
                                 for i in range(m)])

    __sub__ = __add__          # characteristic 2

    def __mul__(self, other):
        """(a phi^i)(b phi^j) = a * b^{2^i} * phi^{i+j}."""
        F = self.F
        if not self.c or not other.c:
            return SkewPoly(F, [])
        out = [0] * (len(self.c) + len(other.c) - 1)
        for i, a in enumerate(self.c):
            if not a:
                continue
            for j, b in enumerate(other.c):
                if b:
                    out[i + j] ^= F.mul(a, F.frob(b, i))
        return SkewPoly(F, out)

    def scale(self, c):
        """left multiplication by the scalar c."""
        return SkewPoly(self.F, [self.F.mul(c, x) for x in self.c])

    def monic(self):
        """left-scale so that the leading coefficient is 1 (the leading coefficient of c*f is c*lead f)."""
        if self.is_zero():
            return self
        return self.scale(self.F.inv(self.lead()))

    # ---- action ---------------------------------------------------------
    def apply(self, v):
        """sum_i c_i v^{2^i} -- the linearised polynomial, i.e. the action of R on F_{2^n}."""
        F = self.F
        s = 0
        for i, c in enumerate(self.c):
            if c:
                s ^= F.mul(c, F.frob(v, i))
        return s

    def root_space(self):
        """basis of {v in F_{2^n} : f(v) = 0}, by F_2-linear algebra on the action (not on the ring)."""
        F = self.F
        piv, ker = {}, []
        # image of each basis vector, tracked with its preimage
        rows = []
        for k in range(F.n):
            rows.append((self.apply(1 << k), 1 << k))
        for img, pre in rows:
            v, p = img, pre
            while v:
                b = v.bit_length() - 1
                if b in piv:
                    v ^= piv[b][0]
                    p ^= piv[b][1]
                else:
                    piv[b] = (v, p)
                    break
            else:
                ker.append(p)
        return ker

    # ---- Euclidean structure --------------------------------------------
    def divmod_right(self, g):
        """(q, r) with self = q*g + r and deg r < deg g.  Division on the RIGHT by g."""
        F = self.F
        if g.is_zero():
            raise ZeroDivisionError("divmod_right by zero")
        r = SkewPoly(F, self.c)
        q = [0] * max(1, self.degree() - g.degree() + 1)
        dg, lg = g.degree(), g.lead()
        while not r.is_zero() and r.degree() >= dg:
            k = r.degree() - dg
            # need c with c * lg^{2^k} = lead(r)
            c = F.mul(r.lead(), F.inv(F.frob(lg, k)))
            q[k] ^= c
            r = r - SkewPoly.monomial(F, c, k) * g
        return SkewPoly(F, q), r

    def gcrd(self, other):
        """greatest common right divisor, returned monic (zero only if both are zero)."""
        a, b = SkewPoly(self.F, self.c), SkewPoly(self.F, other.c)
        while not b.is_zero():
            a, b = b, a.divmod_right(b)[1]
        return a.monic()

    def lclm(self, other):
        """least common left multiple, returned monic.

        Extended right Euclid: keeping r_k = u_k*self + v_k*other, the step after the gcrd gives
        r = 0, i.e. u*self = v*other, and that common value is the LCLM (characteristic 2, so the
        signs vanish).
        """
        F = self.F
        if self.is_zero() or other.is_zero():
            return SkewPoly(F, [])
        r0, r1 = SkewPoly(F, self.c), SkewPoly(F, other.c)
        u0, u1 = SkewPoly.one(F), SkewPoly.zero(F)
        while not r1.is_zero():
            q, r = r0.divmod_right(r1)
            r0, r1 = r1, r
            u0, u1 = u1, u0 - q * u1
        return (u1 * self).monic()
