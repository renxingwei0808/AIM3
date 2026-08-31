"""Algebraic structures of the enhanced Mersenne S-box (Focus Report §2, structures A–G), for any n.

Notation: q=2^n, gcd(e,n)=1, r=(2^e-1)^{-1} mod (q-1), m=e^{-1} mod n.
  S_e(x)      = x^{2^e-1} + x^{-1} = (x+1)^{2^e}/x                       (A)
  Im S_e      = {0} u {y : Tr(y^{-r})=0},  |Im| = q/2                    (B)
  kernel var  a = y^r : B_e(a,x)=a, y a = a^{2^e}; preimages {x, x+a}   (D)
  B_e(u,v)    = u v^{2^e} + u^{2^e} v = u v (M_e(u)+M_e(v))              (polarised form)
  wt(r)=m, wt(-r)=n-m, supp(r) = {i e mod n : i<m}                       (E)
  two-valued inverse G_{e,eps}(y) = a (H_e(1/a) + eps) + 1, a=y^r,        (C)
      H_e = linear section of s -> s^{2^e}+s on the trace-zero hyperplane
"""
from __future__ import annotations
from math import gcd
from .gf2n import GF2n


def weight_data(n: int, e: int):
    q1 = (1 << n) - 1
    r = pow((1 << e) - 1, -1, q1)
    m = pow(e, -1, n)
    return dict(e=e, r=r, m=m, wt_r=bin(r).count("1"), wt_minus_r=bin((-r) % q1).count("1"),
                supp_is_AP=(sorted(i for i in range(n) if (r >> i) & 1) == sorted((i * e) % n for i in range(m))))


def degenerate_class(n: int, e: int, kmax: int = 4):
    """Frobenius-orbit coincidences that empirically inflate degree falls / QE dimensions (experiment E3/E8).
    Returns a list of tags; empty list = generic exponent."""
    tags = []
    if e == 2:
        tags.append("e=2")
    if n % 2 == 0 and e == 1 + n // 2:
        tags.append("e=1+n/2")
    for k in range(2, kmax + 1):
        for s in (1, -1, 2, -2):
            if (k * e - s) % n == 0:
                tags.append(f"{k}e={s:+d} mod n")
    return tags


class SboxStructures:
    def __init__(self, n: int, e: int, F: GF2n | None = None):
        assert gcd(e, n) == 1
        self.n, self.e = n, e
        self.F = F or GF2n(n)
        self.r = self.F.r_inv(e)
        self._H = None

    # --- linear section H_e of L(s) = s^{2^e} + s (image = trace-zero hyperplane) ---
    def _build_section(self):
        """solve L(s) = b by F2-linear algebra; H maps each basis vector of the trace-zero hyperplane to a preimage."""
        F, n, e = self.F, self.n, self.e
        # matrix of L in the polynomial basis: column j = L(X^j)
        cols = [F.frob(1 << j, e) ^ (1 << j) for j in range(n)]
        # Gaussian elimination to express arbitrary b (with Tr b = 0) as L(s): precompute pivots
        piv = {}   # leading bit -> (value, combination of basis indices as bitmask)
        for j, c in enumerate(cols):
            comb = 1 << j
            while c:
                lead = c.bit_length() - 1
                if lead in piv:
                    c ^= piv[lead][0]; comb ^= piv[lead][1]
                else:
                    piv[lead] = (c, comb); break
        self._H = piv

    def H(self, b: int) -> int:
        """a preimage s of b under s -> s^{2^e}+s; requires Tr(b)=0."""
        if self._H is None:
            self._build_section()
        comb = 0
        while b:
            lead = b.bit_length() - 1
            if lead not in self._H:
                raise ValueError("b not in image (Tr(b) != 0)")
            b ^= self._H[lead][0]; comb ^= self._H[lead][1]
        return comb   # bitmask over basis = field element in polynomial basis

    def kernel(self, y: int) -> int:
        return self.F.pow(y, self.r)

    def in_image(self, y: int) -> bool:
        return y == 0 or self.F.trace(self.F.pow(y, -self.r)) == 0

    def two_valued_inverse(self, y: int):
        """both preimages of y != 0 (raises if y not in image); y=0 -> (1,)"""
        F = self.F
        if y == 0:
            return (1,)
        a = self.kernel(y)
        sigma = self.H(F.inv(a))                  # sigma^{2^e}+sigma = 1/a
        x0 = F.mul(a, sigma) ^ 1
        return (x0, x0 ^ a)

    def B(self, u: int, v: int) -> int:
        F, e = self.F, self.e
        return F.mul(u, F.frob(v, e)) ^ F.mul(F.frob(u, e), v)

    def check_random(self, rng, trials: int = 50) -> bool:
        F, e = self.F, self.e
        for _ in range(trials):
            x = rng.randrange(1, F.q)
            y = F.sbox(x, e)
            assert y == F.mul(F.frob(x ^ 1, e), F.inv(x))                       # (A)
            assert self.in_image(y)                                            # (B)
            if y:
                a = self.kernel(y)
                assert self.B(a, x) == a and F.mul(y, a) == F.frob(a, e)        # (D)
                assert self.B(a, x) == F.mul(F.mul(a, x), F.mersenne(a, e) ^ F.mersenne(x, e))  # polarised form
                assert set(self.two_valued_inverse(y)) == {x, x ^ a}            # (C)
        return True


if __name__ == "__main__":
    import random
    rng = random.Random(0)
    for n, exps in [(128, (3, 7, 11)), (192, (11, 23, 47)), (256, (3, 7, 11, 19))]:
        for e in exps:
            S = SboxStructures(n, e)
            S.check_random(rng, 5)
            print(f"n={n} e={e:2d}: structures A/B/C/D OK; {weight_data(n, e)} ; degenerate tags={degenerate_class(n, e)}")
