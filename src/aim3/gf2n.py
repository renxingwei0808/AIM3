"""GF(2^n) arithmetic in the *spec* representation: an integer whose bit i is the coefficient of X^i
(AIMer v3.0 §4.1.1). Works for any n (carry-less multiply + reduction); for n<=16 log/antilog tables are built.
"""
from __future__ import annotations

# AIMer v3.0 §3.1 irreducible polynomials
SPEC_MODULI = {
    128: (1 << 128) | (1 << 7) | (1 << 2) | (1 << 1) | 1,
    192: (1 << 192) | (1 << 7) | (1 << 2) | (1 << 1) | 1,
    256: (1 << 256) | (1 << 10) | (1 << 5) | (1 << 2) | 1,
}
# irreducible polynomials for toy fields (n<=16 identical to legacy/aim3_lab.py so results are comparable)
TOY_MODULI = {
    3: 0xB, 4: 0x13, 5: 0x25, 6: 0x43, 7: 0x83, 8: 0x11B, 9: 0x211, 10: 0x409,
    11: 0x805, 12: 0x1053, 13: 0x201B, 14: 0x402B, 15: 0x8003, 16: 0x1002D,
    17: 0x20009, 18: 0x40081, 19: 0x80027, 20: 0x100009, 21: 0x200005, 22: 0x400003,
    23: 0x800021, 24: 0x100001B, 28: 0x10000009, 32: 0x1000000AF,
}


class GF2n:
    def __init__(self, n: int, modulus: int | None = None, tables: bool | None = None):
        self.n = n
        self.mod = modulus or SPEC_MODULI.get(n) or TOY_MODULI.get(n)
        if self.mod is None:
            raise ValueError(f"no modulus known for n={n}; pass one explicitly")
        assert self.mod >> n == 1, "modulus must have degree exactly n"
        self.q = 1 << n
        self.mask = self.q - 1
        self._tab = None
        if tables is None:
            tables = n <= 16
        if tables:
            self._build_tables()

    # ---- core arithmetic -------------------------------------------------
    def mul_slow(self, a: int, b: int) -> int:
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.mod
        return r & self.mask

    def _build_tables(self):
        g = self._find_generator()
        exp = [0] * (2 * self.q)
        log = [0] * self.q
        v = 1
        for i in range(self.q - 1):
            exp[i] = v
            log[v] = i
            v = self.mul_slow(v, g)
        for i in range(self.q - 1, 2 * self.q):
            exp[i] = exp[i - (self.q - 1)]
        self._tab = (exp, log)

    def _find_generator(self):
        order = self.q - 1
        facs, m, d = set(), order, 2
        while d * d <= m:
            while m % d == 0:
                facs.add(d); m //= d
            d += 1
        if m > 1:
            facs.add(m)
        for g in range(2, self.q):
            if all(self.pow_slow(g, order // p) != 1 for p in facs):
                return g
        raise RuntimeError("no generator found")

    def pow_slow(self, a: int, k: int) -> int:
        r = 1
        while k:
            if k & 1:
                r = self.mul_slow(r, a)
            a = self.mul_slow(a, a)
            k >>= 1
        return r

    def mul(self, a: int, b: int) -> int:
        if self._tab is None:
            return self.mul_slow(a, b)
        if a == 0 or b == 0:
            return 0
        exp, log = self._tab
        return exp[log[a] + log[b]]

    def sq(self, a: int) -> int:
        return self.mul(a, a)

    def pow(self, a: int, k: int) -> int:
        """a^k, k may be negative (exponent taken mod q-1)."""
        if a == 0:
            if k > 0:
                return 0
            if k == 0:
                return 1
            raise ZeroDivisionError("negative power of 0")
        k %= (self.q - 1)
        if self._tab is not None:
            exp, log = self._tab
            return exp[(log[a] * k) % (self.q - 1)]
        return self.pow_slow(a, k)

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError
        return self.pow(a, self.q - 2)

    def frob(self, a: int, k: int = 1) -> int:
        """a^(2^k)"""
        for _ in range(k % self.n):
            a = self.sq(a)
        return a

    def trace(self, a: int) -> int:
        t, c = 0, a
        for _ in range(self.n):
            t ^= c
            c = self.sq(c)
        assert t in (0, 1)
        return t

    # ---- vector view (spec §4.1.1: bit i <-> coefficient of X^i <-> column-vector entry i) ----
    def to_vec(self, a: int) -> list[int]:
        return [(a >> i) & 1 for i in range(self.n)]

    def from_vec(self, v) -> int:
        return sum((bit & 1) << i for i, bit in enumerate(v))

    # ---- AIM3 building blocks ---------------------------------------------
    def mersenne(self, x: int, e: int) -> int:
        return self.pow(x, (1 << e) - 1) if x else 0

    def sbox(self, x: int, e: int) -> int:
        """enhanced Mersenne S-box S[e](x) = x^(2^e-1) + x^(2^n-2), S(0)=0."""
        if x == 0:
            return 0
        return self.pow(x, (1 << e) - 1) ^ self.inv(x)

    def r_inv(self, e: int) -> int:
        """r = (2^e-1)^{-1} mod 2^n-1 (exists iff gcd(e,n)=1)."""
        return pow((1 << e) - 1, -1, self.q - 1)
