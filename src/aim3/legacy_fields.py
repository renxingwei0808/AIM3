"""Extend `legacy/aim3_lab.CONWAY` to n > 16 without touching `legacy/` (CLAUDE.md rule 3).

`legacy/aim3_lab.F` looks its modulus up in a table that stops at n = 16, so every Tier-A script
raises `KeyError: 17` at the sizes this server batch needs.  `src/aim3/gf2n.TOY_MODULI` already
carries irreducible polynomials for n = 17..24 (and they agree with `CONWAY` for every n <= 16, so
the toy instances stay bit-identical to the Tier-A ones).  Importing this module copies the missing
entries across, after **verifying that each added polynomial really is irreducible** — a reducible
modulus would not raise, it would silently produce wrong ranks.

    import aim3.legacy_fields          # noqa: F401  (side effect: CONWAY gains n = 17..24)
"""
from __future__ import annotations

from .gf2n import TOY_MODULI

__all__ = ["extend_conway", "is_irreducible", "ADDED"]

ADDED = {}


def _polymulmod(a: int, b: int, mod: int, n: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> n & 1:
            a ^= mod
    return r


def _polymod(a: int, mod: int) -> int:
    dm = mod.bit_length() - 1
    while a.bit_length() - 1 >= dm and a:
        a ^= mod << (a.bit_length() - 1 - dm)
    return a


def _polygcd(a: int, b: int) -> int:
    while b:
        a, b = b, _polyrem(a, b)
    return a


def _polyrem(a: int, b: int) -> int:
    db = b.bit_length() - 1
    while a and a.bit_length() - 1 >= db:
        a ^= b << (a.bit_length() - 1 - db)
    return a


def is_irreducible(mod: int, n: int) -> bool:
    """Rabin's test for f of degree n over F2: x^(2^n) = x (mod f) and gcd(x^(2^(n/p)) - x, f) = 1."""
    if mod.bit_length() - 1 != n:
        return False
    if not mod & 1:                              # divisible by x
        return False

    def xpow2k(k):                               # x^(2^k) mod f, by k repeated squarings
        v = 2                                    # the polynomial x
        for _ in range(k):
            v = _polymulmod(v, v, mod, n)
        return v

    if xpow2k(n) != 2:
        return False
    m, primes = n, set()
    d = 2
    while d * d <= m:
        while m % d == 0:
            primes.add(d)
            m //= d
        d += 1
    if m > 1:
        primes.add(m)
    for p in primes:
        if _polygcd(mod, xpow2k(n // p) ^ 2) != 1:
            return False
    return True


def extend_conway(max_n: int = 24):
    """Add the missing moduli to `aim3_lab.CONWAY`; return {n: modulus} for what was added."""
    import aim3_lab                               # legacy module, on sys.path via experiments/_common
    added = {}
    for n, mod in sorted(TOY_MODULI.items()):
        if n > max_n or n in aim3_lab.CONWAY:
            continue
        if not is_irreducible(mod, n):
            raise ValueError(f"TOY_MODULI[{n}] = {mod:#x} is not irreducible over F2")
        aim3_lab.CONWAY[n] = mod
        added[n] = mod
    # the shared entries must agree, otherwise Tier-A and server runs describe different fields
    for n, mod in TOY_MODULI.items():
        if n in aim3_lab.CONWAY and aim3_lab.CONWAY[n] != mod and n not in added:
            raise ValueError(f"modulus mismatch at n={n}: legacy {aim3_lab.CONWAY[n]:#x} vs "
                             f"gf2n {mod:#x}")
    ADDED.update(added)
    return added


extend_conway()
