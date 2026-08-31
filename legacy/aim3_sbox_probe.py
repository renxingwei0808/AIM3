#!/usr/bin/env python3
"""
Small-field probe for the AIM3 enhanced Mersenne S-box.

For S_e(x) = x^(2^e-1) + x^(-1), x != 0, this script verifies:
  1. the image has exactly 2^(n-1) elements;
  2. y=0 has the unique nonzero preimage x=1;
  3. each nonzero image value has exactly two preimages;
  4. y != 0 is in the image iff Tr(y^(-r)) = 0, where
       r = (2^e - 1)^(-1) mod (2^n - 1);
  5. the two-valued inverse recovered from a linear section of
       u^(2^e) + u = b
     matches exhaustive enumeration.

This is intended for toy parameters. Exhaustive enumeration is exponential in n.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from math import gcd
from typing import Dict, Iterable, List, Set


@dataclass(frozen=True)
class GF2N:
    n: int
    modulus: int  # Includes the x^n term.

    def __post_init__(self) -> None:
        if self.n <= 0:
            raise ValueError("n must be positive")
        if (self.modulus >> self.n) != 1:
            raise ValueError(
                "modulus must be a degree-n binary polynomial including x^n"
            )
        if (self.modulus & 1) == 0:
            raise ValueError("modulus must have nonzero constant term")

    @property
    def q(self) -> int:
        return 1 << self.n

    @property
    def mask(self) -> int:
        return self.q - 1

    def add(self, a: int, b: int) -> int:
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        a &= self.mask
        b &= self.mask
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.modulus
        return result & self.mask

    def square(self, a: int) -> int:
        return self.mul(a, a)

    def pow(self, a: int, exponent: int) -> int:
        if exponent < 0:
            if a == 0:
                raise ZeroDivisionError("negative power of zero")
            a = self.inv(a)
            exponent = -exponent
        result = 1
        base = a & self.mask
        while exponent:
            if exponent & 1:
                result = self.mul(result, base)
            base = self.square(base)
            exponent >>= 1
        return result

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("zero has no inverse")
        return self.pow(a, self.q - 2)

    def trace(self, a: int) -> int:
        value = 0
        current = a & self.mask
        for _ in range(self.n):
            value ^= current
            current = self.square(current)
        if value not in (0, 1):
            raise ArithmeticError(
                f"trace did not land in GF(2); check irreducibility: {value:#x}"
            )
        return value


def sbox(field: GF2N, x: int, e: int) -> int:
    """AIM3 enhanced Mersenne S-box, including S_e(0)=0."""
    if x == 0:
        return 0
    return field.add(
        field.pow(x, (1 << e) - 1),
        field.inv(x),
    )


def build_section(field: GF2N, e: int) -> Dict[int, int]:
    """
    Choose one representative u for every b in the image of u -> u^(2^e)+u.
    The image is the trace-zero hyperplane when gcd(e,n)=1.
    """
    section: Dict[int, int] = {}
    for u in range(field.q):
        b = field.add(field.pow(u, 1 << e), u)
        section.setdefault(b, u)
    return section


def analyze(field: GF2N, e: int) -> Dict[str, int]:
    if gcd(e, field.n) != 1:
        raise ValueError("AIM3 requires gcd(e,n)=1")

    q_minus_1 = field.q - 1
    d = (1 << e) - 1
    r = pow(d, -1, q_minus_1)

    fibers: Dict[int, List[int]] = defaultdict(list)
    for x in range(1, field.q):
        fibers[sbox(field, x, e)].append(x)

    image: Set[int] = set(fibers)
    counts = Counter(len(values) for values in fibers.values())

    if len(image) != field.q // 2:
        raise AssertionError(
            f"wrong image size: got {len(image)}, expected {field.q // 2}"
        )
    if fibers.get(0) != [1]:
        raise AssertionError(f"zero-output fiber should be [1], got {fibers.get(0)}")
    for y, xs in fibers.items():
        if y != 0 and len(xs) != 2:
            raise AssertionError(f"nonzero image value {y:#x} has {len(xs)} preimages")

    trace_mismatches = 0
    for y in range(1, field.q):
        minus_r = (-r) % q_minus_1
        predicted = field.trace(field.pow(y, minus_r)) == 0
        actual = y in image
        if predicted != actual:
            trace_mismatches += 1

    if trace_mismatches:
        raise AssertionError(f"trace criterion mismatches: {trace_mismatches}")

    section = build_section(field, e)
    inverse_mismatches = 0
    for y in image:
        if y == 0:
            continue
        a = field.pow(y, r)
        b = field.inv(a)  # y^(-r)
        if field.trace(b) != 0:
            raise AssertionError("image value produced nonzero trace")
        u = section[b]

        base = field.mul(a, field.pow(u, 1 << e))
        predicted_preimages = {base, field.add(base, a)}
        actual_preimages = set(fibers[y])
        if predicted_preimages != actual_preimages:
            inverse_mismatches += 1

    if inverse_mismatches:
        raise AssertionError(
            f"two-valued inverse mismatches: {inverse_mismatches}"
        )

    result = {
        "n": field.n,
        "e": e,
        "q": field.q,
        "image_size": len(image),
        "zero_fiber_size": len(fibers[0]),
        "nonzero_image_values": len(image) - 1,
        "trace_mismatches": trace_mismatches,
        "inverse_mismatches": inverse_mismatches,
        "weight_r": r.bit_count(),
        "weight_minus_r": ((-r) % q_minus_1).bit_count(),
    }

    print(
        f"GF(2^{field.n}), e={e}: "
        f"|Im|={result['image_size']}={field.q // 2}, "
        f"fiber histogram={dict(sorted(counts.items()))}, "
        f"trace mismatches={trace_mismatches}, "
        f"inverse mismatches={inverse_mismatches}, "
        f"wt(r)={result['weight_r']}, "
        f"wt(-r)={result['weight_minus_r']}"
    )
    return result


DEFAULT_CASES = [
    # x^5 + x^2 + 1
    (5, 2, 0x25),
    # x^7 + x + 1
    (7, 3, 0x83),
    # AES polynomial x^8 + x^4 + x^3 + x + 1
    (8, 3, 0x11B),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Exhaustively probe the AIM3 enhanced Mersenne S-box."
    )
    parser.add_argument("--n", type=int, help="extension degree")
    parser.add_argument("--e", type=int, help="Mersenne exponent parameter")
    parser.add_argument(
        "--modulus",
        type=lambda value: int(value, 0),
        help="irreducible polynomial as an integer, including x^n (e.g. 0x11b)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    supplied = [args.n is not None, args.e is not None, args.modulus is not None]
    if any(supplied) and not all(supplied):
        raise SystemExit("--n, --e, and --modulus must be supplied together")

    cases: Iterable[tuple[int, int, int]]
    if all(supplied):
        cases = [(args.n, args.e, args.modulus)]
    else:
        cases = DEFAULT_CASES

    for n, e, modulus in cases:
        analyze(GF2N(n=n, modulus=modulus), e=e)


if __name__ == "__main__":
    main()
