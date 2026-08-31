"""Spec-faithful AIM3 (AIMer v3.0, Algorithms 4 and 5) in pure Python, for any n (recommended or toy).

    AIM3(pt, iv) = S[e_l] o Lin[iv] o S[e_0..e_{l-1}] o Lin_in[iv](pt) (+) pt

ENCODING CONVENTIONS THAT THE SPEC TEXT DOES NOT FIX UNAMBIGUOUSLY are isolated in `Conventions`
and MUST be pinned against the official KAT / reference C code in experiment E0 before any
full-size number derived from this file is quoted in a paper:
  * tape_bit_order    : how bit k of the SHAKE tape is read from byte k//8 ('lsb' or 'msb')
  * matrix_apply      : whether the generated A_j is applied as A.v or A^T.v (Alg. 4 says matrices are
                        produced "in their transposed form")
  * tape_matrix_order : whether the n^2 tape bits of one matrix fill T row-major ('row') or
                        column-major ('col'); the reference C stores each squeezed word as a COLUMN
                        of the matrix, which is the 'col' reading
  * lu_order          : whether A = L.U or A = U.L, with L = strictlower(T)+I, U = strictupper(T)+I
Field-element <-> vector conversion follows §4.1.1 exactly: bit i <-> coefficient of X^i <-> entry i.
Everything algebraic (S-box relations, kernel identities) is convention-independent and is checked
by `selftest()`.
"""
from __future__ import annotations
import hashlib
from dataclasses import dataclass
from math import gcd
from .gf2n import GF2n

AIM3_PARAMS = {              # Table 2
    "AIM3-I":   dict(n=128, exps=(3, 7, 11)),
    "AIM3-III": dict(n=192, exps=(11, 23, 47)),
    "AIM3-V":   dict(n=256, exps=(3, 7, 11, 19)),
}


@dataclass(frozen=True)
class Conventions:
    tape_bit_order: str = "lsb"      # 'lsb' | 'msb'
    matrix_apply: str = "A"          # 'A' (x = A.pt) | 'AT' (x = A^T.pt)
    tape_matrix_order: str = "col"   # 'row' | 'col'
    lu_order: str = "LU"             # 'LU' (A = L.U) | 'UL' (A = U.L)


#: all encoding conventions that the spec text leaves open; E00 pins one of them against the KAT.
ALL_CONVENTIONS = [Conventions(b, m, o, l)
                   for b in ("lsb", "msb") for m in ("A", "AT")
                   for o in ("col", "row") for l in ("LU", "UL")]

_BITREV = [int(f"{b:08b}"[::-1], 2) for b in range(256)]


def _shake(lam: int, data: bytes, nbits: int) -> bytes:
    nbytes = (nbits + 7) // 8
    if lam == 128:
        return hashlib.shake_128(data).digest(nbytes)
    return hashlib.shake_256(data).digest(nbytes)


class AIM3:
    def __init__(self, n: int, exps, lam: int | None = None, modulus: int | None = None,
                 conv: Conventions = Conventions()):
        self.n, self.exps, self.ell = n, tuple(exps), len(exps) - 1
        assert all(gcd(e, n) == 1 for e in self.exps), "spec requires gcd(e,n)=1"
        self.lam = lam or (128 if n <= 128 else 256)   # SHAKE128 for L1, SHAKE256 for L3/L5
        self.F = GF2n(n, modulus)
        self.conv = conv

    # ---------------- Algorithm 4: GenerateLinear(iv) ----------------
    def tape_bytes(self, iv: bytes) -> bytes:
        n, ell = self.n, self.ell
        nbits = 2 * ell * n * n + (ell + 1) * n
        tape = _shake(self.lam, iv, nbits)
        if self.conv.tape_bit_order == "msb":
            tape = bytes(_BITREV[b] for b in tape)
        return tape

    def tape_bits(self, iv: bytes):
        """the tape as a list of bits (bit k = k-th bit of the stream under `tape_bit_order`)."""
        n, ell = self.n, self.ell
        nbits = 2 * ell * n * n + (ell + 1) * n
        tape = self.tape_bytes(iv)
        return [(tape[k >> 3] >> (k & 7)) & 1 for k in range(nbits)]

    def generate_linear(self, iv: bytes):
        """(A_0..A_{2l-1}) as lists of row bitmasks, and (b_0..b_l) as field ints.

        A = L.U (or U.L) with L = strictlower(T)+I, U = strictupper(T)+I, where the n x n bit matrix
        T is read off the tape row-major or column-major; see `Conventions`.  Requires 8 | n so that
        each matrix row starts on a byte boundary (true for every AIM3 parameter set and every toy n
        used with the KAT harness; toy n not divisible by 8 falls back to the bitwise path).
        """
        n, ell = self.n, self.ell
        conv = self.conv
        mats = []
        if n % 8 == 0:
            tape = self.tape_bytes(iv)
            rowbytes = n // 8
            for j in range(2 * ell):
                off = j * n * rowbytes
                T = [int.from_bytes(tape[off + r * rowbytes: off + (r + 1) * rowbytes], "little")
                     for r in range(n)]
                mats.append(self._assemble(T))
            cb = 2 * ell * n * rowbytes
            consts = [int.from_bytes(tape[cb + j * rowbytes: cb + (j + 1) * rowbytes], "little")
                      for j in range(ell + 1)]
            return mats, consts
        bits = self.tape_bits(iv)
        for j in range(2 * ell):
            off = j * n * n
            T = [sum(bits[off + r * n + c] << c for c in range(n)) for r in range(n)]
            mats.append(self._assemble(T))
        cb = 2 * ell * n * n
        consts = [sum(bits[cb + j * n + r] << r for r in range(n)) for j in range(ell + 1)]
        return mats, consts

    @staticmethod
    def _transpose(rows, n):
        out = [0] * n
        for i, v in enumerate(rows):
            while v:
                c = (v & -v).bit_length() - 1
                out[c] |= 1 << i
                v &= v - 1
        return out

    def _assemble(self, T):
        """T (list of n row bitmasks) -> the matrix A as a list of n row bitmasks."""
        n, conv = self.n, self.conv
        if conv.tape_matrix_order == "col":
            T = self._transpose(T, n)
        low = [(T[i] & ((1 << i) - 1)) | (1 << i) for i in range(n)]          # strict lower + I
        up = [((T[i] >> (i + 1)) << (i + 1)) | (1 << i) for i in range(n)]    # strict upper + I
        first, second = (low, up) if conv.lu_order == "LU" else (up, low)
        A = []
        for i in range(n):
            acc, v = 0, first[i]
            while v:
                k = (v & -v).bit_length() - 1
                acc ^= second[k]
                v &= v - 1
            A.append(acc)
        return self._transpose(A, n) if conv.matrix_apply == "AT" else A

    @staticmethod
    def apply(Arows, v: int) -> int:
        out = 0
        for i, row in enumerate(Arows):
            out |= (bin(row & v).count("1") & 1) << i
        return out

    # ---------------- Algorithm 5: NoZeroInpSB(pt, iv) ----------------
    def evaluate(self, pt: int, iv: bytes, linear=None, want_internals=False):
        """returns ct, or None (= bottom) if some S-box input is zero. With want_internals, also returns
        (x_0..x_l, y_0..y_l)."""
        mats, consts = linear or self.generate_linear(iv)
        F, ell = self.F, self.ell
        xs, ys = [], []
        for j in range(ell):
            xj = self.apply(mats[j], pt) ^ consts[j]
            xs.append(xj)
            ys.append(F.sbox(xj, self.exps[j]))
        xl = consts[ell]
        for j in range(ell):
            xl ^= self.apply(mats[ell + j], ys[j])
        xs.append(xl)
        ys.append(F.sbox(xl, self.exps[ell]))
        ct = pt ^ ys[ell]
        if any(x == 0 for x in xs):
            ct = None
        return (ct, xs, ys) if want_internals else ct

    # ---------------- key generation (AIMer_keygen loop) ----------------
    def keygen(self, rng):
        n = self.n
        while True:
            pt = rng.getrandbits(n)
            iv = rng.getrandbits(n).to_bytes((n + 7) // 8, "big")
            ct = self.evaluate(pt, iv)
            if ct is not None:
                return pt, iv, ct

    # ---------------- convention-independent self test ----------------
    def selftest(self, rng, trials: int = 3) -> bool:
        F = self.F
        for _ in range(trials):
            pt, iv, ct = self.keygen(rng)
            ct2, xs, ys = self.evaluate(pt, iv, want_internals=True)
            assert ct2 == ct
            for x, y, e in zip(xs, ys, self.exps):
                r = F.r_inv(e)
                # P1: xy = x^{2^e} + 1
                assert F.mul(x, y) == F.frob(x, e) ^ 1
                # kernel variable a = y^r:  B_e(a,x) = a,  y a = a^{2^e}
                if y:
                    a = F.pow(y, r)
                    B = F.mul(a, F.frob(x, e)) ^ F.mul(F.frob(a, e), x)
                    assert B == a and F.mul(y, a) == F.frob(a, e)
                    assert F.sbox(x ^ a, e) == y            # second preimage
                    assert F.trace(F.pow(y, -r)) == 0       # image criterion
            assert ys[-1] == pt ^ ct
        return True


if __name__ == "__main__":
    import random, sys, time
    rng = random.Random(1)
    for name, p in AIM3_PARAMS.items():
        t = time.time()
        A = AIM3(p["n"], p["exps"])
        pt, iv, ct = A.keygen(rng)
        A.selftest(rng, trials=1)
        print(f"{name}: keygen+selftest OK  ({time.time()-t:.1f}s)  ct={ct:0{p['n']//4}x}"[:120])
    print("NOTE: encoding conventions still to be pinned against the official KAT (experiment E0).")
