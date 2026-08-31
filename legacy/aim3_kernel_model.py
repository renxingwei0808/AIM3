#!/usr/bin/env python3
"""
Kernel-variable ("a-model") of the enhanced Mersenne S-box.

For x != 0 and y = S_e(x), let a = y^r with r = (2^e-1)^{-1} mod 2^n-1, i.e.
a is the unique nonzero element of ker(v -> v^{2^e} + y v).

Claims verified here:
  (K1)  B_e(a,x) := a x^{2^e} + a^{2^e} x  =  a          (bilinear in (a,x))
  (K2)  y a = a^{2^e}                                    (bilinear in (a,y))
  (K3)  B_e is symmetric, alternating, and B_e(u,v)=0 iff v in {0,u}
  (K4)  the two preimages of y are x and x + a
  (K5)  dim of Boolean relations of degree <= 2 vanishing on {(x,a)} and {(a,y)}
"""
from __future__ import annotations
import random
from math import gcd
from aim3_lab import F, sbox, relation_dim


def check(n):
    fld = F(n)
    q1 = fld.q - 1
    out = []
    for e in range(2, n - 1):
        if gcd(e, n) != 1:
            continue
        r = pow((1 << e) - 1, -1, q1)
        k1 = k2 = k4 = True
        pts_xa, pts_ay, pts_xat = [], [], []
        fib = {}
        for x in range(1, fld.q):
            y = sbox(fld, x, e)
            fib.setdefault(y, []).append(x)
        for x in range(1, fld.q):
            y = sbox(fld, x, e)
            if y == 0:            # x = 1, the exceptional point
                continue
            a = fld.pw(y, r)
            B = fld.mul(a, fld.pw(x, 1 << e)) ^ fld.mul(fld.pw(a, 1 << e), x)
            if B != a: k1 = False
            if fld.mul(y, a) != fld.pw(a, 1 << e): k2 = False
            if set(fib[y]) != {x, x ^ a}: k4 = False
            pts_xa.append(x | (a << n))
            pts_ay.append(a | (y << n))
            pts_xat.append(x | (a << n) | (y << (2 * n)))
        # K3
        k3 = True
        for _ in range(400):
            u = random.randrange(1, fld.q); v = random.randrange(fld.q)
            B1 = fld.mul(u, fld.pw(v, 1 << e)) ^ fld.mul(fld.pw(u, 1 << e), v)
            B2 = fld.mul(v, fld.pw(u, 1 << e)) ^ fld.mul(fld.pw(v, 1 << e), u)
            if B1 != B2: k3 = False
            if (B1 == 0) != (v in (0, u)): k3 = False
        d_xa = relation_dim(pts_xa, 2 * n, 2)[0]
        d_ay = relation_dim(pts_ay, 2 * n, 2)[0]
        d_xat = relation_dim(pts_xat, 3 * n, 2)[0] if 3 * n <= 33 else None
        out.append((e, k1, k2, k3, k4, d_xa, d_ay, d_xat))
    return out


if __name__ == "__main__":
    random.seed(0)
    print("kernel-variable model  a = S_e(x)^r   (r = (2^e-1)^{-1})")
    print(f"{'n':>3} {'e':>3} | K1 K2 K3 K4 | {'dimQE(x,a)':>11} {'dimQE(a,y)':>11} "
          f"{'dimQE(x,a,y)':>13}")
    for n in (7, 8, 9, 10, 11):
        for (e, k1, k2, k3, k4, d1, d2, d3) in check(n):
            print(f"{n:3d} {e:3d} |  {int(k1)}  {int(k2)}  {int(k3)}  {int(k4)} | "
                  f"{d1:11d} {d2:11d} {str(d3):>13}   (n={n}, 3n={3*n}, 4n={4*n})")
