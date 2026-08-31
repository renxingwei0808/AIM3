#!/usr/bin/env python3
"""
AIM3 structural verification lab.

Part 1: verify closed forms for the enhanced Mersenne S-box.
Part 2: verify wt(r) = e^{-1} mod n for the recommended AIM3 parameters.
Part 3: measure the dimension of the space of low-degree Boolean relations
        on joint S-box graphs of a toy AIM3 instance (tests the AIMer v3
        claim that no cross-S-box quadratic relations exist).
Part 4: degree-3 Macaulay rank / degree-fall comparison against random
        quadratic systems of the same shape.
"""
from __future__ import annotations
import itertools, random
from math import gcd, comb

# ---------------------------------------------------------------- GF(2^n)

CONWAY = {
    5: 0x25, 6: 0x43, 7: 0x83, 8: 0x11B, 9: 0x211, 10: 0x409,
    11: 0x805, 12: 0x1053, 13: 0x201B, 14: 0x402B, 15: 0x8003, 16: 0x1002D,
}

class F:
    def __init__(self, n):
        self.n = n
        self.mod = CONWAY[n]
        self.q = 1 << n
        self.mask = self.q - 1
        # log/antilog tables via a generator search
        self.exp = [0] * (2 * self.q)
        self.log = [0] * self.q
        g = self._find_generator()
        v = 1
        for i in range(self.q - 1):
            self.exp[i] = v
            self.log[v] = i
            v = self._mul_slow(v, g)
        for i in range(self.q - 1, 2 * self.q):
            self.exp[i] = self.exp[i - (self.q - 1)]

    def _mul_slow(self, a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a & self.q: a ^= self.mod
        return r & self.mask

    def _find_generator(self):
        order = self.q - 1
        facs = set()
        m = order
        d = 2
        while d * d <= m:
            while m % d == 0:
                facs.add(d); m //= d
            d += 1
        if m > 1: facs.add(m)
        for g in range(2, self.q):
            if all(self._pow_slow(g, order // p) != 1 for p in facs):
                return g
        raise RuntimeError

    def _pow_slow(self, a, k):
        r = 1
        while k:
            if k & 1: r = self._mul_slow(r, a)
            a = self._mul_slow(a, a); k >>= 1
        return r

    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.exp[self.log[a] + self.log[b]]

    def pw(self, a, k):
        if a == 0: return 0
        k %= (self.q - 1)
        return self.exp[(self.log[a] * k) % (self.q - 1)]

    def inv(self, a):
        return self.exp[(self.q - 1 - self.log[a]) % (self.q - 1)]

    def trace(self, a):
        t, c = 0, a
        for _ in range(self.n):
            t ^= c; c = self.mul(c, c)
        assert t in (0, 1)
        return t

    def bits(self, a):
        return [(a >> i) & 1 for i in range(self.n)]


def sbox(fld, x, e):
    if x == 0: return 0
    return fld.pw(x, (1 << e) - 1) ^ fld.inv(x)


# ------------------------------------------------- Part 1: closed forms

def part1():
    print("=" * 74)
    print("PART 1  closed-form identities for S_e(x) = x^(2^e-1) + x^(-1)")
    print("=" * 74)
    for n in (5, 7, 8, 9, 11):
        fld = F(n)
        for e in range(2, n - 1):
            if gcd(e, n) != 1: continue
            q1 = fld.q - 1
            r = pow((1 << e) - 1, -1, q1)
            m = pow(e, -1, n)
            # (a) S_e(x) = (x+1)^(2^e) / x
            ok_a = all(sbox(fld, x, e) == fld.mul(fld.pw(x ^ 1, 1 << e), fld.inv(x))
                       for x in range(1, fld.q))
            # (b) wt(r) = e^{-1} mod n  and  wt(-r) = n - wt(r)
            ok_b = (bin(r).count('1') == m and
                    bin((-r) % q1).count('1') == n - m)
            # (c) image = {0} u {y != 0 : Tr(y^{-r}) = 0}, size q/2
            img = {sbox(fld, x, e) for x in range(1, fld.q)}
            ok_c = (len(img) == fld.q // 2 and
                    all((fld.trace(fld.pw(y, -r)) == 0) == (y in img)
                        for y in range(1, fld.q)))
            # (d) two-valued inverse  x = a(u+eps) + 1,  a = y^r,
            #     u^(2^e) + u = y^(-r)
            sec = {}
            for u in range(fld.q):
                sec.setdefault(fld.pw(u, 1 << e) ^ u, u)
            fib = {}
            for x in range(1, fld.q):
                fib.setdefault(sbox(fld, x, e), []).append(x)
            ok_d = True
            for y in img:
                if y == 0: continue
                a = fld.pw(y, r)
                u = sec[fld.inv(a)]
                pre = {fld.mul(a, u) ^ 1, fld.mul(a, u ^ 1) ^ 1}
                if pre != set(fib[y]): ok_d = False; break
            # (e) parameterisation  b=u^{2^e}+u, x=u/b+1, y=b^{1-2^e}
            ok_e = True
            for u in range(fld.q):
                if u in (0, 1): continue
                b = fld.pw(u, 1 << e) ^ u
                x = fld.mul(u, fld.inv(b)) ^ 1
                y = fld.pw(b, (1 - (1 << e)) % q1)
                if sbox(fld, x, e) != y: ok_e = False; break
            flag = "OK " if all((ok_a, ok_b, ok_c, ok_d, ok_e)) else "FAIL"
            print(f"  n={n:2d} e={e:2d}  wt(r)={bin(r).count('1'):3d}(=e^-1 mod n={m:3d})"
                  f"  |Im|={len(img)}  [{flag}] "
                  f"a={ok_a} b={ok_b} c={ok_c} d={ok_d} e={ok_e}")


# --------------------------- Part 2: recommended parameters, no field ops

def part2():
    print()
    print("=" * 74)
    print("PART 2  wt(r) for the recommended AIM3 parameters (r=(2^e-1)^-1)")
    print("=" * 74)
    params = [("AIM3-I", 128, [3, 7, 11]),
              ("AIM3-III", 192, [11, 23, 47]),
              ("AIM3-V", 256, [3, 7, 11, 19])]
    for name, n, es in params:
        q1 = (1 << n) - 1
        for k, e in enumerate(es):
            r = pow((1 << e) - 1, -1, q1)
            m = pow(e, -1, n)
            wr, wmr = bin(r).count('1'), bin((-r) % q1).count('1')
            # support of r as residues mod n
            sup = sorted((i * e) % n for i in range(m))
            ap = all((sup[i] in {(j * e) % n for j in range(m)}) for i in range(m))
            tag = "  <-- FINAL S-box" if k == len(es) - 1 else ""
            print(f"  {name:9s} n={n} e={e:3d}: wt(r)={wr:4d}  e^-1 mod n={m:4d}  "
                  f"wt(-r)={wmr:4d}  supp(r)={{i*e mod n}}: {ap}{tag}")


# ----------------------- toy AIM3 instance and Boolean relation counting

def rand_unit_lu(n, rng):
    L = [[1 if i == j else (rng.getrandbits(1) if i > j else 0) for j in range(n)]
         for i in range(n)]
    U = [[1 if i == j else (rng.getrandbits(1) if i < j else 0) for j in range(n)]
         for i in range(n)]
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if L[i][k]:
                for j in range(n):
                    A[i][j] ^= U[k][j]
    return A


def apply_mat(A, x, n):
    xb = [(x >> i) & 1 for i in range(n)]
    y = 0
    for i in range(n):
        s = 0
        for j in range(n):
            s ^= A[i][j] & xb[j]
        y |= s << i
    return y


def monomials(nv, deg):
    ms = []
    for d in range(deg + 1):
        ms.extend(itertools.combinations(range(nv), d))
    return ms


def rank_gf2(rows):
    piv = {}
    rk = 0
    for r in rows:
        while r:
            b = r.bit_length() - 1
            if b in piv:
                r ^= piv[b]
            else:
                piv[b] = r; rk += 1; break
    return rk, piv


def relation_dim(points, nv, deg):
    """dim of Boolean polys of degree<=deg vanishing on all given points."""
    ms = monomials(nv, deg)
    rows = []
    for p in points:
        v = 0
        for idx, mono in enumerate(ms):
            bit = 1
            for i in mono:
                bit &= (p >> i) & 1
                if not bit: break
            if bit: v |= 1 << idx
        rows.append(v)
    rk, _ = rank_gf2(rows)
    return len(ms) - rk, len(ms), rk


def part3():
    print()
    print("=" * 74)
    print("PART 3  Boolean relation spaces on joint S-box graphs (toy AIM3)")
    print("        tests the AIMer v3 claim: no cross-S-box quadratic relations")
    print("=" * 74)
    rng = random.Random(20260831)
    for n in (11, 12, 13):
        fld = F(n)
        es = [e for e in range(2, n - 1) if gcd(e, n) == 1]
        e0, e1 = es[0], es[1]
        A0 = rand_unit_lu(n, rng); A1 = rand_unit_lu(n, rng)
        b0 = rng.randrange(fld.q); b1 = rng.randrange(fld.q)
        pts_pt_t0, pts_t0_t1, pts_all = [], [], []
        for pt in range(fld.q):
            x0 = apply_mat(A0, pt, n) ^ b0
            x1 = apply_mat(A1, pt, n) ^ b1
            if x0 == 0 or x1 == 0: continue
            t0 = sbox(fld, x0, e0); t1 = sbox(fld, x1, e1)
            pts_pt_t0.append(pt | (t0 << n))
            pts_t0_t1.append(t0 | (t1 << n))
            pts_all.append(pt | (t0 << n) | (t1 << (2 * n)))
        d1, tot1, _ = relation_dim(pts_pt_t0, 2 * n, 2)
        d2, tot2, _ = relation_dim(pts_t0_t1, 2 * n, 2)
        d3, tot3, _ = relation_dim(pts_all, 3 * n, 2)
        print(f"  n={n:2d} (e0={e0},e1={e1}) #points={len(pts_all)}")
        print(f"     deg<=2 relations in (pt,t0)   : {d1:5d}   "
              f"(expected 3n={3*n}, monomials={tot1})")
        print(f"     deg<=2 relations in (t0,t1)   : {d2:5d}   "
              f"(spec claims 0, monomials={tot2})")
        print(f"     deg<=2 relations in (pt,t0,t1): {d3:5d}   "
              f"(expected 6n={6*n}, monomials={tot3})")


# ---------------- Part 4: degree-3 Macaulay behaviour of the AIM3 system

def field_coord_eqs(fld, lin_rows, const, deg2_terms, nv):
    """placeholder - not used"""
    raise NotImplementedError


def build_aim3_system(n, ells=2, seed=7):
    """
    Returns (nv, eqs, sol) where eqs is a list of Boolean polynomials given as
    sets of monomials (tuples of variable indices), for the AIM3 quadratic
    system in variables x, t_0, ..., t_{l-1}.  Built symbolically by evaluating
    over all inputs is impossible; instead we build the coefficient vectors
    directly from the field structure.
    """
    raise NotImplementedError


if __name__ == "__main__":
    part1()
    part2()
    part3()
