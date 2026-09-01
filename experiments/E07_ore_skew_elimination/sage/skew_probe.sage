# Sage 10.x cross-check of the Tier-A measurement in ../skew_elimination_probe.py.
# UNTESTED ON THIS MACHINE (no Sage at Tier A) -- validate it against the committed Python numbers in
# results/E07_theory_skew_degrees.csv before trusting anything it prints for larger n (CLAUDE.md, tier B/C rules).
#
#   sage skew_probe.sage ../instance_n11.json
#
# It must reproduce, for the exported instance:
#   Lambda_i(pt) == c_i ;  deg right_gcd(Lambda_0, Lambda_1) == 0 ;
#   deg left_lcm == deg Lambda_0 + deg Lambda_1 ;  and exactly two solvable t_0 (one of them t[0]).
import json, sys

path = sys.argv[1] if len(sys.argv) > 1 else "../instance_n11.json"
D = json.load(open(path))
n, mod, e0, e1 = D["n"], D["modulus"], D["e"][0], D["e"][1]

P.<X> = GF(2)[]
f = P([(mod >> i) & 1 for i in range(n + 1)])          # bit i = coefficient of X^i
K.<w> = GF(2 ** n, modulus=f)
R.<phi> = OrePolynomialRing(K, K.frobenius_endomorphism())

def elt(v):                                            # int -> K, same bit convention
    return K([(v >> i) & 1 for i in range(n)])

def num(a):                                            # K -> int
    return sum(int(c) << i for i, c in enumerate(a.polynomial().list()))

A = [R([elt(c) for c in D["A_skew_coeffs"][i]]) for i in (0, 1)]
b = [elt(v) for v in D["b"]]
t = [elt(v) for v in D["t"]]
c = [elt(v) for v in D["c"]]
pt = elt(D["pt"])

Lam = [(phi ** ei + t[i]) * A[i] for i, ei in enumerate((e0, e1))]
print("Lambda_i(pt) == c_i :", [Lam[i](pt) == c[i] for i in (0, 1)])
print("deg Lambda        :", [Lam[i].degree() for i in (0, 1)])
g = Lam[0].right_gcd(Lam[1]); l = Lam[0].left_lcm(Lam[1])
print("deg gcrd, deg lclm:", g.degree(), l.degree(),
      "(sum of degrees:", Lam[0].degree() + Lam[1].degree(), ")")

# the solvable set C(., t_1): brute force over K, the same quantity the Python script measures
def solvable(t0):
    L0 = (phi ** e0 + t0) * A[0]
    c0 = K(1) + b[0] ** (2 ** e0) + t0 * b[0]
    M = matrix(GF(2), 2 * n, n)
    for k in range(n):
        v = w ** 0 * K([1 if j == k else 0 for j in range(n)])
        for r, val in ((0, L0(v)), (n, Lam[1](v))):
            for i, ci in enumerate(val.polynomial().list() + [0] * n):
                if i < n:
                    M[r + i, k] = ci
    rhs = vector(GF(2), [ (c0.polynomial().list() + [0]*n)[i] for i in range(n) ]
                      + [ (c[1].polynomial().list() + [0]*n)[i] for i in range(n) ])
    try:
        M.solve_right(rhs)
        return True
    except ValueError:
        return False

C = [x for x in K if solvable(x)]
print("|C| =", len(C), " true t_0 in C:", t[0] in C, " as ints:", sorted(num(x) for x in C))
