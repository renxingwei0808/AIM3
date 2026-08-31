#!/usr/bin/env python3
"""E11/E00 — the bit cost of ONE forward AIM3 evaluation (COMPLEXITY_ACCOUNTING section 1).

Replaces the 2^14 placeholder.  Two independent inputs, and they agree:

 * **Spec Table 5** (AIMer v3.0 English specification, p.42): the operation counts of the forward
   evaluation, taken from the addition chains of the AIM3 S-boxes (the spec computes x^{-1} by
   Itoh-Tsujii and evaluates the Mersenne part x^{2^e-1} in parallel with it).

       scheme     F_{2^n} mul (count, depth)   F_{2^n} square (count, depth)   Grover total
       AIM3-I         47, 24                        399, 254                     2^165.6
       AIM3-III       59, 26                        651, 382                     2^232.5
       AIM3-V         73, 28                       1056, 510                     2^298.0

 * **An independent derivation** from the same chains (`itoh_tsujii_mults` below): Itoh-Tsujii needs
   floor(log2(n-1)) + wt(n-1) - 1 multiplications for x^{2^{n-1}-1}, the Mersenne part needs
   floor(log2 e) + wt(e) - 1, and each S-box ends with one multiplication.  That gives 50 / 63 / 77
   multiplications against the spec's 47 / 59 / 73 -- the 3-4 difference is the chain sharing between
   the Itoh-Tsujii tower and the Mersenne tower (both build x^{2^k-1}).  The spec numbers are used.

Cost model (both conventions are reported; constants are stated, not hidden):

    squaring          4n bit ops        (bit-interleaving is free; only the modular reduction costs,
                                         and every spec modulus is a pentanomial/trinomial)
    reduction          4n bit ops
    multiplication    schoolbook: n^2 AND + (n-1)^2 XOR + reduction
                      Karatsuba : 6 * n^{log2 3} + reduction   (the usual AND+XOR estimate)
    n x n binary matrix-vector product   n^2 AND + n(n-1) XOR
    field addition    n XOR

IMPORTANT modelling point for B3: in a single-key exhaustive search the iv is FIXED, so the 2*l
matrices A_j are generated once from SHAKE and reused; the SHAKE tape (2*l*n^2 + (l+1)*n bits) is
therefore NOT part of the per-candidate cost.  What is per-candidate is the S-box chain plus the
2*l matrix-vector products.  Both are reported.

Usage: python3 aim3_eval_cost.py [--csv]     ->  results/E11_theory_eval_cost.csv
"""
import sys, os, time
from math import log2, floor
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record_table

#: spec Table 5, p.42 -- (multiplications, mul depth, squarings, square depth, Grover log2 cost)
SPEC_TABLE5 = {
    "AIM3-I":   dict(n=128, ell=2, exps=(3, 7, 11),      mul=47, mul_depth=24, sqr=399,  sqr_depth=254, grover=165.6),
    "AIM3-III": dict(n=192, ell=2, exps=(11, 23, 47),    mul=59, mul_depth=26, sqr=651,  sqr_depth=382, grover=232.5),
    "AIM3-V":   dict(n=256, ell=3, exps=(3, 7, 11, 19),  mul=73, mul_depth=28, sqr=1056, sqr_depth=510, grover=298.0),
}


def itoh_tsujii_mults(n):
    """multiplications for x^{2^{n-1}-1} by Itoh-Tsujii (then one squaring gives x^{-1})."""
    k = n - 1
    return floor(log2(k)) + bin(k).count("1") - 1


def mersenne_mults(e):
    """multiplications for x^{2^e-1} by the same addition-chain rule."""
    return floor(log2(e)) + bin(e).count("1") - 1


def derived_counts(n, exps):
    """independent check of Table 5's multiplication count (no chain sharing assumed)."""
    return len(exps) * itoh_tsujii_mults(n) + sum(mersenne_mults(e) for e in exps) + len(exps)


def mul_bits(n, mode):
    red = 4 * n
    if mode == "schoolbook":
        return n * n + (n - 1) ** 2 + red
    if mode == "karatsuba":
        return round(6 * n ** log2(3)) + red
    raise ValueError(mode)


def eval_cost(name, mode):
    p = SPEC_TABLE5[name]
    n, ell = p["n"], p["ell"]
    sbox = p["mul"] * mul_bits(n, mode) + p["sqr"] * (4 * n)
    matvec = (2 * ell) * (n * n + n * (n - 1))
    adds = ((ell + 1) + ell + 1) * n            # the b_i, the last-layer sum, the feed-forward
    total = sbox + matvec + adds
    return dict(scheme=name, n=n, ell=ell, exps=";".join(map(str, p["exps"])), mode=mode,
                spec_mul=p["mul"], spec_sqr=p["sqr"], derived_mul_no_sharing=derived_counts(n, p["exps"]),
                bits_mul=mul_bits(n, mode), bits_sbox=sbox, bits_matvec=matvec, bits_add=adds,
                bits_total=total, log2_eval=round(log2(total), 2),
                log2_B3=round(n + log2(total), 2), grover_log2=p["grover"])


if __name__ == "__main__":
    t0 = time.time()
    print("# params: spec Table 5 (AIMer v3.0 English spec, p.42); modes = schoolbook, karatsuba")
    rows = []
    print(f"{'scheme':>9} {'n':>4} {'mul':>4} {'sqr':>5} {'derived':>8} {'mode':>11} "
          f"{'bits/mul':>9} {'bits/eval':>11} {'log2 eval':>10} {'log2 B3':>9}")
    for name in SPEC_TABLE5:
        for mode in ("schoolbook", "karatsuba"):
            r = eval_cost(name, mode)
            rows.append(r)
            print(f"{r['scheme']:>9} {r['n']:>4} {r['spec_mul']:>4} {r['spec_sqr']:>5} "
                  f"{r['derived_mul_no_sharing']:>8} {mode:>11} {r['bits_mul']:>9} "
                  f"{r['bits_total']:>11} {r['log2_eval']:>10.2f} {r['log2_B3']:>9.2f}")
    print()
    print("B3 = 2^n forward evaluations. The old placeholder was 2^14 bit ops per evaluation;")
    print("the measured chains give 2^19.9 .. 2^23.5, i.e. B3 is 6 to 9 bits higher than assumed.")
    if "--csv" in sys.argv:
        record_table("E11_theory_eval_cost.csv", rows, seconds=time.time() - t0)
        print("written results/E11_theory_eval_cost.csv")
