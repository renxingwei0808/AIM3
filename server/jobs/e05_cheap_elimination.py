#!/usr/bin/env python3
"""Batch-2 task 5: eliminate t₁ from the reduced system without computing a resultant.

Server batch 1 measured that `Resultant(Eq1, Eq2, t₁)` finishes at n = 7 in 294.5 s and does **not**
finish at n = 8 inside two hours, leaving one completed size and no slope.  Its §8 proposed two
cheaper routes; this script implements the one that works.

**The structure that makes it cheap.**  Eq1 is *separable*:

    Eq1 = W(t₀) + V(t₁),    W(t₀) = (A₁A₀⁻¹)(G_{e₀,ε₀}(t₀) + b₀) + b₁,   V(t₁) = G_{e₁,ε₁}(t₁)

so its variety is the set of (t₀,t₁) with W(t₀) = V(t₁) — computable by **tabulating V over F_q\\{0\\}
and looking W(t₀) up**, without ever forming an eliminated polynomial.  Each surviving pair is then
tested against Eq2 by field arithmetic.  Total Θ(2ⁿ) field operations and Θ(2ⁿ) memory.

`reduced_system_toy.count_roots` already does exactly this; what was missing was to *time* it as an
elimination method and to report the curve.  This script builds the same instance batch 1 exported
to Magma, times the two phases separately (build, eliminate), and records the root count so that
the answers can be compared with the resultant's where both exist.

**What this does and does not settle.**  It gives the cost of *solving* the reduced system, which
is what track T3 needs; it does not produce an eliminated univariate polynomial, so there is no
`deg` column and no `log2(deg)` slope — that quantity only exists for the resultant route, which
`results/E05_server_resultants.csv` shows is out of reach from n = 8.  Stated plainly rather than
substituted for.

Usage:
  python3 server/jobs/e05_cheap_elimination.py --ns 7,8,...,20 [--seed 1]
                                               [--csv E05_server_cheap_elimination.csv]
"""
import argparse
import math
import os
import random
import sys
import time
from math import gcd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_EXP = os.path.join(ROOT, "experiments")
for p in [os.path.join(ROOT, "src"), os.path.join(ROOT, "legacy"), _EXP,
          os.path.join(_EXP, "E05_two_valued_inverse")]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _common import record                                          # noqa: E402
from reduced_system_toy import build_reduced_system, verify, count_roots   # noqa: E402


def least_squares(xs, ys):
    if len(xs) < 2:
        return None
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default=",".join(str(n) for n in range(7, 21)))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--budget", type=float, default=3600.0,
                    help="per-n wall-clock cap; a size over it is recorded censored and stops the "
                         "escalation, as in batch 1")
    ap.add_argument("--csv", default="E05_server_cheap_elimination.csv")
    a = ap.parse_args()
    ns = [int(v) for v in a.ns.split(",")]
    print(f"# params: ns={ns} seed={a.seed} budget={a.budget}s")
    print("# method: tabulate V(t1) over F_q*, look up W(t0), test each hit against Eq2 -- "
          "Theta(2^n) field ops, no eliminated polynomial is formed")

    pts = []
    for n in ns:
        exps = [e for e in range(2, n) if gcd(e, n) == 1][:3]
        if len(exps) < 3:
            print(f"  n={n}: fewer than 3 admissible exponents, skipped")
            continue
        print(f"\n== n={n} exps={exps}", flush=True)
        t0 = time.time()
        rng = random.Random(1000 + a.seed)
        eq1, eq2, ctx = build_reduced_system(n, exps, a.seed, rng)
        verify(eq1, eq2, ctx)
        t_build = time.time() - t0
        s1, s2 = eq1.stats(), eq2.stats()
        print(f"   build + verify: {t_build:.1f}s; Eq1 {s1['terms']} terms deg_t0={s1['deg_t0']} "
              f"deg_t1={s1['deg_t1']}, Eq2 {s2['terms']} terms deg_t1={s2['deg_t1']}", flush=True)
        if t_build > a.budget:
            print(f"   !! build alone exceeded the {a.budget:.0f}s budget -- censored, stopping")
            record(a.csv, dict(n=n, exps=";".join(map(str, exps)), seed=a.seed,
                               build_seconds=round(t_build, 3), eliminate_seconds="",
                               roots="", spurious="", key_is_root="",
                               eq1_terms=s1["terms"], eq2_terms=s2["terms"],
                               eq1_deg_t1=s1["deg_t1"], eq2_deg_t1=s2["deg_t1"],
                               method="separable table lookup"),
                   seconds=t_build, censored=1, tier="B")
            break

        t1 = time.time()
        roots = count_roots(eq1, eq2, ctx)
        t_elim = time.time() - t1
        key = tuple(ctx["ts"])
        ok = key in roots
        print(f"   eliminate: {t_elim:.2f}s -> {len(roots)} roots in (F*)^2 "
              f"(true key present: {ok}; spurious {len(roots)-1})", flush=True)
        record(a.csv, dict(n=n, exps=";".join(map(str, exps)), seed=a.seed,
                           build_seconds=round(t_build, 3), eliminate_seconds=round(t_elim, 3),
                           roots=len(roots), spurious=len(roots) - 1, key_is_root=int(ok),
                           eq1_terms=s1["terms"], eq2_terms=s2["terms"],
                           eq1_deg_t1=s1["deg_t1"], eq2_deg_t1=s2["deg_t1"],
                           method="separable table lookup"),
               seconds=t_build + t_elim, censored=0, tier="B")
        pts.append((n, t_elim))

    print("\n# fit (COMPLEXITY_ACCOUNTING section 4: >= 3 completed sizes before a slope is quoted)")
    if len(pts) < 3:
        print(f"  only {len(pts)} completed sizes -- no slope")
        return
    xs = [n for n, _ in pts]
    sl = least_squares(xs, [math.log2(max(t, 1e-9)) for _, t in pts])
    print(f"  completed n = {xs}")
    print(f"  slope of log2(elimination seconds) vs n = {sl:.3f}")
    print(f"  (the method is Theta(2^n) by construction, so a slope near 1 is the expected answer; "
          f"it is measured rather than assumed)")
    print("  no log2(deg) slope: this route forms no eliminated polynomial. The resultant route, "
          "which does, is censored from n = 8 (results/E05_server_resultants.csv).")


if __name__ == "__main__":
    main()
