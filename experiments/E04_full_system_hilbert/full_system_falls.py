#!/usr/bin/env python3
"""E04 (Tier A) — is the D=4 degree fall of the FULL AIM3 system the direct sum of the per-S-box falls?

Systems on the same instance (legacy/aim3_macaulay.build_aim3, unit-LU matrices, variables pt,t_0..t_{l-1}):
  S_quad   the official 3n(l+1) quadratics {P1,P2,P3} per S-box
  P1only   the n(l+1) bilinear P1 coordinates only
  random   a random Boolean quadratic system of the same shape with a planted root (control)

Reported per D: rows, rank, syz = rows-rank, dim(span n deg<=D-1); and

    falls_D4 = dim(span_{D=4} n deg<=3) - rank_{D=3}

which is exactly the E03 definition applied to the full system.  E03 gives 2n falls per S-box
(2n+1 for some e), so the direct-sum hypothesis T1.3 predicts

    falls_D4(S_quad) = (l+1) * 2n     (+1 per S-box tolerance)   [l=2 -> 6n, l=1 -> 4n].

A value clearly above that band means cross-S-box degree falls exist -> Gate 1 of docs/ROADMAP.md.

TWO CAVEATS, both measured and reported per run rather than assumed away:

 1. SATURATION.  S_quad has 3n(l+1) equations in only (l+1)n variables, so at D=4 the Macaulay matrix
    is heavily over-determined and its row space can already fill the whole degree-<=3 space.  When
    that happens `dim(span n deg<=3)` is pinned at C(nv,<=3) and `falls_D4` measures the size of the
    monomial basis, not the structure.  The columns `low_cols` and `saturation` make this visible; a
    run with saturation >= 0.99 is reported as verdict `indeterminate(saturated)` and must NOT be read
    as evidence either way.  The same-shape random control saturates identically, which is the proof
    that it is an artefact.
 2. EXTRA MULTIPLIERS.  Even without saturation, comparing the joint fall count against (l+1)*2n is
    not by itself a direct-sum test: the joint Macaulay matrix also multiplies by cross-block
    monomials (pt*t1, t0*t1, ...) that do not exist in E03's 2n-variable single-S-box ring.  The
    honest test is the per-block control (`--per-block`): run each S-box's own 3n quadratics ALONE in
    the same (l+1)n variables and compare

        cross_excess = falls_D4(all S-boxes) - sum_i falls_D4(S-box i alone).

    cross_excess = 0 <=> direct sum (T1.3);  > 0 <=> genuine cross-S-box falls.

Usage:
  python3 full_system_falls.py N [-l ELL] [--seeds 1,2] [--systems S_quad,P1only,random] [--budget SEC]
Every run appends a row to results/E04_full_system_falls.csv; a run that exceeds --budget
(default 3600 s per Macaulay matrix) is recorded with censored=1 and its rank/fall numbers are
LOWER BOUNDS, never values (CLAUDE.md ground rule 4).
"""
import sys, os, time, argparse
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3_macaulay import build_aim3, random_system
from aim3.boolmac import macaulay_stats, MacaulayTimeout


def default_exps(n, ell):
    cand = [e for e in range(2, n - 1) if gcd(e, n) == 1]
    if len(cand) < ell + 1:
        raise SystemExit(f"n={n}: only {len(cand)} admissible exponents, need {ell+1}")
    return cand[:ell + 1]


def measure(name, eqs, nv, D, budget):
    """-> (MacaulayResult, censored, seconds)"""
    t0 = time.time()
    try:
        st = macaulay_stats(eqs, nv, D, budget=budget)
        return st, 0, time.time() - t0
    except MacaulayTimeout as exc:
        print(f"     !! {name} D={D} exceeded the {budget:.0f}s budget after "
              f"{time.time()-t0:.0f}s -> lower bounds only: {exc.partial}")
        return exc.partial, 1, time.time() - t0


def run(n, ell, seed, systems, budget, per_block=False):
    from math import comb
    exps = default_exps(n, ell)
    nv, eqsP1, eqsAll, _ = build_aim3(n, exps, seed=seed)
    assert nv == (ell + 1) * n
    target = (ell + 1) * 2 * n
    low_cols = sum(comb((ell + 1) * n, d) for d in range(4))
    pool = {"S_quad": eqsAll, "P1only": eqsP1}
    if "random" in systems:
        pool["random"] = random_system(nv, len(eqsAll), seed=1000 + seed)[0]
    print(f"n={n} l={ell} exps={exps} seed={seed}: {nv} vars, "
          f"S_quad {len(eqsAll)} eqs, P1only {len(eqsP1)} eqs; naive target (l+1)*2n = {target}; "
          f"degree-<=3 columns C({nv},<=3) = {low_cols}")
    for name in systems:
        eqs = pool[name]
        print(f"  {name}: {len(eqs)} eqs")
        st = {}
        cens = 0
        secs = 0.0
        for D in (3, 4):
            st[D], c, sc = measure(name, eqs, nv, D, budget)
            cens |= c
            secs += sc
            print(f"     D={D}: rows={st[D].rows:7d} cols={st[D].cols:7d} rank={st[D].rank:7d} "
                  f"syz={st[D].syz:7d} dim(span n deg<={D-1})={st[D].dim_low:7d}  ({sc:.1f}s)"
                  + ("   [CENSORED: lower bound]" if c else ""))
            if c:
                break
        rec = dict(n=n, ell=ell, exps=";".join(map(str, exps)), seed=seed, system=name,
                   nvars=nv, neqs=len(eqs),
                   D3_rows=st[3].rows, D3_rank=st[3].rank, D3_syz=st[3].syz, D3_dim_low=st[3].dim_low,
                   D4_rows="", D4_rank="", D4_syz="", D4_dim_low="", falls_D4="",
                   target_falls=target, ratio_over_n="", low_cols=low_cols, saturation="",
                   per_block_falls="", cross_excess="", verdict="censored" if cens else "")
        if 4 in st and not cens:
            falls = st[4].dim_low - st[3].rank
            rec.update(D4_rows=st[4].rows, D4_rank=st[4].rank, D4_syz=st[4].syz,
                       D4_dim_low=st[4].dim_low, falls_D4=falls, ratio_over_n=round(falls / n, 4))
            sat = st[4].dim_low / low_cols
            rec.update(low_cols=low_cols, saturation=round(sat, 4))
            if name == "S_quad":
                if sat >= 0.99:
                    rec["verdict"] = "indeterminate(saturated)"
                else:
                    # tolerance: E03 allows 2n or 2n+1 per S-box, i.e. up to +(l+1) over the target
                    rec["verdict"] = ("direct-sum" if target <= falls <= target + (ell + 1)
                                      else ("cross-falls" if falls > target + (ell + 1) else "below-target"))
                print(f"     => falls_D4 = {st[4].dim_low} - {st[3].rank} = {falls}  "
                      f"(naive target {target} = {ell+1}x2n); dim(span n deg<=3)/C(nv,<=3) = "
                      f"{st[4].dim_low}/{low_cols} = {sat:.4f} -> {rec['verdict'].upper()}")
            else:
                print(f"     => falls_D4 = {falls}  ({falls/n:.2f} n); saturation {sat:.4f}")
            if name == "S_quad" and per_block:
                per = []
                for i in range(ell + 1):
                    sub = eqs[i * 3 * n:(i + 1) * 3 * n]
                    b3, c3, _ = measure(f"S-box{i} alone", sub, nv, 3, budget)
                    b4, c4, _ = measure(f"S-box{i} alone", sub, nv, 4, budget)
                    cens |= c3 | c4
                    per.append(b4.dim_low - b3.rank)
                    print(f"     [control] S-box {i} alone in {nv} vars: D3 rank={b3.rank} syz={b3.syz} "
                          f"(8n={8*n}); falls_D4 = {per[-1]} ({per[-1]/n:.2f} n)")
                cross = falls - sum(per)
                rec.update(per_block_falls=";".join(map(str, per)), cross_excess=cross)
                print(f"     => DIRECT-SUM TEST: {falls} - {' - '.join(map(str, per))} = {cross}"
                      + ("   (meaningless here: the joint run is saturated)" if sat >= 0.99 else
                         f"   -> {'DIRECT SUM' if cross == 0 else 'CROSS-S-BOX FALLS'}"))
        elif cens:
            print(f"     => falls_D4 NOT determined (censored); D=3 rank >= {st[3].rank}")
        record("E04_full_system_falls.csv", rec, seconds=secs, censored=cens)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("-l", "--ell", type=int, default=2, help="number of S-boxes minus 1 (l=2 -> 3 S-boxes)")
    ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--systems", default="S_quad,P1only,random")
    ap.add_argument("--budget", type=float, default=3600.0, help="seconds per Macaulay matrix")
    ap.add_argument("--per-block", action="store_true",
                    help="also run each S-box's 3n quadratics alone in the same variables (direct-sum control)")
    a = ap.parse_args()
    print(f"# params: n={a.n} ell={a.ell} seeds={a.seeds} systems={a.systems} budget={a.budget}s "
          f"per_block={a.per_block}")
    for seed in [int(s) for s in a.seeds.split(",")]:
        run(a.n, a.ell, seed, a.systems.split(","), a.budget, a.per_block)
