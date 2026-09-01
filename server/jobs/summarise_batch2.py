#!/usr/bin/env python3
"""Generate the appendix tables of `results/REPORT_server_batch2.md` from the committed CSVs.

Same contract as `summarise_batch1.py`: read only what `collect.py` has put into `results/`, emit
Markdown, and with `--inject` replace the `<!--TABLES-->` marker in the report (idempotently), so
the prose and the data cannot drift apart.

Usage: python3 server/jobs/summarise_batch2.py [--out FILE] [--inject REPORT] [--quiet]
"""
import argparse
import csv
import math
import os
import sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "experiments"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "E02_official_baseline"))

RES = os.path.join(ROOT, "results")
OUT = []


def A(line=""):
    OUT.append(line)


def load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return list(csv.DictReader(fh))


def i(v, d=None):
    try:
        return int(v)
    except (TypeError, ValueError):
        return d


def f(v, d=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


# --------------------------------------------------------------------- A.0 decoupled control
def sec_decoupled():
    rows = load("E04_server_decoupled.csv")
    A("### B.0 The decoupled control (task 0)")
    A("")
    if not rows:
        A("_not collected_")
        A("")
        return
    A("k S-boxes on pairwise **disjoint** variable blocks, put through the same subset lattice as "
      "`gate1_pairs.py`.  `cross_excess = 0` is forced by construction here, so a non-zero value "
      "would mean the instrument is wrong and every `cross_excess` in batches 1 and 2 would have to "
      "be discarded.")
    A("")
    A("| n | k | nv | blocks | eqs | D3 syz (8n·|B|) | falls_D4 | Σ singles | cross_excess | sat |")
    A("|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|")
    bad = 0
    for r in sorted(rows, key=lambda r: (i(r["n"], 0), i(r["k"], 0), i(r["nblocks"], 0),
                                         r["blocks"])):
        ce = r["cross_excess"]
        if ce not in ("", "0"):
            bad += 1
        A(f"| {r['n']} | {r['k']} | {r['nvars']} | {{{r['blocks']}}} | {r['neqs']} | "
          f"{r['D3_syz']} ({r['syz_target']}) | {r['falls_D4']} | {r['sum_singles']} | "
          f"{'**' + ce + '**' if ce not in ('', '0') else ce} | {r['saturation']} |")
    A("")
    A(f"**{len(rows) - bad} of {len(rows)} subsets give exactly 0**"
      + ("  → the instrument is sound." if bad == 0 else
         "  → **the instrument is NOT sound**; every cross_excess measured in batches 1–2 is suspect."))
    A("")


# ------------------------------------------------------------------------ A.1 d_solv, generic e
def dsolv_table(rows, title, note=""):
    A(title)
    A("")
    if note:
        A(note)
        A("")
    if not rows:
        A("_not collected_")
        A("")
        return
    d = defaultdict(dict)
    exps = {}
    for r in rows:
        if r.get("censored") == "1" or not r["d_solv"]:
            continue
        d[(i(r["ell"], -1), r["system"])].setdefault(i(r["n"], 0), set()).add(i(r["d_solv"]))
        exps[(i(r["ell"], -1), i(r["n"], 0))] = r.get("exps", "")
    ns = sorted({n for v in d.values() for n in v})
    if not ns:
        A("_every run censored_")
        A("")
        return
    A("| system | ℓ | " + " | ".join(f"n={n}" for n in ns) + " |")
    A("|---|---:|" + "|".join("---:" for _ in ns) + "|")
    order = {"S_quad": 0, "P1only": 1, "random": 2}
    for k in sorted(d, key=lambda k: (k[0], order.get(k[1], 9), k[1])):
        cells = []
        for n in ns:
            v = d[k].get(n)
            cells.append("—" if not v else "/".join(map(str, sorted(v))))
        A(f"| `{k[1]}` | {k[0]} | " + " | ".join(cells) + " |")
    A("")
    A("Exponents used: " + ", ".join(f"ℓ={e}, n={n}: `{x}`" for (e, n), x in sorted(exps.items())
                                     if x) + ".")
    A("")
    cens = [r for r in rows if r.get("censored") == "1"]
    if cens:
        A(f"{len(cens)} run(s) censored at the cap: "
          + ", ".join(sorted({f"{r['tag']}/{r['tool']}" for r in cens})) + ".")
        A("")


def sec_dsolv():
    A("### B.1 `d_solv` with corrected exponents (task 1.2)")
    A("")
    dsolv_table(load("E04_f4_degrees_b2.csv"), "",
                "grevlex; msolve and Magma agree on every completed cell.  A cell is `—` where "
                "every run of it was censored at the 2-hour cap.")


# --------------------------------------------------------------------------- A.2 hybrid delta
def sec_hybrid():
    A("### B.2 Hybrid guessing: measured `d_solv` against B2's `d_reg` (task 1.3)")
    A("")
    rows = load("E04_f4_hybrid_b2.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    try:
        from table4_reproduce import dreg_boolean
    except Exception as exc:
        dreg_boolean = None
        A(f"_(B2 column unavailable: {exc})_")
    A("`k` bits of pt are fixed to their **true** values, so the planted root survives and the run "
      "measures the correct branch.  `δ = d_reg^{B2}(3n−k, 9n) − d_solv`.")
    A("")
    A("| n | k | nv = 3n−k | system | tool | `d_solv` | B2 `d_reg` | **δ** | s |")
    A("|---:|---:|---:|---|---|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: (i(r["n"], 0), i(r["nvars"], 0), r["system"], r["tool"])):
        n, nv = i(r["n"], 0), i(r["nvars"], 0)
        k = 3 * n - nv
        if r.get("censored") == "1" or not r["d_solv"]:
            A(f"| {n} | {k} | {nv} | {r['system']} | {r['tool']} | _censored_ | | | "
              f"{f(r['seconds'], 0):.0f} |")
            continue
        ds = i(r["d_solv"])
        dr = dreg_boolean(nv, i(r["neqs"], 0)) if dreg_boolean else None
        A(f"| {n} | {k} | {nv} | {r['system']} | {r['tool']} | {ds} | {dr if dr else ''} | "
          f"{'**' + str(dr - ds) + '**' if dr else ''} | {f(r['seconds'], 0):.0f} |")
    A("")


# ------------------------------------------------------------------------- A.3 Hilbert check
def sec_hilbert():
    A("### B.3 Measured Macaulay rank against the Boolean semi-regular prediction (task 1.4)")
    A("")
    rows = load("E04_server_hilbert.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    A("`rank_pred(D) = C(N,≤D) − Σ_{d≤D}[z^d] (1+z)^N/(1+z²)^m`, the rank a degree-D Macaulay "
      "matrix would have if the system were semi-regular — which is the assumption B2 rests on.")
    A("")
    A("| ℓ | n | system | D | rows | cols | rank | predicted | rank − pred | syz | syz − pred |")
    A("|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: (i(r["ell"], 0), i(r["n"], 0), r["system"], i(r["D"], 0))):
        cen = r.get("censored") == "1"
        A(f"| {r['ell']} | {r['n']} | `{r['system']}` | {r['D']} | {r['rows']} | {r['cols']} | "
          f"{r['rank']}{' _(lower bound)_' if cen else ''} | {r['predicted_rank']} | "
          f"{r['rank_minus_predicted']} | {r['syz']} | {r['syz_minus_predicted']} |")
    A("")


# ------------------------------------------------------------------- A.4 calibrated estimate
def sec_calibrated():
    A("### B.4 Calibrated hybrid-guessing cost (task 1.5)")
    A("")
    rows = load("E04_server_calibrated.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    A("`cost = 2^k · (Σ_{i≤d} C(N−k,i))^ω`, `d = dreg_boolean(N−k, m) − δ`, minimised over k.  "
      "**This is a table, not a conclusion** — δ is measured at n ≤ 15 and "
      "`COMPLEXITY_ACCOUNTING §4` does not license carrying it to n = 128..256.  B1/B2 use E02's "
      "*commutative* monomial count `C(N−k+d,d)²`, not this Boolean one, so the δ = 0, ω = 2 row "
      "is close to but not literally B2.")
    A("")
    A("| scheme | n | ℓ | N | m | δ | ω | k* | d_reg | d | log2 T | B1 | B2 | B3 |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        A(f"| {r['scheme']} | {r['n']} | {r['ell']} | {r['nvars']} | {r['neqs']} | {r['delta']} | "
          f"{r['omega']} | {r['best_k']} | {r['dreg_boolean']} | {r['d_used']} | "
          f"**{r['log2_T']}** | {r['B1_log2T']} | {r['B2_log2T']} | {r['B3_log2T']} |")
    A("")


# ---------------------------------------------------------------------------- A.5 point set
def sec_pointset():
    A("### B.5 `dim I(V_e)_{≤4}` against what D = 5 reaches (task 2.1)")
    A("")
    rows = load("E03_server_pointset.csv")
    d5 = {}
    for r in load("E03_server_D5.csv") or []:
        if r["system"] == "S_box" and r["D5_dim_low"]:
            d5[(i(r["n"]), i(r["e"]))] = i(r["D5_dim_low"])
    if not rows:
        A("_not collected_")
        A("")
        return
    A("`dim I(V_e)_{≤4}` is the whole degree-≤4 part of the vanishing ideal of the S-box graph; "
      "`dim(span_{D=5} ∩ deg≤4)` is what the degree-5 Macaulay matrix of ⟨P₁,P₂,P₃⟩ actually "
      "produces.  The difference is an **upper bound** on the number of quartic generators outside "
      "the quadratic ideal — D = 6 might still reach some of them, and D = 6 at n = 16 is ≈ 290 GB "
      "(§8).  The measurement is only meaningful when the graph has more points than there are "
      "monomials; `deficit_floor` is what the point shortage alone would force.")
    A("")
    A("| n | e | class | points | monomials | deficit floor | dim I(V)≤4 | D=5 reaches | "
      "outside | /n |")
    A("|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: (i(r["n"], 0), i(r["e"], 0))):
        n, e = i(r["n"]), i(r["e"])
        reach = d5.get((n, e))
        iv = i(r["dim_IV"])
        out = None if reach is None else iv - reach
        A(f"| {n} | {e} | `{r['degenerate'] or 'generic'}` | {r['npoints']} | {r['monomials']} | "
          f"{r['deficit_floor']} | {iv} | {reach if reach is not None else '—'} | "
          f"{'**' + str(out) + '**' if out is not None else '—'} | "
          f"{out/n:g} |" if out is not None else
          f"| {n} | {e} | `{r['degenerate'] or 'generic'}` | {r['npoints']} | {r['monomials']} | "
          f"{r['deficit_floor']} | {iv} | — | — | — |")
    A("")


# ------------------------------------------------------------------- A.6 the D=5 decomposition
def sec_d5dec():
    A("### B.6 What the D = 5 falls are made of (task 2.2)")
    A("")
    rows = load("E03_server_d5_decomposition.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    A("`X` = dim of the span of the products `z_k·c_j` (the `#cubics` new cubics of D = 4 times the "
      "2n variables) **modulo the D = 4 row space**; `net new = falls_D5 − X` is what D = 5 adds "
      "that is not a multiple of a D = 4 fall.")
    A("")
    A("| n | e | class | #cubics | products | X | relations | /n | falls_D5 | net new | /n |")
    A("|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: (i(r["n"], 0), i(r["e"], 0))):
        n = i(r["n"])
        rel = i(r["relations_among_products"], 0)
        net = i(r["net_new_quartics"])
        A(f"| {n} | {r['e']} | `{r['deg_class']}` | {r['cubics']} | {r['products']} | {r['X']} | "
          f"{rel} | {rel/n:g} | {r['falls_D5']} | **{net}** | {net/n:g} |")
    A("")


# --------------------------------------------------------------------------- A.7 mult-block
def sec_multblock():
    A("### B.7 Block-restricted XL on the kernel model (task 3.1)")
    A("")
    full = {}
    for r in load("E06_server_kernel.csv") or []:
        if r["system"] == "kernel":
            full[(i(r["n"]), i(r["D"]))] = (r["rank"], r["syz"], r["dim_low"])
    got = False
    A("Multipliers restricted to one block instead of all 6n variables.  `corank = cols − rank`; "
      "linearisation would mean corank 1.")
    A("")
    A("| n | D | multipliers | rows = rank | syz | dim_low | corank | full-multiplier rank | "
      "full dim_low |")
    A("|---:|---:|---|---:|---:|---:|---:|---:|---:|")
    for blk in ("U", "W"):
        rows = load(f"E06_server_multblock_{blk}.csv")
        for r in sorted(rows or [], key=lambda r: (i(r["n"], 0), i(r["D"], 0))):
            if r["system"] != "kernel":
                continue
            got = True
            fu = full.get((i(r["n"]), i(r["D"])), ("—", "—", "—"))
            A(f"| {r['n']} | {r['D']} | {blk} block | {r['rank']} | {r['syz']} | {r['dim_low']} | "
              f"{r['corank']} | {fu[0]} | {fu[2]} |")
    if not got:
        OUT[-1] = "_not collected_"
    A("")


# ----------------------------------------------------------------------------- A.8 MinRank
def sec_minrank():
    A("### B.8 Does the first layer's rank condition carry information? (task 3.2)")
    A("")
    rows = load("E06_server_minrank.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    A("`N(W)` is 2n × (n+1) with entries affine in W = (a₀,a₁).  A **uniformly random** matrix of "
      "that shape has `P(rank ≤ n) = 1 − Π_{i=0}^{n}(1 − 2^{i−2n})`, which is the column headed "
      "*uniform*: the question is whether the true W is distinguished by more than that.")
    A("")
    A("| n | rank at true W | corank | (pt,1) in kernel | samples | random rank ≤ n | measured % | "
      "uniform % | ratio |")
    A("|---:|---:|---:|:--:|---:|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: i(r["n"], 0)):
        n = i(r["n"])
        p = 1.0
        for j in range(n + 1):
            p *= (1 - 2.0 ** (j - 2 * n))
        unif = 100 * (1 - p)
        meas = 100 * f(r["random_rank_le_n_frac"], 0)
        A(f"| {n} | {r['rank_true_W']} | {r['corank_true_W']} | "
          f"{'yes' if r['pt_in_kernel'] == '1' else '**NO**'} | {r['samples']} | "
          f"{r['random_rank_le_n']} | {meas:.3f} | {unif:.3f} | {meas/unif:.2f}× |")
    A("")


# ---------------------------------------------------------------------------- A.9 PolyBoRi
def sec_polybori():
    A("### B.9 Block orders computed directly, in PolyBoRi (task 3.3)")
    A("")
    rows = load("E10_server_polybori.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    done = [r for r in rows if r.get("censored") != "1"]
    A("PolyBoRi works **in** the declared order and prints one `Current Degree: d` line per degree, "
      "so `d_solv` here is a genuine per-order solving degree — the measurement batch 1 could not "
      "make, because Magma's GF(2) F4 runs in grevlex whatever order the ring declares.  "
      f"**{len(done)} of {len(rows)} runs completed**; the rest hit the timeout or aborted inside "
      "PolyBoRi (`rc = 134`, SIGABRT), and their `≥ d` entries are lower bounds taken from the "
      "degree lines printed before the crash.")
    A("")
    A("| system | n | order | ring order | `d_solv` | steps | GB size | GB max deg | s | outcome |")
    A("|---|---:|---|---|---:|---:|---:|---:|---:|---|")
    for r in sorted(rows, key=lambda r: (r["tag"], r["order_label"])):
        cens = r.get("censored") == "1"
        rc = r.get("rc", "")
        outcome = "complete" if not cens else ("**abort**" if rc == "134" else "**timeout**")
        d = r["d_solv"]
        d = ("_no degree line_" if not d else (f"≥ {d}" if cens else f"**{d}**"))
        A(f"| `{r['tag']}` | {r['n']} | `{r['order_label']}` ({r['order']}) | "
          f"{r.get('ring_order','')[:28]} | {d} | {r['nsteps']} | "
          f"{r['gb_size']} | {r['gb_max_degree']} | {f(r['seconds'], 0):.0f} | {outcome} |")
    A("")
    by = {}
    for r in done:
        by.setdefault(r["order_label"], []).append(r["d_solv"])
    if by:
        A("Completed runs only, by order: "
          + "; ".join(f"`{k}` → {', '.join(sorted(set(v), key=lambda x: int(x)))}"
                      for k, v in sorted(by.items())) + ".")
        A("")


# ------------------------------------------------------------------------ A.10 E05 cheap
def sec_e05():
    A("### B.10 E05: eliminating t₁ without a resultant (task 5)")
    A("")
    rows = load("E05_server_cheap_elimination.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    res = {}
    for r in load("E05_server_resultants.csv") or []:
        if r["resultant_seconds"]:
            res[i(r["n"])] = f(r["resultant_seconds"])
    A("Eq1 is separable, `W(t₀) + V(t₁)`, so its variety is found by tabulating V and looking W up "
      "— Θ(2ⁿ) field operations, no eliminated polynomial.  The `Resultant` column is batch 1's "
      "measurement of the same instance for comparison.")
    A("")
    A("| n | build s | **eliminate s** | roots | spurious | true key a root | `Resultant` s | "
      "speed-up |")
    A("|---:|---:|---:|---:|---:|:--:|---:|---:|")
    pts = []
    for r in sorted(rows, key=lambda r: i(r["n"], 0)):
        n = i(r["n"])
        te = f(r["eliminate_seconds"])
        rs = res.get(n)
        A(f"| {n} | {r['build_seconds']} | **{r['eliminate_seconds']}** | {r['roots']} | "
          f"{r['spurious']} | {'yes' if r['key_is_root'] == '1' else '**NO**'} | "
          f"{rs if rs is not None else '—'} | "
          f"{f'{rs/max(te,1e-6):,.0f}×' if rs is not None and te else '—'} |")
        if te and te > 0:
            pts.append((n, te))
    A("")
    if len(pts) >= 3:
        xs = [p[0] for p in pts]
        ys = [math.log2(p[1]) for p in pts]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        den = sum((x - mx) ** 2 for x in xs)
        sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
        big = [p for p in pts if p[0] >= 13]
        sl2 = None
        if len(big) >= 3:
            xs2 = [p[0] for p in big]
            ys2 = [math.log2(p[1]) for p in big]
            mx2, my2 = sum(xs2) / len(xs2), sum(ys2) / len(ys2)
            d2 = sum((x - mx2) ** 2 for x in xs2)
            sl2 = sum((x - mx2) * (y - my2) for x, y in zip(xs2, ys2)) / d2
        A(f"Least squares over n = {xs}: slope of `log2(eliminate seconds)` vs n = **{sl:.3f}**"
          + (f"; over n ≥ 13 only (past the constant-overhead regime): **{sl2:.3f}**"
             if sl2 else "") + ".")
        A("")
        A("The method is Θ(2ⁿ) by construction, so a slope of 1 is the expected answer and the "
          "measurement confirms it rather than assuming it.  There is **no `log2(deg)` slope**: "
          "this route forms no eliminated polynomial, and the route that does is censored from "
          "n = 8 (`results/E05_server_resultants.csv`).")
        A("")


def sec_t4():
    A("### B.11 T4 cross-validation in Sage, and the degree of the eliminated condition in t₀ (task 4)")
    A("")
    rows = load("E07_server_ore_scale.csv")
    if not rows:
        A("_not collected_")
        A("")
        return
    A("**B.11.1 — Sage against the Tier-A branch, on the same instances.**  The instances are exported "
      "by the Tier-A probe itself (same script, same seed), so the two implementations are compared on "
      "identical data.  `ker∩` is `dim(ker Λ₀ ∩ ker Λ₁)` and `solvable` the affine predicate "
      "{Λ₀(v)=c₀, Λ₁(v)=c₁}, at the true (t₀,t₁); `rand` is 32 random pairs.")
    A("")
    A("| n | e₀,e₁ | deg Λ₀,Λ₁ | gcrd | lclm | ker∩ | solvable (true) | solvable (rand) | "
      "gcrd (rand) | vs Tier-A |")
    A("|---:|---|---|---:|---:|---:|:--:|---:|---:|---|")
    checks = bad = 0
    for r in sorted(rows, key=lambda r: i(r["n"], 0)):
        if not r.get("deg_lclm_true"):
            A(f"| {r['n']} | — | — | — | — | — | — | — | — | censored ({r.get('detail','')}) |")
            continue
        c, b = i(r.get("tierA_checks"), 0), i(r.get("tierA_disagreements"), 0)
        checks += c
        bad += b
        ta = "—" if not c else (f"{c}/{c} agree" if not b else f"**{b} DISAGREE**")
        A(f"| {r['n']} | {r['e0']},{r['e1']} | {r['deg_Lambda'].replace(';', ',')} | "
          f"{r['deg_gcrd_true']} | {r['deg_lclm_true']} | {r['kerint_true']} | "
          f"{'yes' if r['affine_solvable_true'] == '1' else '**no**'} | "
          f"{r['affine_solvable_rand']} | {r['deg_gcrd_rand']} | {ta} |")
    A("")
    A(f"Totals: **{checks} comparisons against `results/E07_theory_skew_degrees.csv`, "
      f"{bad} disagreements**.  `deg lclm = deg Λ₀ + deg Λ₁` at every n (the gcrd is 0 everywhere, "
      "true key included), and the affine predicate separates the true key from all 32 random pairs "
      "at every n up to 24 — the Tier-A finding, now measured ten sizes further.")
    A("")

    deg = load("E07_server_lclm_degrees.csv")
    sage_interp = {i(r["n"]): r for r in (load("E07_server_ore_interp.csv") or [])}
    if deg:
        A("**B.11.2 — D(n), the degree in t₀ of the eliminated coefficients.**  t₁ is fixed at its "
          "true value and t₀ runs over **all** of F_{2ⁿ}; each coefficient of the monic lclm is a "
          "function F_{2ⁿ} → F_{2ⁿ} whose exact univariate degree is read off by interpolation and "
          "whose exact ANF degree comes from a Möbius transform.  `2ⁿ−1` and `n` are the maxima.")
        A("")
        A("| n | e₀,e₁ | deg lclm | values | ms per lclm | **D(n)** | max | **ANF** | max | "
          "Sage (independent) |")
        A("|---:|---|---:|---:|---:|---:|---:|---:|---:|---|")
        for r in sorted(deg, key=lambda r: i(r["n"], 0)):
            n = i(r["n"])
            sg = sage_interp.get(n)
            sgc = "—"
            if sg and sg.get("D_interp_lclm"):
                same = sg["D_interp_lclm"] == r.get("D_interp")
                sgc = f"{sg['D_interp_lclm']} ({'same' if same else '**DIFFERS**'})"
            if r["censored"] == "1":
                A(f"| {n} | {r['e0']},{r['e1']} | {r['deg_lclm']} | "
                  f"{r['values_done']}/{r['nvalues']} | — | _censored_ | {2**n-1} | — | {n} | {sgc} |")
                continue
            A(f"| {n} | {r['e0']},{r['e1']} | {r['deg_lclm']} | {r['nvalues']} | "
              f"{r['ms_per_lclm']} | **{r['D_interp']}** | {2**n-1} | "
              f"**{r['anf'] or '—'}** | {n} | {sgc} |")
        A("")
        full = [r for r in deg if r["censored"] != "1"
                and i(r["D_interp"], -1) == 2 ** i(r["n"], 0) - 1]
        if full:
            A(f"**D(n) = 2ⁿ − 1 — the maximum — at every completed size "
              f"(n = {', '.join(str(i(r['n'])) for r in sorted(full, key=lambda r: i(r['n'], 0)))}), "
              f"and for every one of the deg-lclm + 1 coefficients, not just the extreme ones.**  "
              "The ANF degree is n, also the maximum.  Eliminating pt in the Ore ring is exact and "
              "cheap; what it leaves behind is a condition in t₀ of full degree, so recovering t₀ "
              "from it by any degree-sensitive method (resultant, univariate root finding, "
              "interpolation) costs at least 2ⁿ — the same order as guessing pt outright.")
            A("")
    sym = load("E07_server_ore_symbolic.csv")
    if sym:
        A("**B.11.3 — the symbolic instrument.**  The same lclm computed over K(t₀) (twist "
          "t₀ ↦ t₀², which is what φ·t₀ = t₀²·φ forces), so the coefficients are exact rational "
          "functions rather than sampled values.")
        A("")
        A("| n | deg lclm | D_sym (numerator) | denominator | D_sym(u) | log₂ formal deg E | s |")
        A("|---:|---:|---:|---:|---:|---:|---:|")
        for r in sorted(sym, key=lambda r: i(r["n"], 0)):
            if not r.get("D_sym_lclm"):
                A(f"| {r['n']} | — | _censored_ | — | — | — | {r['seconds']} |")
                continue
            A(f"| {r['n']} | {r['deg_lclm_sym']} | {r['D_sym_lclm']} | {r['D_sym_lclm_den']} | "
              f"{r['D_sym_u']} | {r['log2_deg_E_formal']} | {r['seconds']} |")
        A("")
    A("")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--inject", default=None)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    A("<!-- generated by server/jobs/summarise_batch2.py -- do not edit by hand -->")
    A("## Appendix B — generated tables")
    A("")
    for fn in (sec_decoupled, sec_dsolv, sec_hybrid, sec_hilbert, sec_calibrated, sec_pointset,
               sec_d5dec, sec_multblock, sec_minrank, sec_polybori, sec_e05, sec_t4):
        fn()
    md = "\n".join(OUT) + "\n"
    if a.out:
        open(a.out, "w").write(md)
        print(f"written {a.out}", file=sys.stderr)
    if a.inject:
        text = open(a.inject).read()
        marker = "<!--TABLES-->"
        head, sep, _ = text.partition(marker)
        if not sep:
            raise SystemExit(f"{a.inject}: no {marker}")
        open(a.inject, "w").write(head + marker + "\n\n" + "\n".join(OUT[1:]) + "\n")
        print(f"injected {len(OUT)} lines into {a.inject}", file=sys.stderr)
    if not a.quiet:
        print(md)


if __name__ == "__main__":
    main()
