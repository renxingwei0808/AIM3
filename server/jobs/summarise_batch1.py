#!/usr/bin/env python3
"""Turn the collected server-batch-1 CSVs into the Markdown tables of the batch report.

Reads only what `server/jobs/collect.py` has put into `results/`; writes Markdown to stdout (and,
with `--out`, to a file).  Every table is generated, never hand-typed, so the report and the CSVs
cannot drift apart.

Sections, in report order:
  1  backend alignment            server/jobs/verify_backend.py JSON, if present
  2  Gate 1                       E04_server_gate1.csv, E04_server_gate1_pairs.csv,
                                  E12_server_first_layer.csv
  3  the fall-count rule at n=17..20   E03_server_n17_20.csv vs aim3.structures.predicted_falls_D4
  4  E12 / E06 extension          E12_server_cross_relations*.csv, E06_server_kernel.csv
  5  d_solv                       E04_f4_degrees.csv
  6  block orders (E10)           E10_server_block_orders.csv
  7  E05 resultants               E05_server_resultants.csv

Usage: python3 server/jobs/summarise_batch1.py [--out results/REPORT_server_batch1_tables.md]
"""
import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from aim3.structures import predicted_falls_D4, degenerate_class     # noqa: E402

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


def i(v, default=None):
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def f(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# ------------------------------------------------------------------ 1. backend alignment
def sec_backend():
    A("### A.1 Backend alignment")
    A("")
    p = os.path.join(RES, "server", "backend_align.json")
    if not os.path.exists(p):
        A("_`backend_align.json` was not collected; see `server/ENVIRONMENT.md`._")
        A("")
        return
    d = json.load(open(p))
    checks = d["checks"]
    bad = [c for c in checks if not c["ok"]]
    A(f"`server/jobs/verify_backend.py`: **{len(checks)-len(bad)}/{len(checks)} comparisons agreed** "
      f"→ {'**BACKEND ALIGNED**' if d['ok'] else '**NOT ALIGNED**'}.")
    A("")
    A("| check | quantity | python | sage | agree |")
    A("|---|---|---:|---:|:--:|")
    for c in checks:
        A(f"| {c['check']} | `{c['quantity']}` | {c['python']} | {c['sage']} | "
          f"{'yes' if c['ok'] else '**NO**'} |")
    A("")


# ------------------------------------------------------------------------- 2. Gate 1
def sec_gate1():
    A("### A.2 Gate 1 — is the D=4 degree fall a direct sum over the S-boxes?")
    A("")
    rows = load("E04_server_gate1.csv")
    if rows:
        A("#### A.2.1 `full_system_falls.py --per-block` (S_quad)")
        A("")
        A("`saturation = dim(span_{D=4} ∩ deg≤3) / C(nv,≤3)`.  The batch rule admits a run into the "
          "verdict only when saturation < 0.9; the other rows are listed but marked.")
        A("")
        A("| ℓ | n | seed | exps | nv | falls_D4 | per-block falls | Σ | cross_excess | saturation | "
          "in verdict |")
        A("|---:|---:|---:|---|---:|---:|---|---:|---:|---:|:--:|")
        for r in sorted(rows, key=lambda r: (i(r["ell"], 0), i(r["n"], 0), i(r["seed"], 0))):
            if r["system"] != "S_quad" or not r["cross_excess"]:
                continue
            sat = f(r["saturation"])
            per = r["per_block_falls"].split(";") if r["per_block_falls"] else []
            A(f"| {r['ell']} | {r['n']} | {r['seed']} | {r['exps']} | {r['nvars']} | "
              f"{r['falls_D4']} | {' + '.join(per)} | {sum(int(x) for x in per) if per else ''} | "
              f"**{r['cross_excess']}** | {sat:.4f} | "
              f"{'yes' if sat is not None and sat < 0.9 else 'no (saturated)'} |")
        A("")
        A("Controls on the same instances (no per-block decomposition is defined for them):")
        A("")
        A("| ℓ | n | seed | system | falls_D4 | saturation |")
        A("|---:|---:|---:|---|---:|---:|")
        for r in sorted(rows, key=lambda r: (i(r["ell"], 0), i(r["n"], 0), i(r["seed"], 0))):
            if r["system"] == "S_quad" or not r["falls_D4"]:
                continue
            A(f"| {r['ell']} | {r['n']} | {r['seed']} | {r['system']} | {r['falls_D4']} | "
              f"{r['saturation']} |")
        A("")

    rows = load("E04_server_gate1_pairs.csv")
    if rows:
        A("#### A.2.2 The subset lattice, with the random control (`server/jobs/gate1_pairs.py`)")
        A("")
        A("`cross_excess(B) = falls(B) − Σ_{i∈B} falls({i})` for every subset B of the S-box "
          "blocks, in the same nv variables.  The `random` rows are the *same shape* cut into "
          "blocks at the same positions: a non-zero `cross_excess` says something about AIM3 only "
          "if the control's is zero.")
        A("")
        A("| ℓ | n | exps | system | blocks | eqs | D3 syz | falls_D4 | Σ singles | cross_excess | "
          "sat |")
        A("|---:|---:|---|---|---|---:|---:|---:|---:|---:|---:|")
        for r in sorted(rows, key=lambda r: (i(r["ell"], 0), i(r["n"], 0), r["exps"], r["system"],
                                             i(r["nblocks"], 0), r["blocks"])):
            A(f"| {r['ell']} | {r['n']} | {r['exps']} | {r['system']} | {{{r['blocks']}}} | {r['neqs']} | "
              f"{r['D3_syz']} | {r['falls_D4']} | {r['sum_singles']} | "
              f"{'**'+r['cross_excess']+'**' if r['cross_excess'] not in ('', '0') else r['cross_excess']} | "
              f"{r['saturation']} |")
        A("")

    rows = load("E12_server_first_layer.csv")
    if rows:
        A("#### A.2.3 E12 first-layer system (two S-boxes in 3n variables)")
        A("")
        A("| n | exps | nv | mac3 rank | syz (16n) | falls_D4 | S-box 0 | S-box 1 | cross_excess |")
        A("|---:|---|---:|---:|---:|---:|---:|---:|---:|")
        for r in sorted(rows, key=lambda r: i(r["n"], 0)):
            A(f"| {r['n']} | {r['exps']} | {r['nvars']} | {r['mac3_rank']} | "
              f"{r['mac3_syz']} ({r['syz_target']}) | {r['falls_D4']} | {r['falls_sbox0']} | "
              f"{r['falls_sbox1']} | **{r['cross_excess']}** |")
        A("")


# ---------------------------------------------------------- 3. the fall rule outside its fit
def sec_e03():
    A("### A.3 The calibrated fall-count rule at n = 17..20 (outside the range it was fitted on)")
    A("")
    rows = load("E03_server_n17_20.csv")
    if not rows:
        A("_no rows collected yet_")
        A("")
        return
    A("`predicted_falls_D4` (`src/aim3/structures.py`) was calibrated on 36 runs with n ≤ 16.  "
      "Every admissible e (gcd(e,n)=1, 2 ≤ e ≤ n−2) is run here; `censored=1` rows carry lower "
      "bounds and are excluded from the mismatch count.")
    A("")
    A("| n | e | m = e⁻¹ | syz₃ | = 8n | falls_D4 | falls/n | predicted | class | match | "
      "outside ideal | dim I(V)≤3 | s |")
    A("|---:|---:|---:|---:|:--:|---:|---:|---:|---|:--:|---:|---:|---:|")
    mism, done, censored = [], 0, 0
    for r in sorted(rows, key=lambda r: (i(r["n"], 0), i(r["e"], 0))):
        n, e = i(r["n"]), i(r["e"])
        if r.get("censored") == "1" or not r["falls_D4"]:
            censored += 1
            A(f"| {n} | {e} | {pow(e,-1,n)} | {r['syz3']} | "
              f"{'yes' if i(r['syz3']) == 8*n else '**NO**'} | _censored_ | | | | | | | "
              f"{r['seconds']} |")
            continue
        done += 1
        fa, pred = i(r["falls_D4"]), predicted_falls_D4(n, e)
        cls = ";".join(degenerate_class(n, e)) or "generic"
        ok = fa == pred
        if not ok:
            mism.append((n, e, fa, pred, cls))
        A(f"| {n} | {e} | {pow(e,-1,n)} | {r['syz3']} | "
          f"{'yes' if i(r['syz3']) == 8*n else '**NO**'} | {fa} | {fa/n:.2f} | {pred} | `{cls}` | "
          f"{'yes' if ok else '**NO**'} | {r['outside_ideal']} | {r['IV_deg3']} | {r['seconds']} |")
    A("")
    ns = sorted({i(r["n"]) for r in rows})
    A(f"**{len(mism)} mismatch(es) in {done} completed runs** (n = {ns}, {censored} censored).")
    if mism:
        A("")
        A("| n | e | measured | predicted | class |")
        A("|---:|---:|---:|---:|---|")
        for n, e, fa, pred, cls in mism:
            A(f"| {n} | {e} | {fa} | {pred} | `{cls}` |")
    A("")
    bad = [(i(r["n"]), i(r["e"])) for r in rows if r["syz3"] and i(r["syz3"]) != 8 * i(r["n"])]
    A(f"syz₃ = 8n in {'every' if not bad else 'all but ' + str(bad)} run.")
    A("")


# ------------------------------------------------------------------ 4. E12 / E06 extension
def sec_e12_e06():
    A("### A.4 E12 and E06 extensions")
    A("")
    rows = load("E12_server_cross_relations_deg4.csv")
    if rows:
        A("#### A.4.1 `cross_relations.py n 4` — quartic relations on (t₀,t₁)")
        A("")
        A("| n | exps | deg | points | dim on (t₀,t₁) | monomials |")
        A("|---:|---|---:|---:|---:|---:|")
        for r in rows:
            A(f"| {r['n']} | {r['exps']} | {r['deg']} | {r['npoints']} | {r['dim_t0t1']} | "
              f"{r['mono_t0t1']} |")
        A("")
    rows = load("E12_server_cross_relations_explained_deg4.csv")
    if rows:
        A("#### A.4.2 `explain_n16.py --deg 4` — the crude (and at degree 4 incorrect) test")
        A("")
        A("| n | exps | deg | dim total | dim(t₀ only) | dim(t₁ only) | wt(−r₀) | wt(−r₁) | "
          "genuine cross |")
        A("|---:|---|---:|---:|---:|---:|---:|---:|---:|")
        for r in rows:
            A(f"| {r['n']} | {r['exps']} | {r.get('deg','')} | {r['dim_cubic']} | "
              f"{r['dim_t0_only']} | {r['dim_t1_only']} | {r['wt_minus_r0']} | {r['wt_minus_r1']} | "
              f"**{r['cross_cubics']}** |")
        A("")
    rows = load("E12_server_criteria_subtract.csv")
    if rows:
        A("#### A.4.3 Subtracting the *ideal* generated by the image criteria "
          "(`server/jobs/e12_criteria_subtract.py`)")
        A("")
        A("At degree 4 it is not enough to ask whether a relation is supported on t₀ only or on t₁ "
          "only: if g = Tr(t₁^−r₁) is a cubic relation then every product `g·t₀ᵢ` is a quartic "
          "relation touching both blocks that is *not* a cross relation.  `dim C` is the "
          "degree-≤deg part of the ideal generated by the criteria; C ⊆ R, so `dim R − dim C` is "
          "the number of genuine cross-S-box relations.")
        A("")
        A("| n | deg | exps | points | monomials | deg g₀ (wt) | deg g₁ (wt) | dim R | dim C | "
          "genuine cross |")
        A("|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|")
        for r in sorted(rows, key=lambda r: (i(r["n"], 0), i(r["deg"], 0))):
            A(f"| {r['n']} | {r['deg']} | {r['exps']} | {r['npoints']} | {r['monomials']} | "
              f"{r['crit0_degree']} ({r['crit0_wt']}) | {r['crit1_degree']} ({r['crit1_wt']}) | "
              f"{r['dim_R']} | {r['dim_C']} | **{r['genuine_cross']}** |")
        A("")

    rows = load("E06_server_kernel.csv")
    if rows:
        A("#### A.4.4 E06 kernel bilinear model and its random-bilinear control")
        A("")
        A("| n | system | D | nv | eqs | rank | syz | dim_low | falls from D−1 | 2·#eqs | dim_low/n |")
        A("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in sorted(rows, key=lambda r: (i(r["n"], 0), r["system"], i(r["D"], 0))):
            A(f"| {r['n']} | {r['system']} | {r['D']} | {r['nvars']} | {r['neqs']} | {r['rank']} | "
              f"{r['syz']} | {r['dim_low']} | {r['falls_from_prev_D']} | {r['twice_neqs']} | "
              f"{r['dim_low_over_n']} |")
        A("")


# --------------------------------------------------------------------------- 5. d_solv
def sec_dsolv():
    A("### A.5 Real solving degree `d_solv` (F4 step degrees)")
    A("")
    rows = load("E04_f4_degrees.csv")
    if not rows:
        A("_no rows collected yet_")
        A("")
        return
    grev = [r for r in rows if r["order"].startswith("grevlex")]
    A("#### A.5.1 `d_solv` vs n, grevlex, both tools")
    A("")
    by = defaultdict(dict)
    for r in grev:
        if r.get("censored") == "1" or not r["d_solv"]:
            continue
        key = (i(r["ell"], -1), r["system"], i(r["seed"], -1))
        by[key][(i(r["n"], 0), r["tool"])] = i(r["d_solv"])
    ns = sorted({n for v in by.values() for n, _ in v})
    A("| ℓ | system | seed | " + " | ".join(f"n={n}" for n in ns) + " |")
    A("|---:|---|---:|" + "|".join("---:" for _ in ns) + "|")
    for key in sorted(by, key=lambda k: (k[0], k[1], k[2])):
        ell, system, seed = key
        cells = []
        for n in ns:
            ms, mg = by[key].get((n, "msolve")), by[key].get((n, "magma"))
            if ms is None and mg is None:
                cells.append("–")
            elif ms is None or mg is None:
                cells.append(str(ms if ms is not None else mg))
            elif ms == mg:
                cells.append(str(ms))
            else:
                cells.append(f"**{ms}/{mg}**")
        A(f"| {ell} | {system} | {seed} | " + " | ".join(cells) + " |")
    A("")
    A("A cell `a/b` in bold would mean msolve and Magma disagreed; a plain number means they "
      "agreed (or only one tool ran).")
    A("")
    A("#### A.5.2 Largest F4 matrix and cost, grevlex")
    A("")
    A("| n | ℓ | system | seed | tool | d_solv | steps | max rows | max cols | s | peak RSS (MB) |")
    A("|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|")
    for r in sorted(grev, key=lambda r: (i(r["n"], 0), i(r["ell"], 0), r["system"],
                                         i(r["seed"], 0), r["tool"])):
        peak = i(r["peak_rss_kb"])
        peaks = f"{peak/1024:.0f}" if peak else ""
        A(f"| {r['n']} | {r['ell']} | {r['system']} | {r['seed']} | {r['tool']} | "
          f"{r['d_solv'] or '_censored_'} | {r['nsteps']} | {r['max_rows']} | {r['max_cols']} | "
          f"{f(r['seconds'], 0):.1f} | {peaks} |")
    A("")


# ---------------------------------------------------------------------- 6. block orders
def sec_blocks():
    rows = load("E10_server_block_orders.csv")
    if not rows:
        return
    A("### A.6 Block orders (E10)")
    A("")
    A("Magma 2.24's `BooleanPolynomialRing` takes an order *name* only, so the block orders run on "
      "`PolynomialRing(GF(2), nv, \"elim\", S)` with the field equations added by hand.  That ring "
      "is slower than the Boolean one, so the baseline for these rows is `genref` "
      "(= grevlex on the *same* general ring), not the `grevlex` of §5.")
    A("")
    A("| n | ℓ | system | seed | order | ring | d_solv | steps | max rows | s |")
    A("|---:|---:|---|---:|---|---|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: (i(r["n"], 0), i(r["ell"], 0), r["system"],
                                         i(r["seed"], 0), r["order_label"])):
        A(f"| {r['n']} | {r['ell']} | {r['system']} | {r['seed']} | `{r['order_label']}` "
          f"({r['order']}) | {r.get('ring','')} | {r['d_solv'] or '_censored_'} | {r['nsteps']} | "
          f"{r['max_rows']} | {f(r['seconds'], 0):.1f} |")
    A("")


# ------------------------------------------------------------------------ 7. E05 resultants
def sec_e05():
    rows = load("E05_server_resultants.csv")
    if not rows:
        return
    A("### A.7 E05 — resultant of the two-valued-inverse reduced system")
    A("")
    A("| n | exps | ε | Eq1 terms | Eq2 terms | deg_t1(Eq1) | deg_t1(Eq2) | resultant s | "
      "deg | terms | roots | key is a root | factors | peak RSS (GB) | censored |")
    A("|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|:--:|---:|---:|:--:|")
    pts = []
    for r in sorted(rows, key=lambda r: i(r["n"], 0)):
        peak = i(r["peak_rss_kb"], 0) / 2**20
        A(f"| {r['n']} | {r['exps']} | {r['eps']} | {r['eq1_terms']} | {r['eq2_terms']} | "
          f"{r['eq1_deg_t1']} | {r['eq2_deg_t1']} | {r['resultant_seconds']} | "
          f"{r['resultant_degree']} | {r['resultant_terms']} | {r['n_roots']} | "
          f"{r['key_is_root_of_resultant']} | {r['n_factors']} | {peak:.1f} | "
          f"{'yes' if r.get('censored') == '1' else 'no'} |")
        if r.get("censored") != "1":
            t, d = f(r["resultant_seconds"]), i(r["resultant_degree"], 0)
            if t and t > 0 and d > 0:
                pts.append((i(r["n"]), t, d))
    A("")
    if len(pts) >= 2:
        xs = [p[0] for p in pts]

        def slope(ys):
            mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
            den = sum((x - mx) ** 2 for x in xs)
            return None if den == 0 else sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
        st = slope([math.log2(p[1]) for p in pts])
        sd = slope([math.log2(p[2]) for p in pts])
        A(f"Least squares over the completed sizes n = {xs}:")
        A("")
        A(f"* slope of `log2(resultant seconds)` vs n = **{st:.3f}**")
        A(f"* slope of `log2(resultant degree)` vs n = **{sd:.3f}**")
        A("")
        A(f"ROADMAP §3.3 closes track T3 when the time slope is ≥ 2: "
          f"**{'met' if st >= 2 else 'not met'}**"
          + ("." if len(xs) >= 4 else
             "  — but with fewer than 4 completed sizes this is not yet an extrapolation "
             "(COMPLEXITY_ACCOUNTING §4)."))
        A("")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--inject", default=None,
                    help="replace the <!--TABLES--> marker in this file with the generated tables "
                         "(the marker is re-inserted, so the injection is idempotent)")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    A("<!-- generated by server/jobs/summarise_batch1.py -- do not edit by hand -->")
    A("## Appendix — generated tables")
    A("")
    sec_backend()
    sec_gate1()
    sec_e03()
    sec_e12_e06()
    sec_dsolv()
    sec_blocks()
    sec_e05()
    md = "\n".join(OUT) + "\n"
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(md)
        print(f"written {a.out}", file=sys.stderr)
    if a.inject:
        text = open(a.inject).read()
        marker = "<!--TABLES-->"
        head, sep, _ = text.partition(marker)
        if not sep:
            raise SystemExit(f"{a.inject}: no {marker} marker to inject into")
        body = "\n".join(OUT[1:])                 # drop only the generated-file comment
        with open(a.inject, "w") as fh:
            fh.write(head + marker + "\n\n" + body + "\n")
        print(f"injected {len(OUT)} lines into {a.inject}", file=sys.stderr)
    if not a.quiet:
        print(md)


if __name__ == "__main__":
    main()
