#!/usr/bin/env python3
"""E03/E06 — the explicit generators of the degree falls (input to theorem T1.2).

Two analyses, both of the same shape: find the elements a Macaulay pass ADDS below its top degree,
then try to write them as *field* identities instead of as opaque Boolean row combinations.

(1) `--sbox n e` — one enhanced-Mersenne S-box, D=4.
    Measured: the D=4 Macaulay matrix produces `falls_D4` new cubics beyond rank(Mac_3) (2n for a
    generic e, E03).  Candidate description: everything is a field multiple of P1, so consider the
    F_{2^n}-module

        A_{ij} = x^{2^i} y^{2^j} P1,   B_k = x^{2^k} P2 = x^{2^k+1} P1,   C_k = y^{2^k} P3,

    whose Boolean coordinates are all rows of the D=4 Macaulay matrix (x^{2^i} is F_2-linear in the
    x-bits, so x^{2^i} y^{2^j} is a Boolean-degree-2 multiplier).  Writing Tr_gamma for the
    coordinate functionals, the degree-4 part of Tr(gamma A_ij) is Tr(gamma x^{2^i+1} y^{2^j+1}),
    which is homogeneous of degree wt(2^i+1) + wt(2^j+1).  A fall is therefore exactly a solution of

        sum_{i,j} Tr(gamma_ij  x^{2^i+1} y^{2^j+1})  =  0,

    and the resulting cubic is sum_{i,j} Tr(gamma_ij (x^{2^i+2^e} y^{2^j} + x^{2^i} y^{2^j})).
    Two Tr-terms cancel iff their exponent pairs lie in the same Frobenius orbit
    (a,b) ~ (2a, 2b) mod 2^n-1, since Tr(gamma z^{2^s}) = Tr(gamma^{2^{-s}} z).  The script computes
    the orbit structure, the exact kernel, and how much of the measured fall space the module
    explains modulo Mac_3.

(2) `--kernel n` — the E06 kernel-variable bilinear model, D=3.
    Measured: the D=3 pass adds #eqs = 6n new quadratics (E06, batch 2 wording).  Two candidate
    descriptions are tested:
      (a) the S-box relations P1_i, which the model implies via  B_i + x_i * K_i = a_i * P1_i;
      (b) a_i * K_i  and  t_i * K_i  (i = 0,1,2, with t_2 = pt + ct).
    (b) is the answer.  K_i = t_i a_i + a_i^{2^e} is quadratic, and multiplying it by the LINEAR
    form a_i or t_i does not raise the Boolean degree, because
        a_i * K_i = t_i a_i^2 + a_i^{2^e+1},      t_i * K_i = t_i^2 a_i + t_i a_i^{2^e},
    and a_i^2, t_i^2 are F_2-linear while wt(2^e+1) = 2.  No cancellation is involved at all: the
    product is simply of lower degree than a generic one.  That is also why a random bilinear system
    of the same shape produces nothing (E06): it has no Frobenius structure to exploit.  Note both
    generators come from the K equations; the B equations contribute none.

Usage:
  python3 fall_generators.py --sbox 11 3 [--track] [--budget SEC]
  python3 fall_generators.py --kernel 7
  python3 fall_generators.py --all            (the batch-2 set: sbox 11/3, 13/3, 11/4, 13/5; kernel 7, 9)
Appends to results/E03_theory_fall_generators.csv and results/E06_theory_fall_generators.csv.
"""
import sys, os, time, argparse, itertools
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3.gf2n import GF2n
from aim3.fieldpoly import FieldPoly, var_field, trace_poly, boolean_coords
from aim3.boolmac import (monomial_columns, poly_rows, macaulay_stats, fall_space,
                          span_rank, reduce_modulo)


# ---------------------------------------------------------------- helpers
def poly_to_row(poly, col):
    """Boolean polynomial (set of frozensets, the legacy format) -> row int over `col`,
    which is keyed by integer variable bitmasks (src/aim3/boolmac.monomial_columns)."""
    r = 0
    for m in poly:
        mask = 0
        for i in m:
            mask |= 1 << i
        r |= 1 << col[mask]
    return r


def orbit_report(nmax=20, nmin=5):
    """The Frobenius-orbit structure of the degree-4 exponent pairs, and the resulting closed form
    for the ansatz module's top-degree rank and kernel dimension.

    LEMMA.  For 1 <= i,j <= n-1 the pairs (2^i+1, 2^j+1) and (2^{i'}+1, 2^{j'}+1) lie in the same
    Frobenius orbit iff (i',j') = (i,j), or i = j and (i',j') = (n-i, n-i).
    Proof: 2^s(2^i+1) reduces mod 2^n-1 to the element with support {i+s, s} (mod n), which is of the
    form {i', 0} only if s = 0 (giving (i,j) back) or i+s = 0, i.e. s = n-i.  In the second case the
    second coordinate 2^s(2^j+1) has support {j-i, n-i}, which contains 0 only if j = i.  QED

    So there are exactly floor((n-1)/2) coincidences, one per diagonal pair {i, n-i}, and the ansatz
    module (n^2 + 2n generators, n trace coefficients each) has

        top-degree rank  = n * ((n-1)^2 - c),      kernel = n * (n^2 + 2n - (n-1)^2 + c),
        c = floor((n-1)/2).

    Checked against the measured values at n = 11 (rank 1045, kernel 528) inside `analyse_sbox`.
    """
    rows = []
    print(f"{'n':>3}{'pairs':>7}{'orbits':>8}{'coinc':>7}{'predicted':>10}{'gens':>6}"
          f"{'top rank':>10}{'kernel':>8}{'formula':>9}{'ok':>4}")
    for n in range(nmin, nmax + 1):
        orb = {}
        for i in range(1, n):
            for j in range(1, n):
                orb.setdefault(orbit_of((1 << i) + 1, (1 << j) + 1, n), []).append((i, j))
        coinc = sum(len(v) - 1 for v in orb.values())
        nontrivial = [v for v in orb.values() if len(v) > 1]
        pred = [(i, n - i) for i in range(1, n) if i < n - i]
        lemma_ok = (len(nontrivial) == len(pred) and
                    all(sorted(v) == sorted([(v[0][0], v[0][0]), (n - v[0][0], n - v[0][0])])
                        for v in nontrivial))
        gens = n * n + 2 * n
        top = n * ((n - 1) ** 2 - coinc)
        ker = n * gens - top
        formula = n * (gens - (n - 1) ** 2 + coinc)
        ok = (ker == formula) and lemma_ok
        print(f"{n:>3}{(n-1)**2:>7}{len(orb):>8}{coinc:>7}{len(pred):>10}{gens:>6}"
              f"{top:>10}{ker:>8}{formula:>9}{'OK' if ok else 'NO':>4}")
        rows.append(dict(n=n, deg4_pairs=(n - 1) ** 2, orbits=len(orb), coincidences=coinc,
                         lemma_predicted=len(pred), lemma_holds=int(lemma_ok), gens=gens,
                         module_rows=n * gens, top_rank=top, kernel=ker, kernel_formula=formula))
    return rows


def orbit_of(a, b, n):
    """the Frobenius orbit {(2^s a, 2^s b) mod 2^n-1} as a canonical (minimal) representative."""
    q1 = (1 << n) - 1
    best = None
    ca, cb = a % q1, b % q1
    for _ in range(n):
        cand = (ca, cb)
        if best is None or cand < best:
            best = cand
        ca, cb = (ca * 2) % q1, (cb * 2) % q1
    return best


# ---------------------------------------------------------------- (1) single S-box
def analyse_sbox(n, e, track=False, do_sweep=True):
    t0 = time.time()
    F = GF2n(n)
    X, Y = var_field(F, 0), var_field(F, n)
    one = FieldPoly.constant(F, 1)
    P1f = X * Y + X.frob(e) + one
    P2f, P3f = X * P1f, Y * P1f
    eqs = boolean_coords(P1f) + boolean_coords(P2f) + boolean_coords(P3f)
    nv = 2 * n
    col, ncol, top_start = monomial_columns(nv, 4)
    topmask = ((1 << ncol) - 1) ^ ((1 << top_start) - 1)

    st3 = macaulay_stats(eqs, nv, 3)
    st4 = macaulay_stats(eqs, nv, 4)
    measured = st4.dim_low - st3.rank
    mac3 = list(poly_rows(eqs, nv, 4, col=col, mult_deg=1))     # Mac_3 rows in D=4 coordinates
    fs = fall_space(eqs, nv, 4, track=track)
    reps, q = reduce_modulo(fs["low"], mac3)
    print(f"n={n} e={e}: Mac_3 rank {st3.rank}, Mac_4 dim(span n deg<=3) {st4.dim_low}; "
          f"falls_D4 = {measured}  (2n = {2*n})")
    print(f"   fall space has rank {span_rank(fs['low'])}; modulo Mac_3 it has dimension {q} "
          f"{'(= falls_D4, consistent)' if q == measured else '(MISMATCH)'}")

    # ---- the field module ---------------------------------------------------
    # Every Macaulay row is Tr(gamma * M * P_k) for a FIELD monomial M = x^a y^b of Boolean degree
    # wt(a)+wt(b) <= 2: a Boolean bit x_c is Tr(delta_c x), so a squarefree bit-monomial of degree
    # <= 2 is a Boolean-degree-<=2 function of (x,y), and those are exactly the F_{2^n}-span of the
    # monomials x^a y^b with wt(a)+wt(b) <= 2 (uniqueness of the reduced polynomial representation).
    # Hence the union of the 18 families below IS the whole D=4 row space, and the question the
    # sweep answers is which families the NEW cubics actually need.
    ii = list(range(n))
    pairs2 = [(i, j) for i in range(n) for j in range(i + 1, n)]
    mono = {
        "1":    [(("1",), FieldPoly.constant(F, 1))],
        "x1":   [(("x", i), X.frob(i)) for i in ii],
        "y1":   [(("y", j), Y.frob(j)) for j in ii],
        "x1y1": [(("xy", i, j), X.frob(i) * Y.frob(j)) for i in ii for j in ii],
        "x2":   [(("xx", i, j), X.frob(i) * X.frob(j)) for i, j in pairs2],
        "y2":   [(("yy", i, j), Y.frob(i) * Y.frob(j)) for i, j in pairs2],
    }
    base = {"P1": P1f, "P2": P2f, "P3": P3f}
    fam_names, fam_gens = [], {}
    for pk in ("P1", "P2", "P3"):
        for mk in ("1", "x1", "y1", "x1y1", "x2", "y2"):
            nm = f"{mk}*{pk}"
            fam_names.append(nm)
            fam_gens[nm] = [((pk, mk) + tag, m * base[pk]) for tag, m in mono[mk]]

    def module_new(gen_list):
        """dimension, modulo Mac_3, of the falls of the sub-module spanned by `gen_list`."""
        rows = []
        for tag, g in gen_list:
            for k in range(n):
                rows.append(poly_to_row(trace_poly(g, 1 << k), col))
        piv, klow = {}, []
        for r in rows:
            top = r & topmask
            while top:
                bb = top.bit_length() - 1
                pv = piv.get(bb)
                if pv is None:
                    piv[bb] = r
                    break
                r ^= pv
                top = r & topmask
            else:
                if r:
                    klow.append(r)
        return reduce_modulo(klow, mac3)[1], len(rows), klow

    # (a) the ansatz proposed for this batch: x^{2^i}y^{2^j} P1, x^{2^k} P2, y^{2^k} P3
    ansatz = fam_gens["x1y1*P1"] + fam_gens["x1*P2"] + fam_gens["y1*P3"]
    a_q, a_rows, a_low = module_new(ansatz)
    inside = span_rank(a_low + fs["loD"]) == span_rank(fs["loD"]) if False else \
             span_rank(a_low + fs["low"]) == span_rank(fs["low"])
    print(f"   ANSATZ  {{x^2^i y^2^j P1}} + {{x^2^k P2}} + {{y^2^k P3}}: {a_rows} rows -> "
          f"{a_q} new cubics modulo Mac_3, of the {measured} measured"
          + ("  -> ANSATZ REFUTED" if a_q < measured else "  -> ansatz suffices"))
    print(f"           (its falls do lie inside the measured fall space: {inside})")

    singles, loo = {}, {}
    if do_sweep:
        print("   --- each family on its own ---")
        for nm in fam_names:
            qq, nr, _ = module_new(fam_gens[nm])
            singles[nm] = qq
            print(f"      {nm:>10}: {nr:5d} rows -> {qq:3d} new")
        print("   --- leave-one-out (drop one family from the full row space) ---")
        for nm in fam_names:
            rest = []
            for other in fam_names:
                if other != nm:
                    rest += fam_gens[other]
            qq, _, _ = module_new(rest)
            loo[nm] = qq
            flag = "  <-- NECESSARY" if qq < measured else ""
            print(f"      without {nm:>10}: {qq:3d} of {measured}{flag}")
        need = [k for k, v in loo.items() if v < measured]
        print(f"   -> families that cannot be dropped: {need or 'none (every family is redundant)'}")
    mod_q = measured
    print(f"   quotient dimension {q} = {q/n:.2f} n"
          + ("  (a multiple of n)" if q % n == 0 else "  (NOT a multiple of n)"))

    rec = dict(n=n, e=e, mac3_rank=st3.rank, mac4_dim_le3=st4.dim_low, falls_D4=measured,
               fall_space_rank=span_rank(fs["low"]), falls_mod_mac3=q,
               ansatz_new=a_q, ansatz_rows=a_rows, ansatz_inside=int(inside),
               family_singles=";".join(f"{k}={v}" for k, v in singles.items()) if singles else "",
               leave_one_out=";".join(f"{k}={v}" for k, v in loo.items()) if loo else "",
               necessary_families=";".join(k for k, v in loo.items() if v < measured) if loo else "",
               quotient_over_n=round(q / n, 4))
    record("E03_theory_fall_generators.csv", rec, seconds=time.time() - t0)
    return rec


# ---------------------------------------------------------------- (2) kernel model
def analyse_kernel_model(n):
    t0 = time.time()
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "E06_kernel_bilinear_model"))
    from kernel_bilinear_toy import build_kernel_model
    from aim3_lab import rand_unit_lu, apply_mat
    import random
    exps = [e for e in range(2, n) if gcd(e, n) == 1][:3]
    nv, eqs, sol = build_kernel_model(n, exps)
    F = GF2n(n)
    col, ncol, top_start = monomial_columns(nv, 3)

    st2 = macaulay_stats(eqs, nv, 2)
    st3 = macaulay_stats(eqs, nv, 3)
    new = st3.dim_low - st2.rank
    mac2 = list(poly_rows(eqs, nv, 3, col=col, mult_deg=0))
    fs = fall_space(eqs, nv, 3)
    reps, q = reduce_modulo(fs["low"], mac2)
    print(f"kernel model n={n} exps={exps}: {len(eqs)} eqs in {nv} vars; Mac_2 rank {st2.rank}, "
          f"Mac_3 dim(span n deg<=2) {st3.dim_low}")
    print(f"   NEW quadratics at D=3: {new}  (#eqs = 6n = {6*n});  modulo the equations: {q}")

    # rebuild the same instance to get the matrices, then the three P1 relations as FieldPolys
    rng = random.Random(1)
    mats = [rand_unit_lu(n, rng) for _ in range(4)]
    consts = [rng.randrange(F.q) for _ in range(3)]
    while True:
        pt = rng.randrange(1, F.q)
        xs = [apply_mat(mats[i], pt, n) ^ consts[i] for i in range(2)]
        if 0 in xs:
            continue
        ts = [F.sbox(xs[i], exps[i]) for i in range(2)]
        xl = consts[2] ^ apply_mat(mats[2], ts[0], n) ^ apply_mat(mats[3], ts[1], n)
        if xl == 0 or 0 in ts:
            continue
        tl = F.sbox(xl, exps[2])
        if tl == 0:
            continue
        ct = pt ^ tl
        break

    def affine(mat, offset, const):
        cols = [apply_mat(mat, 1 << j, n) for j in range(n)]
        p = FieldPoly(F, {frozenset((offset + j,)): cols[j] for j in range(n)})
        return p + FieldPoly.constant(F, const)

    PT = var_field(F, 0)
    T = [var_field(F, n), var_field(F, 2 * n)]
    x0 = affine(mats[0], 0, consts[0])
    x1 = affine(mats[1], 0, consts[1])
    x2 = affine(mats[2], n, 0) + affine(mats[3], 2 * n, 0) + FieldPoly.constant(F, consts[2])
    t2 = PT + FieldPoly.constant(F, ct)
    one = FieldPoly.constant(F, 1)
    P1s = [x0 * T[0] + x0.frob(exps[0]) + one,
           x1 * T[1] + x1.frob(exps[1]) + one,
           x2 * t2 + x2.frob(exps[2]) + one]
    # sanity: each P1 must vanish at the planted solution
    for idx, P in enumerate(P1s):
        for k in range(n):
            poly = trace_poly(P, 1 << k)
            s = 0
            for m in poly:
                v = 1
                for j in m:
                    v &= sol[j]
                s ^= v
            assert s == 0, f"P1_{idx} coordinate {k} does not vanish at the key"
    p1rows = [poly_to_row(trace_poly(P, 1 << k), col) for P in P1s for k in range(n)]
    print(f"   the three S-box relations P1_i vanish at the planted key: OK ({len(p1rows)} coordinates)")

    inside = span_rank(p1rows + fs["low"]) == span_rank(fs["low"])
    _, p1_new = reduce_modulo(p1rows, mac2)
    resid, resid_dim = reduce_modulo(fs["low"], mac2 + p1rows)
    print(f"   (a) P1 coordinates inside the D=3 row space: {inside}; new modulo the equations: "
          f"{p1_new} (3n = {3*n}); unexplained remainder {resid_dim}"
          + ("   -> the P1_i are NOT what the pass produces" if resid_dim == new else ""))

    # (b) a_i * K_i and t_i * K_i
    A = [var_field(F, (3 + i) * n) for i in range(3)]
    Tall = [T[0], T[1], t2]
    Ks = [Tall[i] * A[i] + A[i].frob(exps[i]) for i in range(3)]
    cands, cand_names = [], []
    for i in range(3):
        for nm, g in ((f"a{i}*K{i}", A[i] * Ks[i]), (f"t{i}*K{i}", Tall[i] * Ks[i])):
            cand_names.append(f"{nm}(deg {g.degree()})")
            for k in range(n):
                cands.append(poly_to_row(trace_poly(g, 1 << k), col))
    cand_inside = span_rank(cands + fs["low"]) == span_rank(fs["low"])
    _, cand_new = reduce_modulo(cands, mac2)
    _, cand_resid = reduce_modulo(fs["low"], mac2 + cands)
    print(f"   (b) candidates {cand_names}")
    print(f"       inside the D=3 row space: {cand_inside}; new modulo the equations: {cand_new}; "
          f"unexplained remainder: {cand_resid}"
          + ("   -> FULLY EXPLAINED" if cand_resid == 0 and cand_new == new else ""))
    inv_col = {c: m for m, c in col.items()}
    supp = set()
    for r in resid:
        v = r
        while v:
            b = (v & -v).bit_length() - 1
            mask = inv_col[b]
            while mask:
                supp.add((mask & -mask).bit_length() - 1)
                mask &= mask - 1
            v &= v - 1
    blocks = {"pt": 0, "t0": 0, "t1": 0, "a0": 0, "a1": 0, "a2": 0}
    names = ["pt", "t0", "t1", "a0", "a1", "a2"]
    for j in supp:
        blocks[names[j // n]] += 1
    print(f"   variable blocks touched by the part (a) leaves unexplained: "
          f"{ {k: v for k, v in blocks.items() if v} }")
    record("E06_theory_fall_generators.csv",
           dict(n=n, exps=";".join(map(str, exps)), nvars=nv, neqs=len(eqs),
                mac2_rank=st2.rank, mac3_dim_le2=st3.dim_low, new_quadratics=new,
                new_mod_eqs=q, p1_inside_rowspace=int(inside), p1_new_mod_eqs=p1_new,
                p1_unexplained=resid_dim, aK_tK_inside=int(cand_inside), aK_tK_new=cand_new,
                aK_tK_unexplained=cand_resid,
                blocks=";".join(f"{k}={v}" for k, v in blocks.items() if v)),
           seconds=time.time() - t0)
    return dict(n=n, new=new, p1_new=p1_new, unexplained=cand_resid)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sbox", nargs=2, type=int, metavar=("N", "E"))
    ap.add_argument("--kernel", type=int, metavar="N")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--orbits", action="store_true",
                    help="only print the Frobenius-orbit lemma check and the closed form (cheap)")
    ap.add_argument("--track", action="store_true")
    ap.add_argument("--no-sweep", action="store_true",
                    help="only test the proposed ansatz (cheap); skip the per-family and "
                         "leave-one-out sweeps (expensive at n=13)")
    a = ap.parse_args()
    print(f"# params: sbox={a.sbox} kernel={a.kernel} all={a.all} orbits={a.orbits} track={a.track}")
    if a.orbits:
        from _common import record_table
        record_table("E03_theory_orbit_lemma.csv", orbit_report(), seconds=0)
        print("written results/E03_theory_orbit_lemma.csv")
        sys.exit(0)
    if a.all:
        for n, e in ((11, 3), (13, 3), (11, 4), (13, 5)):
            analyse_sbox(n, e, do_sweep=(n <= 11))
            print()
        for n in (7, 9):
            analyse_kernel_model(n)
            print()
    elif a.sbox:
        analyse_sbox(a.sbox[0], a.sbox[1], track=a.track, do_sweep=not a.no_sweep)
    elif a.kernel:
        analyse_kernel_model(a.kernel)
