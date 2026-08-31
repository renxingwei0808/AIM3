#!/usr/bin/env python3
"""E01 — does the linear-layer DISTRIBUTION change the algebra?  (control for E03/E04/E12)

Every toy instance in batch 1 used the spec's unit-LU matrices (A = L*U, L unit lower, U unit upper,
both from the SHAKE tape).  `docs/ROADMAP.md` section 4.3 requires the matrix distribution to be a
control variable, not an assumption.  This script rebuilds the E12 first-layer system -- two
first-layer S-boxes, 6n quadratics {P1,P2,P3} in the 3n variables (pt,t0,t1) -- with three
distributions and compares the D=3 rank / syzygies and the D=4 degree falls:

    unit-LU        A = L*U as in AIM3 Algorithm 4                       (the spec distribution)
    GL_n           uniform random invertible n x n matrix over F_2
    density        random matrix with the SAME expected density as unit-LU, conditioned invertible

`build_first_layer_with` mirrors `experiments/E12_cross_sbox_relations/first_layer_ideal.py`
exactly; with `--check` it is verified to reproduce E12's unit-LU numbers bit for bit, so any
difference below is due to the distribution and not to a different construction.

Usage: python3 matrix_distribution_control.py N [--seeds 1,2] [--no-D4] [--check] [--budget SEC]
Appends to results/E01_theory_matrix_distribution.csv.
"""
import sys, os, time, random, argparse
from math import gcd, comb
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import record
from aim3_lab import F, sbox, rand_unit_lu, apply_mat
from aim3_macaulay import AffExpr, prod_coords, lin_coords, add_polys
from aim3.boolmac import macaulay_stats, MacaulayTimeout


# ---------------------------------------------------------------- matrix distributions
def mat_density(A):
    n = len(A)
    return sum(sum(row) for row in A) / (n * n)


def is_invertible(A, n):
    rows = [sum(A[i][j] << j for j in range(n)) for i in range(n)]
    piv = {}
    for r in rows:
        while r:
            b = r.bit_length() - 1
            if b in piv:
                r ^= piv[b]
            else:
                piv[b] = r
                break
        else:
            return False
    return len(piv) == n


def rand_gl(n, rng):
    """uniform random invertible matrix (rejection; the density of GL_n is ~0.289 for large n)."""
    while True:
        A = [[rng.getrandbits(1) for _ in range(n)] for _ in range(n)]
        if is_invertible(A, n):
            return A


def rand_density(n, rng, p):
    """random matrix with entry density p, conditioned invertible."""
    scale = 10 ** 6
    thr = int(p * scale)
    while True:
        A = [[1 if rng.randrange(scale) < thr else 0 for _ in range(n)] for _ in range(n)]
        if is_invertible(A, n):
            return A


DISTRIBUTIONS = ("unit-LU", "GL_n", "density")


def make_matrices(dist, n, rng, p=None):
    if dist == "unit-LU":
        return [rand_unit_lu(n, rng), rand_unit_lu(n, rng)]
    if dist == "GL_n":
        return [rand_gl(n, rng), rand_gl(n, rng)]
    if dist == "density":
        return [rand_density(n, rng, p), rand_density(n, rng, p)]
    raise ValueError(dist)


# ---------------------------------------------------------------- the first-layer system
def build_first_layer_with(n, exps2, A, b):
    """The two first-layer S-boxes in 3n variables (pt, t0, t1) with the GIVEN matrices/constants.
    Identical in structure to E12's build_first_layer; only the matrix source differs."""
    fld = F(n)
    nv = 3 * n

    def colimgs(M):
        return [apply_mat(M, 1 << j, n) for j in range(n)]

    eqsAll = []
    for i in range(2):
        ci = colimgs(A[i])
        u = AffExpr(fld, {j: ci[j] for j in range(n)}, b[i])
        t = AffExpr(fld, {(i + 1) * n + k: 1 << k for k in range(n)}, 0)
        t2 = AffExpr(fld, {(i + 1) * n + k: fld.pw(1 << k, 2) for k in range(n)}, 0)
        e = exps2[i]
        P1 = add_polys(prod_coords(u, t), lin_coords(u.frob(e)), lin_coords(AffExpr(fld, {}, 1)))
        P2 = add_polys(prod_coords(u.frob(1), t), prod_coords(u.frob(e), u), lin_coords(u))
        P3 = add_polys(prod_coords(u, t2), prod_coords(u.frob(e), t), lin_coords(t))
        eqsAll += P1 + P2 + P3
    return nv, eqsAll


def valid_key(n, exps2, A, b, rng):
    fld = F(n)
    for _ in range(10000):
        pt = rng.randrange(1, fld.q)
        xs = [apply_mat(A[i], pt, n) ^ b[i] for i in range(2)]
        if 0 in xs:
            continue
        ts = [sbox(fld, xs[i], exps2[i]) for i in range(2)]
        sol = [0] * (3 * n)
        for k in range(n):
            sol[k] = (pt >> k) & 1
            for i in range(2):
                sol[(i + 1) * n + k] = (ts[i] >> k) & 1
        return sol
    raise RuntimeError("no valid key found")


def run(n, seed, dists, do_D4, budget):
    fld = F(n)
    es = [e for e in range(3, n - 1) if gcd(e, n) == 1][:2]      # identical to E12
    nv = 3 * n
    low_cols = sum(comb(nv, d) for d in range(4))
    # measure the unit-LU density once, to parameterise the "same density" control
    dens = sum(mat_density(rand_unit_lu(n, random.Random(1000 + i))) for i in range(20)) / 20
    out = []
    for dist in dists:
        t0 = time.time()
        rng = random.Random(seed)
        A = make_matrices(dist, n, rng, p=dens)
        b = [rng.randrange(fld.q), rng.randrange(fld.q)]
        _, eqs = build_first_layer_with(n, es, A, b)
        sol = valid_key(n, es, A, b, rng)
        for p in eqs:
            s = 0
            for m in p:
                v = 1
                for j in m:
                    v &= sol[j]
                s ^= v
            assert s == 0, f"{dist}: system does not vanish at the key"
        st3 = macaulay_stats(eqs, nv, 3)
        rec = dict(n=n, exps=";".join(map(str, es)), seed=seed, distribution=dist,
                   density=round(sum(mat_density(M) for M in A) / 2, 4),
                   nvars=nv, neqs=len(eqs), mac3_rows=st3.rows, mac3_rank=st3.rank,
                   mac3_syz=st3.syz, syz_target=16 * n, mac3_dim_low=st3.dim_low,
                   mac4_dim_le3="", falls_D4="", falls_over_n="", low_cols=low_cols, saturation="")
        line = (f"n={n} seed={seed} {dist:>8} (density {rec['density']:.3f}): D=3 rank {st3.rank} "
                f"syz {st3.syz} (2*8n={16*n})")
        cens = 0
        if do_D4:
            try:
                st4 = macaulay_stats(eqs, nv, 4, budget=budget)
                falls = st4.dim_low - st3.rank
                rec.update(mac4_dim_le3=st4.dim_low, falls_D4=falls,
                           falls_over_n=round(falls / n, 4),
                           saturation=round(st4.dim_low / low_cols, 4))
                line += (f"; D=4 falls {falls} ({falls/n:.2f} n, 4n={4*n}) "
                         f"saturation {st4.dim_low/low_cols:.3f}")
            except MacaulayTimeout as exc:
                cens = 1
                line += f"; D=4 CENSORED after {budget:.0f}s (rank >= {exc.partial.rank})"
        print("   " + line)
        record("E01_theory_matrix_distribution.csv", rec, seconds=time.time() - t0, censored=cens)
        out.append(rec)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--no-D4", action="store_true")
    ap.add_argument("--budget", type=float, default=3600.0)
    ap.add_argument("--check", action="store_true",
                    help="verify that the unit-LU path reproduces E12/first_layer_ideal.py exactly")
    a = ap.parse_args()
    print(f"# params: n={a.n} seeds={a.seeds} D4={not a.no_D4} budget={a.budget}s")
    if a.check:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "E12_cross_sbox_relations"))
        from first_layer_ideal import build_first_layer
        es = [e for e in range(3, a.n - 1) if gcd(e, a.n) == 1][:2]
        rng = random.Random(1)
        A = [rand_unit_lu(a.n, rng), rand_unit_lu(a.n, rng)]
        b = [rng.randrange(F(a.n).q), rng.randrange(F(a.n).q)]
        _, mine = build_first_layer_with(a.n, es, A, b)
        nv0, _, theirs, _, _ = build_first_layer(a.n, es, 1)
        same = len(mine) == len(theirs) and all(x == y for x, y in zip(mine, theirs))
        print(f"   cross-check against E12/first_layer_ideal.build_first_layer: "
              f"{'IDENTICAL' if same else 'DIFFERENT'}")
        assert same
    for seed in [int(s) for s in a.seeds.split(",")]:
        run(a.n, seed, DISTRIBUTIONS, not a.no_D4, a.budget)
