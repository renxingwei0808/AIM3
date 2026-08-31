import os, sys, random
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from aim3 import AIM3, GF2n
from aim3.structures import SboxStructures, weight_data, degenerate_class

def test_field_and_sbox_identities():
    F = GF2n(11)
    for x in range(1, 200):
        for e in (2, 3, 4):
            assert F.sbox(x, e) == F.mul(F.frob(x ^ 1, e), F.inv(x))

def test_weight_theorem_recommended():
    for n, exps, expect in [(128, (3, 7, 11), (43, 55, 35)), (192, (11, 23, 47), (35, 167, 143)), (256, (3, 7, 11, 19), (171, 183, 163, 27))]:
        for e, w in zip(exps, expect):
            d = weight_data(n, e)
            assert d["wt_r"] == w == d["m"] and d["wt_minus_r"] == n - w and d["supp_is_AP"]
            assert degenerate_class(n, e) == []

def test_toy_aim3_selftest():
    assert AIM3(11, (2, 3, 4)).selftest(random.Random(0), 3)

def test_structures_full_size():
    assert SboxStructures(128, 3).check_random(random.Random(1), 3)

def test_boolmac_matches_legacy():
    """src/aim3/boolmac.py (one-pass, streamed) must agree with legacy/aim3_macaulay.macaulay."""
    sys.path.insert(0, os.path.join(ROOT, "legacy"))
    from aim3_macaulay import build_aim3, macaulay
    from aim3_sbox_syz import sbox_system
    from aim3.boolmac import macaulay_stats
    nv, eqsP1, eqsAll, _ = build_aim3(7, [2, 3, 4], seed=1)
    for eqs, V in ((eqsAll, nv), (eqsP1, nv)):
        for D in (2, 3):
            rk, low, nrow, ncol = macaulay(eqs, V, D)
            f = macaulay_stats(eqs, V, D)
            assert (rk, low, nrow, ncol) == (f.rank, f.dim_low, f.rows, f.cols)
    _, P1, P2, P3 = sbox_system(9, 4)
    for D in (3, 4):
        rk, low, nrow, ncol = macaulay(P1 + P2 + P3, 18, D)
        f = macaulay_stats(P1 + P2 + P3, 18, D)
        assert (rk, low, nrow, ncol) == (f.rank, f.dim_low, f.rows, f.cols)

def test_kat_conventions_pinned():
    """AIM3 with the DEFAULT Conventions must reproduce the official AIMer v3.0 KAT records
    (experiments/E00_reference_and_facts/kat_vectors.json, extracted by kat_pin.py)."""
    import json
    from aim3 import AIM3, AIM3_PARAMS
    path = os.path.join(ROOT, "experiments", "E00_reference_and_facts", "kat_vectors.json")
    data = json.load(open(path))
    by_n = {p["n"]: p["exps"] for p in AIM3_PARAMS.values()}
    for name, s in data["sets"].items():
        A = AIM3(s["n"], by_n[s["n"]])
        for r in s["records"]:
            assert A.evaluate(int(r["pt"], 16), bytes.fromhex(r["iv"])) == int(r["ct"], 16), \
                f"{name}: KAT mismatch — the pinned Conventions no longer reproduce the official vectors"


def test_export_matches_boolmac():
    """src/aim3/boolmac_export.py must describe the SAME matrix boolmac.py streams, and the column
    reversal the Sage backend performs must reproduce boolmac's rank / dim_low exactly.

    The Sage/M4RI engine itself is validated against this one on real instances by
    server/jobs/verify_backend.py; this test locks the file format and the column order, which is
    the part that can silently go wrong without any error being raised.
    """
    import tempfile
    from math import gcd
    sys.path.insert(0, os.path.join(ROOT, "legacy"))
    from aim3_macaulay import build_aim3
    from aim3.boolmac import macaulay_stats, poly_rows, monomial_columns
    from aim3.boolmac_export import export_macaulay, read_macaulay, poly_row_indices

    def rank_like_sage(path):
        """What rank_backend.sage does: load column j at position ncols-1-j, echelonise left to
        right, and split the pivots at the degree-D boundary."""
        meta, rows = read_macaulay(path)
        nc, ts = meta["ncols"], meta["top_start"]
        piv = {}
        for a in rows:
            r = 0
            for j in a:
                r |= 1 << (nc - 1 - j)
            while r:
                b = (r & -r).bit_length() - 1          # leftmost column = lowest reversed bit
                p = piv.get(b)
                if p is None:
                    piv[b] = r
                    break
                r ^= p
        rank = len(piv)
        rank_top = sum(1 for b in piv if b < nc - ts)
        return meta, rank, rank - rank_top

    for n, ell, D in [(5, 1, 3), (5, 1, 4), (5, 2, 3), (7, 1, 3)]:
        exps = [e for e in range(2, n) if gcd(e, n) == 1][:ell + 1]
        nv, eqsP1, eqsAll, _ = build_aim3(n, exps)
        for eqs in (eqsP1, eqsAll):
            st = macaulay_stats(eqs, nv, D)
            col, _, _ = monomial_columns(nv, D)
            streamed = list(poly_rows(eqs, nv, D, col=col))
            exported = [sum(1 << c for c in idx) for idx in poly_row_indices(eqs, nv, D, col=col)]
            assert streamed == exported, f"row mismatch at n={n} ell={ell} D={D}"
            with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as fh:
                path = fh.name
            try:
                export_macaulay(eqs, nv, D, path)
                meta, rank, dim_low = rank_like_sage(path)
                assert (meta["nrows"], meta["ncols"]) == (st.rows, st.cols)
                assert rank == st.rank, f"rank {rank} != {st.rank} at n={n} ell={ell} D={D}"
                assert dim_low == st.dim_low, f"dim_low {dim_low} != {st.dim_low}"
            finally:
                for p in (path, path + ".meta.json"):
                    if os.path.exists(p):
                        os.remove(p)


def test_legacy_field_table_extended():
    """src/aim3/legacy_fields.py must add irreducible moduli for n = 17..24 to the legacy table
    without changing any entry it shares with it (otherwise Tier-A and server instances at n <= 16
    would describe different fields)."""
    sys.path.insert(0, os.path.join(ROOT, "legacy"))
    import aim3.legacy_fields as lf
    import aim3_lab
    from aim3.gf2n import TOY_MODULI
    for n in (17, 18, 19, 20):
        assert n in aim3_lab.CONWAY, f"n={n} still missing from the legacy field table"
        assert lf.is_irreducible(aim3_lab.CONWAY[n], n)
    for n, mod in TOY_MODULI.items():
        if n in aim3_lab.CONWAY:
            assert aim3_lab.CONWAY[n] == mod, f"modulus disagreement at n={n}"
