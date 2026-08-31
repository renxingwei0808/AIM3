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
