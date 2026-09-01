"""Unit tests for src/aim3/ore.py (Ore / skew polynomials over GF(2^n)).

The four properties checked here are exactly the ones track T4 (ROADMAP section 3.4) relies on:
composition = multiplication, the LCLM is a common left multiple, the GCRD's root space is the
intersection of the root spaces, and the degree identity deg lclm = deg f + deg g - deg gcrd.

All four are statements about the RING R = F_{2^n}<phi>, where phi^n != 1 and degrees add.  Only
(c) crosses over to the ACTION on F_{2^n} (where R acts through R/(phi^n - 1)); it is stated for
the root spaces of f, g and gcrd(f, g) and holds because the Bezout identity gcrd = u*f + v*g is an
identity in R and therefore also holds after applying both sides to a field element.
"""
import os, sys, random
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from aim3.gf2n import GF2n
from aim3.ore import SkewPoly, linearise_values, from_bit_matrix


def _rand(F, rng, d):
    """random SkewPoly of degree exactly d (d >= 0)."""
    c = [rng.randrange(F.q) for _ in range(d)] + [rng.randrange(1, F.q)]
    return SkewPoly(F, c)


def test_ore_mul_is_composition():
    """(a) multiplication in R is composition of the linearised maps, and it is associative."""
    rng = random.Random(20260901)
    for n in (5, 7, 8, 11):
        F = GF2n(n)
        for _ in range(20):
            f, g, h = (_rand(F, rng, rng.randrange(0, 5)) for _ in range(3))
            assert (f * g) * h == f * (g * h), f"associativity failed at n={n}"
            for _ in range(10):
                v = rng.randrange(F.q)
                assert (f * g).apply(v) == f.apply(g.apply(v)), f"composition failed at n={n}"
            # the action is F_2-linear
            u, v = rng.randrange(F.q), rng.randrange(F.q)
            assert f.apply(u ^ v) == f.apply(u) ^ f.apply(v)


def test_ore_divmod_right():
    """self = q*g + r with deg r < deg g (the precondition for the Euclidean algorithm below)."""
    rng = random.Random(1)
    for n in (5, 7, 11):
        F = GF2n(n)
        for _ in range(30):
            f = _rand(F, rng, rng.randrange(0, 7))
            g = _rand(F, rng, rng.randrange(0, 5))
            q, r = f.divmod_right(g)
            assert q * g + r == f, f"divmod_right identity failed at n={n}"
            assert r.degree() < g.degree(), f"remainder degree {r.degree()} >= {g.degree()}"


def test_ore_lclm_is_right_divisible():
    """(b) lclm(f, g) is right-divisible by both f and g."""
    rng = random.Random(2)
    for n in (5, 7, 11):
        F = GF2n(n)
        for _ in range(25):
            f = _rand(F, rng, rng.randrange(1, 5))
            g = _rand(F, rng, rng.randrange(1, 5))
            L = f.lclm(g)
            assert not L.is_zero()
            for d in (f, g):
                assert L.divmod_right(d)[1].is_zero(), f"lclm not right-divisible at n={n}"


def test_ore_gcrd_root_space_is_intersection():
    """(c) root_space(gcrd(f, g)) = root_space(f) cap root_space(g), checked EXHAUSTIVELY over
    F_{2^n} for n <= 10.

    Random f, g are coprime with overwhelming probability, which would make the test vacuous, so
    half the cases are built as f = a*h, g = b*h with a planted common right factor h.
    """
    rng = random.Random(3)

    def span(basis):
        s = {0}
        for b in basis:
            s |= {x ^ b for x in s}
        return s

    for n in (5, 7, 9, 10):
        F = GF2n(n)
        planted = 0
        for trial in range(12):
            if trial % 2:
                f, g = _rand(F, rng, rng.randrange(1, 4)), _rand(F, rng, rng.randrange(1, 4))
            else:
                h = _rand(F, rng, rng.randrange(1, 3))
                f = _rand(F, rng, rng.randrange(1, 3)) * h
                g = _rand(F, rng, rng.randrange(1, 3)) * h
            d = f.gcrd(g)
            got = span(d.root_space())
            want = {v for v in range(F.q) if f.apply(v) == 0 and g.apply(v) == 0}
            assert got == want, f"gcrd root space != intersection at n={n} (trial {trial})"
            if len(want) > 1:
                planted += 1
        assert planted > 0, f"n={n}: every case had a trivial intersection -- the test was vacuous"


def test_ore_degree_identity():
    """(d) deg lclm = deg f + deg g - deg gcrd (an identity in the ring R, not in its action)."""
    rng = random.Random(4)
    for n in (5, 7, 11, 13):
        F = GF2n(n)
        nontrivial = 0
        for trial in range(24):
            if trial % 2:
                f, g = _rand(F, rng, rng.randrange(1, 6)), _rand(F, rng, rng.randrange(1, 6))
            else:
                h = _rand(F, rng, rng.randrange(1, 4))
                f = _rand(F, rng, rng.randrange(1, 4)) * h
                g = _rand(F, rng, rng.randrange(1, 4)) * h
            L, d = f.lclm(g), f.gcrd(g)
            assert L.degree() == f.degree() + g.degree() - d.degree(), \
                f"degree identity failed at n={n}: {L.degree()} vs {f.degree()}+{g.degree()}-{d.degree()}"
            if d.degree() > 0:
                nontrivial += 1
        assert nontrivial > 0, f"n={n}: gcrd was trivial in every case"


def test_ore_from_bit_matrix_matches_legacy_apply_mat():
    """from_bit_matrix must reproduce legacy/aim3_lab.apply_mat exactly on the whole field, so that
    a unit-LU matrix of the AIM3 instance and its SkewPoly are the same F_2-linear map."""
    sys.path.insert(0, os.path.join(ROOT, "legacy"))
    from aim3_lab import apply_mat, rand_unit_lu
    rng = random.Random(5)
    for n in (5, 7, 9, 11):
        F = GF2n(n)
        M = rand_unit_lu(n, rng)
        P = from_bit_matrix(F, M)
        assert P.degree() <= n - 1
        for v in range(min(F.q, 512)):
            assert P.apply(v) == apply_mat(M, v, n), f"from_bit_matrix disagrees at n={n}, v={v}"
        # a unit-LU matrix is invertible, so the map is a bijection: trivial root space
        assert P.root_space() == []


def test_ore_linearise_values_moore_solve():
    """linearise_values must interpolate an arbitrary F_2-linear map through the Moore matrix."""
    rng = random.Random(6)
    for n in (5, 8, 11):
        F = GF2n(n)
        basis = [1 << k for k in range(n)]
        vals = [rng.randrange(F.q) for _ in range(n)]
        h = SkewPoly(F, linearise_values(F, basis, vals))
        for b, want in zip(basis, vals):
            assert h.apply(b) == want, f"Moore solve failed at n={n}"
