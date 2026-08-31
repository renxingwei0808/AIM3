"""aim3 - reference implementation and algebraic-structure toolkit for the AIM3 one-way function
(AIMer v3.0, Section 3 / Section 4.1.4). Pure Python, no dependencies.

Modules
-------
gf2n        : GF(2^n) arithmetic for arbitrary n (spec moduli for n=128/192/256, table of toy moduli)
aim3_ref    : spec-faithful AIM3 (GenerateLinear = Alg.4, NoZeroInpSB = Alg.5) with explicit encoding conventions
structures  : kernel variable a=y^r, linear section H_e, two-valued inverse, weight theorem, exponent audit helpers
"""
from .gf2n import GF2n, SPEC_MODULI, TOY_MODULI
from .aim3_ref import AIM3, Conventions, ALL_CONVENTIONS, AIM3_PARAMS
