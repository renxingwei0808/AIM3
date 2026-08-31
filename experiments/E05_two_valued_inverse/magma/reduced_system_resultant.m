// UNTESTED SKELETON: two-valued-inverse reduced system for toy AIM3 (l=2) over K=GF(2^n), branch bits (eps0,eps1,eps2).
//   Unknowns t0,t1 in K.  x0 = G_{e0,eps0}(t0),  pt = A0^{-1}(x0 + b0),
//   Eq1:  A1(pt) + b1 = G_{e1,eps1}(t1)
//   Eq2:  b2 + B0(t0) + B1(t1) = G_{e2,eps2}(pt + ct)
// Measure: degrees, number of terms, timing of Resultant(Eq1,Eq2,t1) for n = 8..24, and the number of spurious roots.
// G_{e,eps}(y) := y^r * H(y^{-r}) + eps*y^r + 1, where H is an F2-linear section of s -> s^(2^e)+s (compute by linear algebra).
// Compare against the AIM2-style reduced system (degree 2^e - 1) to quantify "non-permutation patch = X bits".
print "skeleton only - see README";
