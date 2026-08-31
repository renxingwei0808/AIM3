# UNTESTED SKELETON (Sage 10.x). Skew polynomials over K = GF(2^n) with Frobenius phi.
# n = 8; K.<w> = GF(2^n); phi = K.frobenius_endomorphism(); R.<X> = OrePolynomialRing(K, phi)
# A random linearized polynomial A(x) = sum a_j x^{2^j}  <->  sum a_j X^j in R.
# Eq_i:  Lambda_i = (X^e * A_i) + t_i * A_i  (with t_i a K-scalar for the probe)  -> common root structure via right gcd:
#   g = Lambda_0.right_gcd(Lambda_1);  l = Lambda_0.left_lcm(Lambda_1)
# Probe: for random t_0,t_1 measure deg g, deg l; then treat t_0 symbolically (K[t_0]) and measure coefficient degrees of the
# LCLM as a function of n. Record in results/E07_lclm_degrees.csv.
print("skeleton only")
