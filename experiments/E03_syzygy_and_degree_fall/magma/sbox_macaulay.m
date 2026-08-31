// UNTESTED SKELETON (write on the server with Magma; verify against results/E03_sbox_degree_fall.csv for n<=13 first)
// Degree-D Macaulay rank of the single enhanced-Mersenne S-box system {P1,P2,P3} in the Boolean ring, for larger n.
// Usage: magma n:=20 e:=3 D:=3 sbox_macaulay.m
n := StringToInteger(n); e := StringToInteger(e); D := StringToInteger(D);
K<w> := GF(2, n);
B := BooleanPolynomialRing(2*n, "grevlex");       // vars 1..n = x bits, n+1..2n = y bits
basis := [w^i : i in [0..n-1]];                    // polynomial basis  (matches spec bit i <-> X^i when w is the class of X)
X := &+[ B.(i+1)*basis[i+1] : i in [0..n-1] ];     // NOTE: this is an element of B tensor K; Magma needs the coordinate trick below
// --- coordinate expansion helper: field-valued affine expression -> n Boolean polynomials ---
// Represent a field-valued Boolean-affine expression as a list of K-coefficients indexed by Boolean monomials.
// To keep this skeleton short, we build P1,P2,P3 by evaluating trace forms:  coordinate_j(F) = Tr(lambda_j * F)
// where {lambda_j} is the dual basis of {basis}.  Implement:  dual := DualBasis(basis) (Magma: use Trace pairing).
// P1(x,y) = x*y + x^(2^e) + 1 ; P2 = x*P1 ; P3 = y*P1  (as in spec §6.3.1)
// TODO(server): implement coordinate expansion with a K-linear-combination-of-monomials data structure,
// then build the Macaulay matrix M_D (rows = monomial multiples of degree <= D-2, columns = monomials of degree <= D)
// and print Rank(M_D), NumberOfRows(M_D), the syzygy count, and dim(rowspace ∩ degree<=D-1) via echelon on the
// degree-D columns first. Cross-check: n=11,e=3 must give rank 671, rows 759, syz 88 (8n), D4 falls 22 (2n).
print "skeleton only - see README";
