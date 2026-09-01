# E07 — Ore / skew-polynomial elimination (Magic-Pot generalisation, track T4)

**Idea.** Over F_{2^n}, each AIM3 relation is *Frobenius-linear* in one block and *scalar-affine* in the other:
first-layer  Eq_i: (φ^{e_i}∘A_i + t_i·A_i)(pt) = c_i(t_i)   — an element Λ_i(t_i) of the Ore ring R = F_{2^n}[t_i]⟨φ⟩ (φa = a²φ)
acting on pt;  last layer: Λ₂(pt)(t₀,t₁) = c₂(pt).  Magic Pot's F_{2^n}[pt]-linear algebra becomes module computations over R.

**Core questions.**
1. Common-root elimination of pt between Λ₀ and Λ₁ via LCLM / GCRD (Ore's Euclidean algorithm) and the *skew resultant*: what
   is the degree, in the conjugates t₀^{2^k}, t₁^{2^k}, of the resulting scalar condition? (Conjecture: ≈ n in each — the Dickson
   determinant of x^{2^e}+t·x is a norm-like polynomial; then Frobenius consistency must be imposed.)
2. Row-degree growth of weak-Popov / approximant bases over R with coefficients in F_{2^n}[t]: polynomial or exponential in n?
3. Is there a scalar s ∈ F_{2^n} (candidates: t₂ = pt+ct, a₂, or a normal-basis coordinate) in which *all* equations are
   polynomial of moderate degree? (Obstruction: every choice meets Frobenius closure somewhere; a negative theorem is publishable.)

**Tools.** `src/aim3/ore.py` (pure-Python Ore ring, Tier A; unit-tested in `tests/test_ore.py`), Sage
`OrePolynomialRing(K, K.frobenius_endomorphism(), 'phi')` (`sage/skew_probe.sage`), Magma `TwistedPolynomialRing` on the server.
**Kill criterion.** If for n = 8..20 the LCLM/approximant row degrees grow like 2^{cn} with c ≥ 1/2, close the attack line and
write the impossibility argument instead (still a T4 deliverable).

---

## Tier-A run (batch 3): `skew_elimination_probe.py`

```
python3 skew_elimination_probe.py --selftest --probe --degree --control --n 6 7 8 9 10 11 12 13 14 --trials 32
python3 skew_elimination_probe.py --joint  --n 6 7 8 9 10 11
python3 skew_elimination_probe.py --export instance_n11.json --n 11 --no-csv
```
Output: `results/E07_theory_skew_degrees.csv` (one row per (n, measurement)).  Runtime ≈ 2 min in total on a laptop.

**Exponents.** `COMPLEXITY_ACCOUNTING.md` rule 5 restricts toy instances to exponents of the *generic* class
(`degenerate_class(n,e) == []`).  For n ≤ 10 and n = 12 the only generic exponent is n−1, so both S-boxes get the same e
there; that is not an AIM3 parameter shape (the spec uses distinct e_i) but it keeps the instance out of the degenerate
classes, which is what a degree measurement needs.  n = 11, 13, 14 use two distinct generic exponents.

### Set-up

With `x_i = A_i(pt) + b_i`, `t_i = y_i = S[e_i](x_i)`, the S-box relation `x_i y_i = x_i^{2^{e_i}} + 1` is, in R,

    Λ_i(t_i) = (φ^{e_i} + t_i)·A_i ,      Eq_i:  Λ_i(t_i)(pt) = c_i(t_i) = 1 + b_i^{2^{e_i}} + t_i·b_i .

`Λ_i(t_i)(pt) = c_i(t_i)` is **affine** in pt, not linear.  Eliminating pt with the LCLM cofactors
`L = u·Λ₀ = v·Λ₁` gives the scalar condition

    E(t₀, t₁) := u(c₀(t₀)) + v(c₁(t₁)) = 0 ,

which is the skew resultant of T4.1.  Alongside it the script measures the *linear-algebra* predicate

    C(t₀, t₁) := [ the affine F_2 system {Λ₀(v) = c₀, Λ₁(v) = c₁} in v has a solution ] ,

which is the statement "there exists a pt consistent with both equations".

### Correction to the batch-3 hypothesis

The brief conjectured "at the true (t₀,t₁) the GCRD has degree ≥ 1, because pt and pt+a lie in the intersection of the two
root spaces; at random (t₀,t₁) the GCRD is 0".  **That is false, and not by accident.**  `ker(φ^{e}+t) = {0, t^{r}}`
(r = (2^e−1)^{-1} mod 2^n−1) is *one*-dimensional, so `ker Λ_i = A_i^{-1}({0, a_i})` is a 1-dimensional F_2 space; two
1-dimensional subspaces of an n-dimensional space meet only in 0 unless they coincide, which needs `a₁ = A₁A₀^{-1}(a₀)` —
probability ≈ 2^{-n}, and *nothing about the true key forces it*, because the true key satisfies the **affine** equations,
not the homogeneous ones.  Measured (`--probe`, n = 6..14, 32 random (t₀,t₁) each): `deg gcrd = 0` and
`dim(ker Λ₀ ∩ ker Λ₁) = 0` at the true key in every single case, exactly as at random.  The statistic that *does* separate
the true key is the affine one: `C = 1` at the true (t₀,t₁) in every n, and 0 for 287 of the 288 random pairs tried
(9 values of n x 32 trials; the one exception is at n = 9, consistent with the ≈2·2^{-n} chance of an accidental intersection).

### Core measurement: the degree of the eliminated condition

`t₁` is fixed at its true value and `t₀` runs over **all** of F_{2^n}, so both degrees below are computed by exhaustion
and are **exact values, not interpolation lower bounds**:

* Boolean (ANF) degree, by a Möbius transform over the 2^n values of t₀ — exact;
* univariate degree over F_{2^n}, from `c_j = Σ_{a≠0} E(a)·a^{-j}` (every map F_{2^n}→F_{2^n} is a polynomial of degree
  ≤ 2^n−1) — exact.

Both instruments are checked against closed-form degrees by `--selftest` before any measurement is trusted.

| n | e₀,e₁ | deg Λ_i | deg lclm (in R) | Boolean deg E | univariate deg E | Boolean deg C | \|C\| | \|Z(E)\| |
|---|-------|---------|-----------------|---------------|------------------|---------------|-----|--------|
| 6 | 5,5 | 10,10 | 20 | 5 | 62 (of 63) | 5 | 2 | 2 |
| 7 | 6,6 | 12,12 | 24 | 7 | 127 | 6 | 2 | 2 |
| 8 | 7,7 | 14,14 | 28 | 8 | 255 | 7 | 2 | 2 |
| 9 | 8,8 | 16,16 | 32 | 9 | 511 | 8 | 2 | 2 |
| 10 | 9,9 | 18,18 | 36 | 10 | 1023 | 9 | 2 | 2 |
| 11 | 3,4 | 13,14 | 27 | 11 | 2047 | 10 | 2 | 2 |
| 12 | 11,11 | 22,22 | 44 | 12 | 4095 | 11 | 2 | 2 |
| 13 | 3,5 | 15,17 | 32 | 13 | 8191 | 12 | 2 | 2 |
| 14 | 3,5 | 16,18 | 34 | 14 | 16383 | 13 | 2 | 2 |

`Z(E) = C` in every run (mutual inclusion checked both ways), i.e. the condition derived in R — where φ^n ≠ 1 — loses
nothing against the action on F_{2^n}: no separate Frobenius-consistency step was needed here.  Exactly one t₀ per instance
makes `gcrd(Λ₀,Λ₁)` nontrivial in R; it is not the true one and does not affect the degrees.

`--control` repeats the measurement with a **wrong** t₁: then |C| ∈ {0, 2} (0 for 3 of the 9 n, 2 for the rest) and the true
t₀ is never in C.

### Conclusion

**The eliminated condition has maximal degree.**  D(n) = 2^n − 1 in the univariate sense and n (the maximum for n bits) in
the Boolean sense — measured exactly at n = 6..14, so the answer to the T4.3 trichotomy is **exponential**, not O(n) or O(n²).
For the LCLM elimination this is not an artefact of the toy sizes but a structural fact, and the structure proves it for
every n:

> **Lemma [PROVED].** Fix t₁ ≠ 0 with gcd(e₁,n) = 1 and A₁ invertible.  Eq₁ is an affine equation whose homogeneous part
> `v ↦ Λ₁(v)` has a 1-dimensional kernel, so its solution set in pt is empty or a coset {pt, pt+w} of size 2 with
> `A₁(w) = a₁ = t₁^{r₁}`.  Each of the (at most two) candidate pt determines `x₀ = A₀(pt)+b₀` and hence at most one
> `t₀ = S[e₀](x₀)`.  Therefore `|C(·, t₁)| ≤ 2`.  The indicator of a set {α, β} ⊆ F_2^n with α ≠ β has ANF degree exactly
> n−1 (the coefficient of the weight-(n−1) monomial omitting bit i is [α_i ≠ β_i]).  Hence `deg_ANF C = n−1` for every n.

The experiment confirms the lemma's hypotheses hold in practice: `|C| = 2` at every n = 6..14 and the second element is
**exactly** the predicted `S[e₀](A₀(pt+w)+b₀)` (column "C = predicted 2-set: 1" in the CSV detail field), and with a wrong
t₁ the count is 0 or 2 as the lemma allows.

So the skew resultant of the first layer is informationally equivalent to *guessing t₀ and testing it*: it is the indicator of
a 2-element set, not a low-degree algebraic relation.  T4.1's conjecture "total degree ≈ n" is numerically right and
operationally empty — n is the **largest** possible Boolean degree on n bits.  Jointly in (t₀,t₁) the picture is the same:
`|C|` is 2^n − 2 or 2^n − 3 out of 2^{2n} (n = 6..11) with Boolean degree 2n or 2n−1 of a possible 2n.

**Status.** [PROVED] the ≤ 2 lemma and the n−1 ANF degree of C; [VERIFIED n≤14] the exact degrees of E, the |C| = 2 count and
Z(E) = C; [OPEN] whether some *other* elimination order (last layer first, or an approximant basis over
F_{2^n}[t]⟨φ⟩ that keeps t symbolic throughout) avoids the collapse — the measurement above closes the LCLM route only, and
ROADMAP §3.4's kill criterion is stated for approximant-basis row degrees, which are not measured here.
**No extrapolation is made from these curves**: every number in the table is an exhaustive count at that n.

### Cross-validation input for the server (Sage / Magma)

`--export FILE.json` writes one instance (`instance_n11.json` is committed, n = 11, e = (3,4)):

| key | meaning |
|-----|---------|
| `n`, `modulus` | field F_2[X]/(modulus); bit i of an element is the coefficient of X^i (AIMer v3.0 §4.1.1) |
| `e` | the two S-box exponents (e₀, e₁) |
| `A_bitmatrix` | A_i as an n×n bit matrix, `y_i = XOR_j A[i][j]·x_j` (legacy `apply_mat` convention) |
| `A_skew_coeffs` | the same maps as Ore polynomials: `A_i(v) = Σ_j a_j v^{2^j}` |
| `b`, `x`, `t`, `c` | the constants, the true S-box inputs/outputs, and `c_i(t_i)` |
| `pt` | the planted key |

A Sage cross-check must reproduce: `Λ_i = (phi^e_i + t_i)*A_i` with `Λ_i(pt) = c_i`; `Λ₀.right_gcd(Λ₁)` of degree 0;
`deg Λ₀.left_lcm(Λ₁) = deg Λ₀ + deg Λ₁`; and the same two solvable t₀ values.  `sage/skew_probe.sage` loads the JSON.
