# Complexity accounting rules (binding for every number that goes into the paper)

1. **Units.** State every cost in one of: bit operations, F_{2^n}-multiplications, or AIM3 evaluations — and give the conversion
   used. **The conversion is now fixed** (batch 2, `experiments/E11_hybrid_estimates/aim3_eval_cost.py`,
   `results/E11_theory_eval_cost.csv`); the 2^{14} placeholder is retired.

   Spec Table 5 (AIMer v3.0 English specification, p.42) gives the operation counts of one **forward**
   AIM3 evaluation, taken from the addition chains of its S-boxes (x^{-1} by Itoh–Tsujii, evaluated in
   parallel with the Mersenne part x^{2^e−1}):

   | scheme | F_{2^n} mul (count, depth) | F_{2^n} square (count, depth) | Grover total (spec) |
   |---|---|---|---|
   | AIM3-I | 47, 24 | 399, 254 | 2^{165.6} |
   | AIM3-III | 59, 26 | 651, 382 | 2^{232.5} |
   | AIM3-V | 73, 28 | 1056, 510 | 2^{298.0} |

   Independent check: the same chain rule (Itoh–Tsujii needs ⌊log₂(n−1)⌋+wt(n−1)−1 multiplications,
   the Mersenne part ⌊log₂ e⌋+wt(e)−1, plus one multiplication per S-box) gives 50 / 63 / 77
   multiplications with no sharing, against the spec's 47 / 59 / 73 — the 3–4 difference is the sharing
   between the Itoh–Tsujii and Mersenne x^{2^k−1} towers. The spec numbers are the ones used.

   **Bit-operation model** (constants stated, not hidden): squaring and modular reduction cost 4n bit
   ops each (bit-interleaving is free and every spec modulus is a trinomial/pentanomial); a field
   multiplication costs n²+(n−1)²+4n (schoolbook) or 6·n^{log₂3}+4n (Karatsuba); an n×n binary
   matrix–vector product costs n²+n(n−1); a field addition costs n. In a **single-key** search the iv is
   fixed, so the 2ℓ matrices are generated from SHAKE once and the tape is *not* a per-candidate cost;
   the per-candidate cost is the S-box chain plus the 2ℓ matrix–vector products.

   | scheme | bits / evaluation (schoolbook) | (Karatsuba) | log₂ |
   |---|---|---|---|
   | AIM3-I | 1 887 791 | 976 414 | **20.85 / 19.90** |
   | AIM3-III | 5 167 931 | 2 312 744 | **22.30 / 21.14** |
   | AIM3-V | 11 473 993 | 4 816 758 | **23.45 / 22.20** |

   So one AIM3 evaluation costs 2^{19.9}–2^{23.5} bit operations, **6 to 9 bits more than the retired
   2^{14} placeholder**. Every earlier statement that converted between AIM3 evaluations and bit
   operations with 2^{14} understated the bit-operation side by that much.
2. **Baselines.** Compare against B1 (official 262.6/381.0/651.5), B2 (Boolean-series 218.4/312.5/529.1) and B3 (exhaustive
   2^n evaluations). A result "beats B1 but not B3" is an *improved algebraic attack*, not a security break; say so explicitly.
   **B3 in bit operations** (2^n × the table above): AIM3-I 2^{148.9} / 2^{147.9}, AIM3-III 2^{214.3} / 2^{213.1},
   AIM3-V 2^{279.5} / 2^{278.2} (schoolbook / Karatsuba). B1 and B2 are already bit-operation counts (they
   count field operations at ω=2 over F₂), so the three are directly comparable only in that unit — state
   which convention a row uses. Note B1 < B3 for all three parameter sets, i.e. the official algebraic
   estimate is already below exhaustive search in this unit; that is a property of the baselines, not a
   result of this project.
3. **ω.** Report with ω = 2 (as the designers) **and** ω = 2.81; never mix.
4. **Extrapolation discipline.** Extrapolate only when the same model shows a stable law across ≥ 4 values of n and ≥ 2 exponents;
   censored (timed-out / out-of-memory) runs are lower bounds only; never fit an exponential to a single small-n point; separate
   field ops / bit ops / memory; separate single-key and multi-target models.
5. **Toy exponents must be generic.** Every toy instance used for a *quantitative* claim must take its
   exponents from the generic class, i.e. `aim3.structures.degenerate_class(n, e) == []`, unless the
   experiment's purpose is precisely to study a degenerate class — in which case say so in the row.
   The default exponent list `[e for e in range(2, n) if gcd(e,n)==1][:l+1]` starts at **e = 2**, which
   is the most degenerate class there is (7n falls instead of 2n, and the only class with cubics of
   the ideal unreachable at D=4), so a script that uses the default is not automatically compliant.
   **Consequence already on the record:** the `d_solv` table of `results/REPORT_server_batch1.md` §5
   was produced with the default list and therefore contains e = 2 in every instance; it is a valid
   *measurement of those instances* but **must not be used for extrapolation** or compared against
   B1's `d_reg`. The server branch re-runs it with generic exponents.

6. **Spurious solutions.** Every model must report its number of solutions on toy instances (2-to-1 S-boxes create extra
   solutions at small n); an attack cost includes the filtering.
7. **Guessing.** Hybrid costs are min over k of 2^k · (solver cost with k bits fixed) *plus* the cost of detecting wrong guesses.
8. **Provenance.** Every table row links to a results CSV with commit hash, tier, machine, and wall-clock time.
