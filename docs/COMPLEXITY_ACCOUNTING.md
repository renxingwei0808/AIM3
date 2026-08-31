# Complexity accounting rules (binding for every number that goes into the paper)

1. **Units.** State every cost in one of: bit operations, F_{2^n}-multiplications, or AIM3 evaluations — and give the conversion
   used (fix cost(AIM3-I evaluation) in E00 by counting the addition chains of Table 5, spec p.42; placeholder 2^{14} bit ops).
2. **Baselines.** Compare against B1 (official 262.6/381.0/651.5), B2 (Boolean-series 218.4/312.5/529.1) and B3 (exhaustive
   2^n evaluations). A result "beats B1 but not B3" is an *improved algebraic attack*, not a security break; say so explicitly.
3. **ω.** Report with ω = 2 (as the designers) **and** ω = 2.81; never mix.
4. **Extrapolation discipline.** Extrapolate only when the same model shows a stable law across ≥ 4 values of n and ≥ 2 exponents;
   censored (timed-out / out-of-memory) runs are lower bounds only; never fit an exponential to a single small-n point; separate
   field ops / bit ops / memory; separate single-key and multi-target models.
5. **Spurious solutions.** Every model must report its number of solutions on toy instances (2-to-1 S-boxes create extra
   solutions at small n); an attack cost includes the filtering.
6. **Guessing.** Hybrid costs are min over k of 2^k · (solver cost with k bits fixed) *plus* the cost of detecting wrong guesses.
7. **Provenance.** Every table row links to a results CSV with commit hash, tier, machine, and wall-clock time.
