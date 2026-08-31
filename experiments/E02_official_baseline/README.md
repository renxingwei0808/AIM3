# E02 — Exact reproduction of the official Table 4 and the comparison baselines

**Result (done).** `table4_reproduce.py` reproduces AIMer v3.0 Table 4 to 0.1 bit:

| scheme | official (k, d_reg, log₂T) | formula that reproduces it |
|---|---|---|
| AIM3-I   | 57, 19, 262.6 | d_reg from **(1−z²)^m/(1−z)^{N−k}** (generic, non-Boolean series); cost 2^k·C(N−k+d,d)² |
| AIM3-III | 102, 25, 381.0 | same |
| AIM3-V   | 150, 44, 651.5 | same |

The spec *states* the Boolean series (1+z)^N/((1−z)(1+z²)^m) (its eq. (2)) but the numbers come from the generic series.
With the Boolean series and the same cost model the optimum is **218.4 / 312.5 / 529.1** (k = 5/2/9).

**Baselines to be used everywhere in the paper.**
* **B1** = official Table 4 numbers (262.6 / 381.0 / 651.5).
* **B2** = Boolean-series semi-regular numbers (218.4 / 312.5 / 529.1).
* **B3** = exhaustive search: 2^n AIM3 evaluations (bit-op cost to be fixed in E00; placeholder 2^{n+14}).
* Any claimed attack states its cost in the same unit and against all three.

**Why this matters for a "new techniques" paper.** It is not a correction of the design (we do not claim the designers erred);
it fixes *which number* a new attack must beat and shows the semi-regular model is sensitive to the series choice — motivating
the structured Hilbert function of T1.
