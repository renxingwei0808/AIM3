# E01 — Enhanced-Mersenne S-box structures A–G (exhaustive small-field verification)

**Goal.** Re-run and extend the exhaustive verifications of the Focus Report §2 (closed form, image = trace-zero criterion,
kernel variable and alternating bilinear form B_e, two-valued inverse, inverse-exponent weight theorem) and keep them as
regression tests for everything downstream.

**Scripts (legacy, unchanged):** `legacy/aim3_sbox_probe.py`, `legacy/aim3_lab.py` (parts 1–3), `legacy/aim3_kernel_model.py`.
**New:** `src/aim3/structures.py` does the same checks for *any* n by random sampling (uses an F₂-linear section H_e computed by
linear algebra instead of exhaustive search), including the polarised identity B_e(u,v) = u·v·(M_e(u)+M_e(v)).

**Run.**
```
python3 legacy/aim3_lab.py            # structures A/B/C/E + cross-S-box quadratic relations (n=11..13)
python3 legacy/aim3_kernel_model.py   # structure D, dim QE(x,a)=n
python3 -m aim3.structures            # from src/: random checks at n=128/192/256 with the recommended exponents
```
**Status.** All verified. Open item: prove dim QE{(x,a)} = n for generic e (measured n=7..11); characterise the exceptional
exponents (2e ≡ ±1 mod n gives 2n) — feeds E08.
