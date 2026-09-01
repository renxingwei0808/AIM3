# Exported systems for the F4 step-degree runs (task 6)

This directory holds the **descriptors** of the 81 systems that `run_f4_degrees.py` was pointed at:
one `<tag>.json` per system with its shape (variables, equations, term count, maximum degree), its
provenance (experiment, n, ell, seed, exponents, builder, host, commit) and the planted solution,
plus `index.json` listing them all.

The **inputs themselves** — `<tag>.ms` (msolve) and `<tag>.magma.m` (Magma) — are 14 MB and stay on
the server under `/home/luck/rxw/AIM3_remote/out/E04_systems/`, per batch rule 0.4: they are
regenerable inputs, not results. Regenerate them exactly with

    python3 experiments/E04_full_system_hilbert/export_systems.py --out DIR \
        --ns 7,9,11,13,15 --ells 1,2 --seeds 1,2
    for n in 7 9 11 13; do
      python3 experiments/E06_kernel_bilinear_model/kernel_bilinear_toy.py $n 2 --export DIR
    done

The generators are deterministic in (n, ell, seed), so a regenerated file is byte-identical and the
`solution` field in each descriptor still applies.
