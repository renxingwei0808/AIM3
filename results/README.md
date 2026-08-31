# results/
CSV files written by the experiment scripts. Columns common to all: `n, e/exps, ell, seed, tier, host, commit, seconds`.
* `E02_baselines.csv` — official Table 4 reproduction and baselines B1/B2 (done).
* `E03_sbox_degree_fall.csv` — per-(n,e) syzygy and degree-fall counts (append-only).
* `E08_audit_{128,192,256}.csv` — exponent audit tables (done).
Large matrices / logs from the server go to `results/server/<date>_<exp>/` (git-ignored if > 5 MB; keep the summary CSV).
