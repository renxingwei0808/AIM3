# E10 — Block-order / elimination-order Gröbner experiments (FreeLunch-style, track T2)

Test whether a variable order aligned with the twisted-bilinear block structure lowers the F4 step degrees compared to plain
grevlex: orders {pt > T}, {T > pt}, {W > U} (kernel model), and product orders; also "natural equations already close to a
Gröbner basis" checks (reduce S_quad w.r.t. itself under each order).  Magma skeleton in `magma/`. Inputs: the same toy
instances as E04 (export from Python). Metrics: max step degree, #matrix rows at the top step, time, memory. Tier C.

## Server log — 2026-08-31 (server batch 1, branch `claude/server-batch-1`)

`magma/aim3_blockorder_f4.m` is superseded by
`experiments/E04_full_system_hilbert/magma/aim3_f4_degrees.m`, which takes the order as a parameter
(`grevlex`, `lex`, `elim:K`, `elimlast:K`, `gen:grevlex`) and is driven by `run_f4_degrees.py` with
the symbolic labels `blockPT` / `blockT` (pt block first / t block first) and `blockU` / `blockW`
(the kernel model's U and W blocks).  Rows: `results/E10_server_block_orders.csv`.

**Two facts about Magma 2.24 that this experiment must be read through.**

1. `BooleanPolynomialRing(n, order)` takes an order *name* only — an `"elim"` order raises
   *"Bad argument types"*.  Block orders therefore have to run on
   `PolynomialRing(GF(2), n, "elim", S)` with the field equations `x²+x` added by hand, which is a
   different (slower) engine.  The only fair baseline for those rows is `gen:grevlex` — grevlex on
   the *same* ring with the same field equations — not the Boolean-ring `grevlex` of E04.  Every row
   records which ring was used.
2. The declared elimination order **is** honoured in the returned Gröbner basis (checked directly:
   the same toy system gives a 30-element basis under grevlex and a 15-element one under
   `elim`), but `SetVerbose("Faugere", 2)` reports
   `Order: Graded Reverse Lexicographical (bit vector word)` in *both* cases: Magma's GF(2) F4 runs
   in grevlex and converts afterwards.  **So a step-degree log under a block order is a grevlex log,
   and `d_solv` cannot be compared across orders with this tool.**  What is comparable is the total
   time and the shape of the resulting basis, and those are what the report uses.

A third fact limits the whole experiment on these instances: the toy systems have a **single**
solution over F₂, and the reduced Gröbner basis of such an ideal is the same in every monomial
order.  A block-order comparison can only bite where the solution set is larger or where the
intermediate F4 profile is visible — neither holds here.  Details:
`results/REPORT_server_batch1.md` §5.
