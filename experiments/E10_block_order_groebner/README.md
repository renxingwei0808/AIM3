# E10 — Block-order / elimination-order Gröbner experiments (FreeLunch-style, track T2)

Test whether a variable order aligned with the twisted-bilinear block structure lowers the F4 step degrees compared to plain
grevlex: orders {pt > T}, {T > pt}, {W > U} (kernel model), and product orders; also "natural equations already close to a
Gröbner basis" checks (reduce S_quad w.r.t. itself under each order).  Magma skeleton in `magma/`. Inputs: the same toy
instances as E04 (export from Python). Metrics: max step degree, #matrix rows at the top step, time, memory. Tier C.
