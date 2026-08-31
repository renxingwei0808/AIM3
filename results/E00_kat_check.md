# E00 — KAT pinning of the AIM3 encoding conventions

Source: official AIMer reference-implementation repository (`samsungsds-opensource/AIMer`, commit with message *update AIMer V3.0 reference implementation*), file `AIMer.zip`, directory `KAT/`.  The bundled `Reference_Implementation/aim3.c` carries the AIM3 exponents (3,7,11)/(11,23,47)/(3,7,11,19), i.e. it is the v3.0 primitive.

Key layout read from `sign.c::crypto_sign_keypair`: `pk = iv || ct`, `sk = pt || iv || ct`; field elements are little-endian bytes (`field_common.c::gf_from_bytes`).

The `f` and `s` variants of each level share the same KAT seed and therefore the same AIM3 records (checked for all six directories); only the `f` files are listed.  The first three records of each set are committed to `experiments/E00_reference_and_facts/kat_vectors.json`, so `python3 experiments/E00_reference_and_facts/kat_pin.py --vectors` reproduces this table offline.

| parameter set | n | KAT record 0 reproduced by | conventions that match |
|---|---|---|---|
| AIM3-I | 128 | ct = `d46c65bf72e9968d032713a4e90eb1ab` | `lsb/A/col/LU`, `lsb/AT/row/LU` |
| AIM3-III | 192 | ct = `a198e107d6204c6e9d40f7bd54b8c52d4638136a1ea4d1dd` | `lsb/A/col/LU`, `lsb/AT/row/LU` |
| AIM3-V | 256 | ct = `86d5d269a4738b7915dd52a348b1047304cd210991046092d7a50394c4af1a0a` | `lsb/A/col/LU`, `lsb/AT/row/LU` |

Convention axes (`src/aim3/aim3_ref.py::Conventions`): `tape_bit_order` (lsb/msb) / `matrix_apply` (A/AT) / `tape_matrix_order` (col/row) / `lu_order` (LU/UL).

* AIM3-I: `pt=dd1a6bdbe4106d0caa9476b0a035997c` `iv=91282214654cb55e7c2cacd53919604d` `ct=d46c65bf72e9968d032713a4e90eb1ab`
* AIM3-III: `pt=0348b1cc251ad82fdd1a6bdbe4106d0caa9476b0a035997c` `iv=8626ed79d451140800e03b59b956f8210e556067407d13dc` `ct=a198e107d6204c6e9d40f7bd54b8c52d4638136a1ea4d1dd`
* AIM3-V: `pt=2d7f73369973cd2d0348b1cc251ad82fdd1a6bdbe4106d0caa9476b0a035997c` `iv=8626ed79d451140800e03b59b956f8210e556067407d13dc90fa9e8b872bfb8f` `ct=86d5d269a4738b7915dd52a348b1047304cd210991046092d7a50394c4af1a0a`
