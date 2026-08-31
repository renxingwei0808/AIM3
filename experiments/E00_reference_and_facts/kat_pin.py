#!/usr/bin/env python3
"""E00 — pin the AIM3 encoding conventions against an official AIMer KAT response file.

The AIMer KAT files are NIST `PQCsignKAT_<skbytes>.rsp` files.  For AIMer,
    pk = iv || ct,      sk = pt || iv || ct,
each of AIM3_NUM_BYTES_FIELD = n/8 bytes, and a field element is read from bytes little-endian
(`gf_from_bytes` in the reference implementation: 64-bit words, each little-endian, in increasing
order -- which for a whole n-bit element is exactly `int.from_bytes(buf, "little")`).

The script takes one or more .rsp files, extracts (pt, iv, ct) from the first `count = 0` record, and
evaluates AIM3 under every convention in `aim3.ALL_CONVENTIONS` until ct is reproduced.

The first three records of each official KAT are also committed to `kat_vectors.json` in this
directory, so the pinning is reproducible without network access:
  python3 kat_pin.py --vectors           (re-runs the search against the committed vectors)

Usage:
  python3 kat_pin.py PQCsignKAT_48.rsp [more.rsp ...] [--records K] [--md]
  python3 kat_pin.py --vectors [--md]
    --records K  also cross-check the winning convention on the first K records (default 3)
    --md         (re)write results/E00_kat_check.md

Conventions searched (see src/aim3/aim3_ref.py):
  tape_bit_order    lsb | msb
  matrix_apply      A | AT
  tape_matrix_order col | row
  lu_order          LU | UL
"""
import sys, os, time, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _common import RESULTS
from aim3 import AIM3, Conventions, ALL_CONVENTIONS, AIM3_PARAMS

PARAM_BY_N = {p["n"]: (name, p["exps"]) for name, p in AIM3_PARAMS.items()}


def parse_rsp(path, limit=None):
    """-> (n, [(pt, iv_bytes, ct), ...]) reading `sk = ...` lines."""
    recs = []
    n = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("sk ="):
                continue
            sk = bytes.fromhex(line.split("=", 1)[1].strip())
            nb = len(sk) // 3
            n = nb * 8
            pt = int.from_bytes(sk[:nb], "little")
            iv = sk[nb:2 * nb]
            ct = int.from_bytes(sk[2 * nb:3 * nb], "little")
            recs.append((pt, iv, ct))
            if limit and len(recs) >= limit:
                break
    if n is None:
        raise SystemExit(f"{path}: no `sk =` line found — is this a PQCsignKAT .rsp file?")
    return n, recs


VECTORS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kat_vectors.json")


def load_vectors():
    """-> [(label, n, [(pt, iv_bytes, ct), ...]), ...] from the committed kat_vectors.json."""
    data = json.load(open(VECTORS))
    out = []
    for name, s in data["sets"].items():
        n = s["n"]
        out.append((name, n, [(int(r["pt"], 16), bytes.fromhex(r["iv"]), int(r["ct"], 16))
                              for r in s["records"]]))
    return out


def pin_records(label, n, recs, records=3):
    recs = recs[:max(records, 1)]
    if n not in PARAM_BY_N:
        raise SystemExit(f"{label}: field size n={n} is not an AIM3 parameter set {sorted(PARAM_BY_N)}")
    name, exps = PARAM_BY_N[n]
    pt, iv, ct = recs[0]
    print(f"{label}: {name} n={n} exps={exps}, {len(recs)} record(s)")
    print(f"   pt={pt:0{n//4}x}\n   iv={iv.hex()}\n   ct={ct:0{n//4}x}")
    winners = []
    for conv in ALL_CONVENTIONS:
        t0 = time.time()
        got = AIM3(n, exps, conv=conv).evaluate(pt, iv)
        ok = got == ct
        print(f"   {conv.tape_bit_order:>3}/{conv.matrix_apply:>2}/{conv.tape_matrix_order}/"
              f"{conv.lu_order}: {'MATCH' if ok else 'no'}   ({time.time()-t0:.1f}s)")
        if ok:
            winners.append(conv)
    kept = []
    for conv in winners:
        A = AIM3(n, exps, conv=conv)
        extra = all(A.evaluate(p, i) == c for p, i, c in recs[1:])
        print(f"   -> {conv} reproduces record 0"
              + (f" and records 1..{len(recs)-1}" if extra and len(recs) > 1 else
                 (" but FAILS a later record" if not extra else "")))
        if extra:
            kept.append(conv)
    return name, n, exps, recs[0], kept


def pin(path, records=3):
    n, recs = parse_rsp(path, limit=max(records, 1))
    return pin_records(os.path.basename(path), n, recs, records)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("rsp", nargs="*")
    ap.add_argument("--vectors", action="store_true",
                    help="use the committed kat_vectors.json instead of .rsp files")
    ap.add_argument("--records", type=int, default=3)
    ap.add_argument("--md", action="store_true")
    a = ap.parse_args()
    if a.vectors or not a.rsp:
        print(f"# params: committed kat_vectors.json records={a.records}")
        out = [pin_records(lbl, n, recs, a.records) for lbl, n, recs in load_vectors()]
    else:
        print(f"# params: files={a.rsp} records={a.records}")
        out = [pin(p, a.records) for p in a.rsp]
    print()
    lines = ["# E00 — KAT pinning of the AIM3 encoding conventions", "",
             "Source: official AIMer reference-implementation repository "
             "(`samsungsds-opensource/AIMer`, commit with message *update AIMer V3.0 reference "
             "implementation*), file `AIMer.zip`, directory `KAT/`.  The bundled "
             "`Reference_Implementation/aim3.c` carries the AIM3 exponents (3,7,11)/(11,23,47)/"
             "(3,7,11,19), i.e. it is the v3.0 primitive.", "",
             "Key layout read from `sign.c::crypto_sign_keypair`: `pk = iv || ct`, "
             "`sk = pt || iv || ct`; field elements are little-endian bytes "
             "(`field_common.c::gf_from_bytes`).", "",
             "The `f` and `s` variants of each level share the same KAT seed and therefore the same "
             "AIM3 records (checked for all six directories); only the `f` files are listed.  The first "
             "three records of each set are committed to "
             "`experiments/E00_reference_and_facts/kat_vectors.json`, so "
             "`python3 experiments/E00_reference_and_facts/kat_pin.py --vectors` reproduces this table "
             "offline.", "",
             "| parameter set | n | KAT record 0 reproduced by | conventions that match |",
             "|---|---|---|---|"]
    for name, n, exps, (pt, iv, ct), winners in out:
        w = ", ".join(f"`{c.tape_bit_order}/{c.matrix_apply}/{c.tape_matrix_order}/{c.lu_order}`"
                      for c in winners) or "**none of the 16**"
        lines.append(f"| {name} | {n} | ct = `{ct:0{n//4}x}` | {w} |")
    lines += ["", "Convention axes (`src/aim3/aim3_ref.py::Conventions`): "
              "`tape_bit_order` (lsb/msb) / `matrix_apply` (A/AT) / `tape_matrix_order` (col/row) / "
              "`lu_order` (LU/UL).", ""]
    for name, n, exps, (pt, iv, ct), winners in out:
        lines.append(f"* {name}: `pt={pt:0{n//4}x}` `iv={iv.hex()}` `ct={ct:0{n//4}x}`")
    print("\n".join(lines))
    if a.md:
        os.makedirs(RESULTS, exist_ok=True)
        with open(os.path.join(RESULTS, "E00_kat_check.md"), "w") as f:
            f.write("\n".join(lines) + "\n")
        print("\nwritten results/E00_kat_check.md")
