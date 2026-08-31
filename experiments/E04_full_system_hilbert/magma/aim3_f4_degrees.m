// F4 step-degree log for the Boolean systems exported by src/aim3/polysys_export.py.
//
// Magma's `load` only accepts a string LITERAL, so the system file cannot be named on the command
// line.  `experiments/E04_full_system_hilbert/run_f4_degrees.py` therefore writes a three-line
// wrapper and runs that:
//
//     load "<SYS>.magma.m";           // defines nvars, sysname, MONS (monomial index lists, 1-based)
//     ord := "grevlex"; out := "<OUT>.json";
//     load "<...>/aim3_f4_degrees.m";
//
// Run it by hand the same way; this file expects nvars/sysname/MONS/ord/out to be assigned already.
// The RING IS BUILT HERE, not in the generated file, so the same exported system can be run under
// several monomial orders -- that is what E10's block-order comparison needs:
//
//   ord := "grevlex"            the E04 reference order
//   ord := "lex"
//   ord := "elim:K"             block order with the FIRST K variables in the leading block
//   ord := "elimlast:K"         block order with the LAST  K variables in the leading block
//   ord := "gen:grevlex"        grevlex on the same general ring the block orders need
//
// IMPORTANT (Magma 2.24): `BooleanPolynomialRing` accepts a plain order NAME only -- an "elim"
// order raises "Bad argument types".  Block orders therefore run on
// `PolynomialRing(GF(2), nvars, "elim", S)` with the field equations x_i^2 + x_i added by hand.
// That ring is slower than the dedicated Boolean one, so a block-order timing must be compared
// with `gen:grevlex` (same ring, same field equations), never with `grevlex`.  The JSON records
// which ring was used in "ring".
//
// For the S_quad layout (variables 1..n are pt, the rest are t_0..t_{l-1}) that gives
// "pt block first" = elim:n and "t block first" = elimlast:(nvars-n); for the E06 kernel model
// (U = 1..3n, W = 3n+1..6n) it gives U>W = elim:3n and W>U = elimlast:3n.  See
// experiments/E10_block_order_groebner/magma/aim3_blockorder_f4.m.
//
// stdout carries the SetVerbose("Faugere",2) log; `src/aim3/polysys_export.parse_magma_log`
// extracts the per-step degrees from it and `d_solv` is their maximum.  This file itself writes
// only the summary JSON, so a killed run still leaves a parsable log.
//
// Cross-validated against msolve at n=7 by experiments/E04_full_system_hilbert/run_f4_degrees.py.

if not assigned MONS then
    print "aim3_f4_degrees.m expects nvars, sysname, MONS to be loaded first";
    print "(see experiments/E04_full_system_hilbert/run_f4_degrees.py)";
    quit;
end if;
if not assigned out then out := "f4.json"; end if;
if not assigned ord then ord := "grevlex"; end if;
if assigned nthreads then
    SetNthreads(StringToInteger(nthreads));
end if;

parts := Split(ord, ":");
if #parts eq 1 then
    ring := "BooleanPolynomialRing";
    B := BooleanPolynomialRing(nvars, parts[1]);
    x := [B.i : i in [1..nvars]];
    F := [ &+[B| &*[B| x[i] : i in m] : m in poly] : poly in MONS ];
else
    ring := "PolynomialRing(GF(2)) + field equations";
    if parts[1] eq "gen" then
        B := PolynomialRing(GF(2), nvars, parts[2]);
    elif parts[1] eq "elim" then
        B := PolynomialRing(GF(2), nvars, "elim", [1 .. StringToInteger(parts[2])]);
    elif parts[1] eq "elimlast" then
        k := StringToInteger(parts[2]);
        B := PolynomialRing(GF(2), nvars, "elim", [nvars - k + 1 .. nvars]);
    else
        printf "unknown order %o\n", ord;
        quit;
    end if;
    x := [B.i : i in [1..nvars]];
    F := [ &+[B| &*[B| x[i] : i in m] : m in poly] : poly in MONS ]
         cat [ x[i]^2 + x[i] : i in [1 .. nvars] ];
end if;

printf "### system %o : %o polynomials in %o variables, order %o, ring %o\n",
       sysname, #F, nvars, ord, ring;
va, vb, vc := GetVersion();
printf "### magma V%o.%o-%o, nthreads %o\n", va, vb, vc, GetNthreads();

SetVerbose("Faugere", 2);
t0 := Cputime();
m0 := GetMemoryUsage();
G := GroebnerBasis(F);
secs := Cputime(t0);
SetVerbose("Faugere", 0);

peak := GetMaximumMemoryUsage();
trivial := (#G eq 1) and (G[1] eq B!1);
degs := [ #Monomials(g) eq 0 select 0 else Degree(g) : g in G ];
nsol := -1;
if assigned variety and not trivial then
    nsol := #Variety(Ideal(G));
end if;

printf "### GB size %o, max degree %o, %o s, peak memory %o bytes, trivial %o, solutions %o\n",
       #G, Max(degs), secs, peak, trivial, nsol;

fh := Open(out, "w");
fprintf fh, "{\n";
fprintf fh, " \"system\": \"%o\",\n", sysname;
fprintf fh, " \"order\": \"%o\",\n", ord;
fprintf fh, " \"ring\": \"%o\",\n", ring;
fprintf fh, " \"nvars\": %o,\n", nvars;
fprintf fh, " \"neqs\": %o,\n", #F;
fprintf fh, " \"gb_size\": %o,\n", #G;
fprintf fh, " \"gb_max_degree\": %o,\n", Max(degs);
fprintf fh, " \"seconds\": %o,\n", secs;
fprintf fh, " \"peak_bytes\": %o,\n", peak;
fprintf fh, " \"trivial_gb\": %o,\n", trivial select "true" else "false";
fprintf fh, " \"n_solutions\": %o,\n", nsol;
fprintf fh, " \"tool\": \"magma\",\n";
fprintf fh, " \"tool_version\": \"%o.%o-%o\"\n", va, vb, vc;
fprintf fh, "}\n";
delete fh;

quit;
