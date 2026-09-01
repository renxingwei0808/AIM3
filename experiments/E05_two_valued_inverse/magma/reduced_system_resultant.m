// E05 / track T3 (Tier C) — eliminate t1 from the two-valued-inverse reduced system.
//
// Input: a file written by `reduced_system_toy.py --export DIR`, defining
//     n, modulus, exps, seed, eps, key, pt, ct, EQ1, EQ2
// where EQk is a sequence of <exponent of t0, exponent of t1, coefficient> and the coefficient is
// an integer in the AIMer v3.0 §4.1.1 bit convention (bit i = coefficient of w^i).
//
// Magma's `load` only accepts a string literal, so this file expects the data to be loaded first.
// `experiments/E05_two_valued_inverse/run_resultants.py` writes the wrapper:
//
//     load "<DIR>/E05_n<n>_s<seed>.magma.m";
//     out := "<...>.json"; dofactor := 1; doroots := 1;
//     load "<...>/magma/reduced_system_resultant.m";
//
// What is measured (ROADMAP §3.3 / T3.2): wall-clock and memory of Resultant(Eq1, Eq2, t1), the
// degree and the number of terms of the eliminated polynomial, how many of its roots there are and
// whether the true key is one of them, and — optionally — whether it factors (the T3.3 lead).
// The termination criterion for T3 is the slope of log2(time) against n; >= 2 closes the track.

if not assigned EQ1 then
    print "reduced_system_resultant.m expects n, modulus, EQ1, EQ2, key to be loaded first";
    print "(see experiments/E05_two_valued_inverse/run_resultants.py)";
    quit;
end if;
if not assigned out then out := "e05_resultant.json"; end if;
if not assigned dofactor then dofactor := 0; end if;
if not assigned doroots then doroots := 1; end if;

P<X> := PolynomialRing(GF(2));
K<w> := ext<GF(2) | P!Intseq(modulus, 2)>;
R<t0, t1> := PolynomialRing(K, 2);

FromInt := func< c | &+[K| w^i : i in [0 .. n-1] | ((c div 2^i) mod 2) eq 1] >;
Build := func< L | &+[R| FromInt(tr[3]) * t0^tr[1] * t1^tr[2] : tr in L] >;

eq1 := Build(EQ1);
eq2 := Build(EQ2);
k0 := FromInt(key[1]);
k1 := FromInt(key[2]);

printf "### E05 n=%o exps=%o seed=%o eps=%o\n", n, exps, seed, eps;
printf "### Eq1: %o terms, deg_t0=%o deg_t1=%o\n", #Terms(eq1), Degree(eq1, t0), Degree(eq1, t1);
printf "### Eq2: %o terms, deg_t0=%o deg_t1=%o\n", #Terms(eq2), Degree(eq2, t0), Degree(eq2, t1);

v1 := Evaluate(eq1, [k0, k1]);
v2 := Evaluate(eq2, [k0, k1]);
keyok := (v1 eq 0) and (v2 eq 0);
printf "### true key is a root of both equations: %o\n", keyok;
if not keyok then
    print "### ABORT: the exported system does not vanish at the true key";
    quit;
end if;

SetMemoryLimit(0);
t := Cputime();
Res := Resultant(eq1, eq2, t1);
secs := Cputime(t);
peak := GetMaximumMemoryUsage();
zero := Res eq 0;
degR := zero select -1 else Degree(Res, t0);
nterms := zero select 0 else #Terms(Res);
printf "### Resultant: %o s, degree %o, %o terms, zero=%o, peak %o bytes\n",
       secs, degR, nterms, zero, peak;

nroots := -1; keyisroot := false; rootsecs := -1.0;
if doroots eq 1 and not zero then
    U<u> := PolynomialRing(K);
    ru := UnivariatePolynomial(Res);
    tr := Cputime();
    rts := Roots(ru);
    rootsecs := Cputime(tr);
    nroots := &+[Integers()| r[2] : r in rts];
    keyisroot := exists{ r : r in rts | r[1] eq k0 };
    printf "### Roots: %o (with multiplicity) in %o s; the true t0 is a root: %o\n",
           nroots, rootsecs, keyisroot;
end if;

nfactors := -1; factordegs := []; factorsecs := -1.0;
if dofactor eq 1 and not zero then
    U2<v> := PolynomialRing(K);
    tf := Cputime();
    fac := Factorization(UnivariatePolynomial(Res));
    factorsecs := Cputime(tf);
    nfactors := &+[Integers()| f[2] : f in fac];
    factordegs := Sort([Degree(f[1]) : f in fac]);
    printf "### Factorization: %o factors in %o s, degrees %o\n",
           nfactors, factorsecs, factordegs[1 .. Min(20, #factordegs)];
end if;

fh := Open(out, "w");
fprintf fh, "{\n";
fprintf fh, " \"n\": %o,\n", n;
fprintf fh, " \"seed\": %o,\n", seed;
fprintf fh, " \"exps\": \"%o\",\n", exps;
fprintf fh, " \"eps\": \"%o\",\n", eps;
fprintf fh, " \"eq1_terms\": %o,\n", #Terms(eq1);
fprintf fh, " \"eq2_terms\": %o,\n", #Terms(eq2);
fprintf fh, " \"eq1_deg_t1\": %o,\n", Degree(eq1, t1);
fprintf fh, " \"eq2_deg_t1\": %o,\n", Degree(eq2, t1);
fprintf fh, " \"key_is_root_of_system\": %o,\n", keyok select "true" else "false";
fprintf fh, " \"resultant_seconds\": %o,\n", secs;
fprintf fh, " \"resultant_degree\": %o,\n", degR;
fprintf fh, " \"resultant_terms\": %o,\n", nterms;
fprintf fh, " \"resultant_is_zero\": %o,\n", zero select "true" else "false";
fprintf fh, " \"peak_bytes\": %o,\n", peak;
fprintf fh, " \"n_roots\": %o,\n", nroots;
fprintf fh, " \"root_seconds\": %o,\n", rootsecs;
fprintf fh, " \"key_is_root_of_resultant\": %o,\n", keyisroot select "true" else "false";
fprintf fh, " \"n_factors\": %o,\n", nfactors;
fprintf fh, " \"factor_seconds\": %o,\n", factorsecs;
fprintf fh, " \"factor_degrees\": \"%o\",\n", factordegs[1 .. Min(50, #factordegs)];
va, vb, vc := GetVersion();
fprintf fh, " \"tool\": \"magma\",\n";
fprintf fh, " \"tool_version\": \"%o.%o-%o\"\n", va, vb, vc;
fprintf fh, "}\n";
delete fh;

quit;
