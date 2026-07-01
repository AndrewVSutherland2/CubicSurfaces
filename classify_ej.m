/* classify_ej.m -- classify the 350 Elsenhans-Jahnel database surfaces by
   W(E6) conjugacy class of the 27-line Galois action, via mod-q Frobenius
   cycle types, and record their Delta_Cl (per-class size benchmarks).

   The database is NOT in the paper's class order (established earlier), so we
   match each surface against all 350 classes: a class is compatible when every
   observed cycle type occurs in it; among compatible classes we rank by the
   Chebotarev log-likelihood of the observed type sequence.  Classes with
   identical type sets are indistinguishable by this test and are reported as
   ties (to be separated by the exact 27-line Galois group only if needed).

   One chunk of a round-robin sweep.  Env: CHUNK_ID, NUM_CHUNKS, OUTFILE,
   PROGFILE.  Output per surface:
     idx|maxcoefdigits|dcl_digits|nprimes|ncand|best_label|best_order|ratio2|tied_labels|type_counts
   Resumable via PROGFILE. */
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "cubic_surface_nosplit.m";
load "cubic_surface_pentareduce.m";
Pxyzw<x,y,z,w> := PolynomialRing(Rationals(), 4);
load "ElsenhansJahnelDatabase.txt";

gets := func< n | GetEnvironmentValue(n) >;
geti := func< n, dd | (s eq "" select dd else StringToInteger(s)) where s := GetEnvironmentValue(n) >;
CHUNK := geti("CHUNK_ID", 0); NCH := geti("NUM_CHUNKS", 1);
OUTFILE := gets("OUTFILE"); PROGFILE := gets("PROGFILE");

/* per-class cycle-type data straight from the S27 generators */
cyc := func< g | Sort(&cat[ [ c[1] : j in [1..c[2]] ] : c in CycleStructure(g) ]) >;
labels := []; typesOf := []; weightOf := []; orderOf := [];
for ln in Split(Read("WE6subgroups.txt")) do
    p := Split(ln, ":");
    gstr := p[2];
    H := eval ("PermutationGroup< 27 | " cat gstr[2..#gstr-1] cat " >");
    w := AssociativeArray();
    for cl in ConjugacyClasses(H) do
        t := cyc(cl[3]);
        error if &+t ne 27, "cycle type does not sum to 27";
        w[t] := (IsDefined(w, t) select w[t] else 0) + cl[2];
    end for;
    Append(~labels, p[1]);
    Append(~typesOf, Keys(w));
    Append(~weightOf, w);
    Append(~orderOf, #H);
end for;
printf "chunk %o/%o: class data built (%o classes)\n", CHUNK, NCH, #labels;

done := {};
try for ln in Split(Read(PROGFILE)) do
    pp := Split(ln, "|"); if #pp ge 1 then Include(~done, StringToInteger(pp[1])); end if;
end for; catch e ; end try;

mine := [ i : i in [1..#ElsenhansJahnelDatabase] | (i-1) mod NCH eq CHUNK and i notin done ];
printf "chunk %o/%o: %o surfaces to process\n", CHUNK, NCH, #mine;

for idx in mine do
    F := Pxyzw ! ElsenhansJahnelDatabase[idx];
    inv, disc := ClebschSalmonInvariants(F);
    dcl := DeltaClOfInvariants(inv);
    maxc := Max([ #Sprint(Abs(Numerator(c))) : c in Coefficients(F) ]);

    seen := AssociativeArray();   // type -> count
    nprimes := 0; lastn := 0; stable := 0; final := [];
    for q in PrimesInInterval(11, 20000) do
        ct := FrobCycleTypeModQ(F, q, 25);
        if ct eq [] then continue; end if;
        seen[ct] := (IsDefined(seen, ct) select seen[ct] else 0) + 1;
        nprimes +:= 1;
        if nprimes ge 12 and nprimes mod 4 eq 0 then
            ok := [ i : i in [1..#labels] | Keys(seen) subset typesOf[i] ];
            if #Keys(seen) eq lastn then stable +:= 4; else stable := 0; end if;
            lastn := #Keys(seen);
            /* enough evidence: type set stable for 24 good primes */
            if stable ge 24 or nprimes ge 120 then final := ok; break; end if;
        end if;
    end for;
    if final eq [] then final := [ i : i in [1..#labels] | Keys(seen) subset typesOf[i] ]; end if;

    /* rank compatible classes by Chebotarev log-likelihood */
    lik := [];
    for i in final do
        ll := 0.0;
        for t in Keys(seen) do
            ll +:= seen[t] * Log(weightOf[i][t] / orderOf[i]);
        end for;
        Append(~lik, <ll, i>);
    end for;
    Sort(~lik, func< a,b | b[1] - a[1] gt 0 select 1 else (b[1] - a[1] lt 0 select -1 else 0) >);
    besti := lik[1][2];
    ratio2 := #lik ge 2 select Sprintf("%.2o", lik[1][1] - lik[2][1]) else "inf";
    tied := [ labels[l[2]] : l in lik | l[1] eq lik[1][1] ];
    tcounts := Sort([ <t, seen[t]> : t in Keys(seen) ]);

    PrintFile(OUTFILE, Sprintf("%o|%o|%o|%o|%o|%o|%o|%o|%o|%o",
        idx, maxc, #Sprint(Numerator(dcl)), nprimes, #final,
        labels[besti], orderOf[besti], ratio2, tied, tcounts));
    PrintFile(PROGFILE, Sprintf("%o", idx));
    printf "chunk %o: EJ[%o] -> %o (ord %o, %o cand, %o primes)\n",
        CHUNK, idx, labels[besti], orderOf[besti], #final, nprimes;
end for;
printf "chunk %o DONE\n", CHUNK;
quit;
