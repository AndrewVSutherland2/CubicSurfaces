/* reemit_nosplit.m -- re-emit the no-split database rows with (K5, lambda)-
   reduced models (cubic_surface_pentareduce.m).

   For each class in database_seed_nosplit.txt: regenerate candidate moduli
   points with CubicSurfaceNoSplittingField, screen the candidates by how much
   of E = I40 survives trial division (the unfactored part of E is exactly
   what blocks the ideal-level cube-stripping), pentahedrally reduce the most
   factorable few, certify the winner's 27-line action with the full cycle-type
   certificate, and emit a replacement row.  One chunk of a round-robin sweep.

   Env: CHUNK_ID, NUM_CHUNKS, NOLTMP, OUTFILE, PROGFILE.  Output per class:
     label:source_coeffs:orbit_sizes:verdict:digits:cubic
   ("# KEEP label ..." comment rows mark classes where no reduction succeeded;
   the merge step then retains the old row.)  Resumable via PROGFILE. */
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "cubic_surface_nosplit.m";
load "database_pipeline.m";
load "cubic_surface_pentareduce.m";
Q := RationalField(); Pq<t> := PolynomialRing(Q);
Pxyzw<x,y,z,w> := PolynomialRing(Q, 4);

gets := func< n | GetEnvironmentValue(n) >;
geti := func< n, dd | (s eq "" select dd else StringToInteger(s)) where s := GetEnvironmentValue(n) >;
CHUNK := geti("CHUNK_ID", 0); NCH := geti("NUM_CHUNKS", 1);
NOLTMP := gets("NOLTMP"); OUTFILE := gets("OUTFILE"); PROGFILE := gets("PROGFILE");
if NOLTMP eq "" then NOLTMP := "."; end if;
GPPROG := GetCurrentDirectory() cat "/pentareduce.gp";
MAXSURF := geti("MAXSURF", 4);        // candidate moduli points per class
NREDUCE := geti("NREDUCE", 3);        // how many of them to pentahedrally reduce
SKIPLE  := geti("SKIPLE", 0);         // skip classes whose stored model already has <= this many digits

U  := BuildUniversalCobleData(: Print := false);
RD := BuildResolvent27Data(U : Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
byLabel := AssociativeArray(); for s in subs do byLabel[s[1]] := s; end for;

storedSrc := AssociativeArray(); storedOrb := AssociativeArray();
for ln in Split(Read("database_seed_nosplit.txt")) do
    if #ln eq 0 or ln[1] eq "#" then continue; end if;
    p := Split(ln, ":");
    if SKIPLE gt 0 then
        dg := Max([ #x : x in Split(p[#p], "+-*^ wxyz") | x ne "" ]);
        if dg le SKIPLE then continue; end if;
    end if;
    storedSrc[p[1]] := p[2]; storedOrb[p[1]] := p[3];
end for;

labels := Sort([ k : k in Keys(storedSrc) ]);
mine := [ labels[i] : i in [1..#labels] | (i-1) mod NCH eq CHUNK ];
done := {};
try for ln in Split(Read(PROGFILE)) do
    pp := Split(ln, ":"); if #pp ge 1 and #pp[1] ge 1 and pp[1][1] ne "#" then Include(~done, pp[1]); end if;
end for; catch e ; end try;
printf "chunk %o/%o: %o labels to process\n", CHUNK, NCH, #[ z : z in mine | z notin done ];

/* digits of the trial-division-resistant part of E: the screening score */
EScore := function(inv)
    n := Abs(Numerator(inv[5]) * Denominator(inv[5]));
    if n le 1 then return 0; end if;
    _, comps := TrialDivision(n, 10^6);
    return #comps eq 0 select 0 else &+[ #Sprint(c) : c in comps ];
end function;

for label in mine do
    if label in done then continue; end if;
    s := byLabel[label]; G := s[2];
    coeffs := eval storedSrc[label];
    f := &+[ coeffs[k+1]*t^k : k in [0..#coeffs-1] ];

    ok := false; r := 0; prec := 1200;
    for attempt in [1..3] do
        try
            r := CubicSurfaceNoSplittingField(G, f : Universal := U, Resolvent27 := RD,
                     TmpDir := NOLTMP, Prec := prec, MaxSurfaces := MAXSURF,
                     InvariantsOnly := true, Print := false);
            ok := true; break;
        catch e
            prec *:= 2;
        end try;
    end for;
    if not ok then
        PrintFile(OUTFILE, Sprintf("# KEEP %o (regeneration failed)", label));
        PrintFile(PROGFILE, label);
        printf "chunk %o: %o REGEN-FAIL\n", CHUNK, label;
        continue;
    end if;

    cands := r`all_candidates;
    scored := Sort([ <EScore(cands[i][3]), #Sprint(cands[i][1]), i> : i in [1..#cands] ]);
    best := Pxyzw!0; bestdig := 0; tried := 0; nleft := -1;
    for sc in scored do
        if tried ge NREDUCE then break; end if;
        tried +:= 1;
        cand := cands[sc[3]];
        rok, res, msg := PentahedralReduce(cand[3] : TmpDir := NOLTMP, GPProg := GPPROG);
        if rok and (bestdig eq 0 or res`digits lt bestdig) then
            best := res`surface; bestdig := res`digits; nleft := #res`leftover;
        end if;
        if bestdig gt 0 and bestdig le 24 then break; end if;
    end for;

    if bestdig eq 0 then
        PrintFile(OUTFILE, Sprintf("# KEEP %o (no reduction succeeded)", label));
        PrintFile(PROGFILE, label);
        printf "chunk %o: %o REDUCE-FAIL\n", CHUNK, label;
        continue;
    end if;

    cert := LinesGaloisCertificate(best, G, RD`phi27 : Print := false);
    PrintFile(OUTFILE, Sprintf("%o:%o:%o:%o:%o:%o",
        label, storedSrc[label], Sprint(r`orbit_sizes), cert`verdict, bestdig, best));
    PrintFile(PROGFILE, label);
    printf "chunk %o: %o ord %o -> %o digits (%o candidates, %o leftover) %o\n",
        CHUNK, label, Order(G), bestdig, #cands, nleft, cert`verdict;
end for;
printf "chunk %o DONE\n", CHUNK;
quit;
