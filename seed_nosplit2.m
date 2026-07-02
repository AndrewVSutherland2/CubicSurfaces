/* seed_nosplit2.m -- batch driver for the classes the first no-split sweep
   skipped: group order up to 26000 and minimal degree up to 27 (the original
   seed_nosplit.m capped at order 2000 / degree 10; WE6fields.txt already has
   polynomials for degrees 12, 16, 18, 24, 27).  Unlike the first sweep this
   emits (K5, lambda)-reduced models directly (cubic_surface_pentareduce.m),
   screening the candidate moduli points by the trial-division-resistant part
   of E = I40, and certifies with the full cycle-type certificate.

   One chunk of a round-robin sweep.  Env: CHUNK_ID, NUM_CHUNKS, NOLTMP,
   OUTFILE, PROGFILE.  Output per class:
     label:source_coeffs:orbit_sizes:verdict:digits:cubic
   ("# WALL label" = realization failed, "# NOPOLY label" = no field on file.)
   Resumable via PROGFILE. */
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

U  := BuildUniversalCobleData(: Print := false);
RD := BuildResolvent27Data(U : Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
polyOf := AssociativeArray();
for r in [ Split(x, ":") : x in Split(Read("WE6fields.txt")) ] do
    dt := <StringToInteger(r[1]), StringToInteger(r[2])>;
    if not IsDefined(polyOf, dt) then polyOf[dt] := r[3]; end if;
end for;
/* per-label overrides (gap_fields.txt): fields bootstrapped from the EJ
   surfaces / LMFDB / compositum constructions for the classes WE6fields.txt
   does not cover.  Any f with Gal(f) isomorphic to the class group works. */
overrides := AssociativeArray();
try for l in Split(Read("gap_fields.txt")) do
    p := Split(l, ":");
    if #p eq 2 and #p[1] ge 1 and p[1][1] ne "#" then overrides[p[1]] := p[2]; end if;
end for; catch e ; end try;

done := {};
for fn in ["database_seed.txt", "database_seed_nosplit.txt", PROGFILE] do
    try for l in Split(Read(fn)) do
        p := Split(l, ":");
        if #p ge 1 and #p[1] ge 1 and p[1][1] ne "#" then Include(~done, p[1]); end if;
    end for; catch e ; end try;
end for;

cand := [];
for si in [1..#subs] do s := subs[si];
    if s[1] in done or s[3] eq 1 or s[3] gt 27 or Order(s[2]) gt 52000
       or (not IsDefined(polyOf, <s[3],s[4]>) and not IsDefined(overrides, s[1]))
       then continue; end if;
    Append(~cand, <Order(s[2]), si>);
end for;
Sort(~cand, func< a,b | a[1]-b[1] >);
myidx := [ cand[i][2] : i in [1..#cand] | (i-1) mod NCH eq CHUNK ];
printf "chunk %o/%o: %o classes to process\n", CHUNK, NCH, #myidx;

EScore := function(inv)
    n := Abs(Numerator(inv[5]) * Denominator(inv[5]));
    if n le 1 then return 0; end if;
    _, comps := TrialDivision(n, 10^6);
    return #comps eq 0 select 0 else &+[ #Sprint(c) : c in comps ];
end function;

for si in myidx do
    s := subs[si]; label := s[1]; G := s[2]; dd := s[3]; tt := s[4];
    if label in done then continue; end if;
    polystr := IsDefined(overrides, label) select overrides[label] else polyOf[<dd,tt>];
    coeffs := eval polystr;
    f := &+[ coeffs[k+1]*t^k : k in [0..#coeffs-1] ];

    ok := false; r := 0; prec := 1200;
    for attempt in [1..3] do
        try
            r := CubicSurfaceNoSplittingField(G, f : Universal := U, Resolvent27 := RD,
                     TmpDir := NOLTMP, Prec := prec, MaxSurfaces := 4,
                     SplitPrimeBound := 3000000, InvariantsOnly := true, Print := false);
            ok := true; break;
        catch e
            prec *:= 2;
        end try;
    end for;
    if not ok then
        PrintFile(OUTFILE, Sprintf("# WALL %o (d=%o ord=%o)", label, dd, Order(G)));
        PrintFile(PROGFILE, label);
        printf "chunk %o: %o ord %o WALL\n", CHUNK, label, Order(G);
        continue;
    end if;

    cands := r`all_candidates;
    scored := Sort([ <EScore(cands[i][3]), #Sprint(cands[i][1]), i> : i in [1..#cands] ]);
    best := Pxyzw!0; bestdig := 0; tried := 0; nleft := -1;
    for sc in scored do
        if tried ge 3 then break; end if;
        tried +:= 1;
        cnd := cands[sc[3]];
        rok, res, msg := PentahedralReduce(cnd[3] : TmpDir := NOLTMP, GPProg := GPPROG);
        if rok and (bestdig eq 0 or res`digits lt bestdig) then
            best := res`surface; bestdig := res`digits; nleft := #res`leftover;
        end if;
        if bestdig gt 0 and bestdig le 24 then break; end if;
    end for;
    if bestdig eq 0 then          // reduction failed everywhere: build the raw model
        best := Pxyzw ! penta_pol_2_gl(ClebschInvariantsToPentahedralPolynomial(cands[scored[1][3]][3]));
        bestdig := Max([ #Sprint(Abs(Numerator(c))) : c in Coefficients(best) ]);
    end if;

    cert := LinesGaloisCertificate(best, G, RD`phi27 : Print := false);
    PrintFile(OUTFILE, Sprintf("%o:%o:%o:%o:%o:%o",
        label, polystr, Sprint(r`orbit_sizes), cert`verdict, bestdig, best));
    PrintFile(PROGFILE, label);
    printf "chunk %o: %o d=%o ord %o -> %o digits %o\n",
        CHUNK, label, dd, Order(G), bestdig, cert`verdict;
end for;
printf "chunk %o DONE\n", CHUNK;
quit;
