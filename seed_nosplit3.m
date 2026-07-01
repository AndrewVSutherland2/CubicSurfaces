/* seed_nosplit3.m -- third pass over the classes the batch could not certify.
   Reads a newline-list of labels from FAILFILE; for each, walks a (HeightBound,
   Prec) ladder and captures the failing exception so height-limited cases get
   rescued and genuinely size-walled cases get classified.  Same output format
   as seed_nosplit.m, plus a trailing "# WALL label : <exception>" on total failure. */
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "cubic_surface_nosplit.m";
load "database_pipeline.m";
Q := RationalField(); Pq<t> := PolynomialRing(Q);
gets := func< n | GetEnvironmentValue(n) >;
geti := func< n, d | (s eq "" select d else StringToInteger(s)) where s := GetEnvironmentValue(n) >;
CHUNK  := geti("CHUNK_ID", 0);  NCH := geti("NUM_CHUNKS", 1);
NOLTMP := gets("NOLTMP");  OUTFILE := gets("OUTFILE");  PROGFILE := gets("PROGFILE"); FAILFILE := gets("FAILFILE");
if NOLTMP eq "" then NOLTMP := "."; end if;

U    := BuildUniversalCobleData(: Print := false);
RD   := BuildResolvent27Data(U : Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
byLabel := AssociativeArray();
for s in subs do byLabel[s[1]] := s; end for;
polyOf := AssociativeArray();
for r in [ Split(x, ":") : x in Split(Read("WE6fields.txt")) ] do
    dt := <StringToInteger(r[1]), StringToInteger(r[2])>;
    if not IsDefined(polyOf, dt) then polyOf[dt] := r[3]; end if;
end for;

done := {};
try for l in Split(Read(PROGFILE)) do
    p := Split(l, ":"); if #p ge 1 and #p[1] ge 1 and p[1][1] ne "#" then Include(~done, p[1]); end if;
end for; catch e ; end try;

labels := [ l : l in Split(Read(FAILFILE)) | #l gt 0 ];
mine := [ labels[i] : i in [1..#labels] | (i-1) mod NCH eq CHUNK and labels[i] notin done ];
printf "chunk %o/%o: %o labels\n", CHUNK, NCH, #mine;

ladder := [ 2400, 4800 ];
for label in mine do
    s := byLabel[label]; G := s[2]; dd := s[3]; tt := s[4];
    coeffs := eval polyOf[<dd,tt>];
    f := &+[ coeffs[k+1]*t^k : k in [0..#coeffs-1] ];
    ok := false; surf := 0; osz := []; orders := {}; consistent := false; errtxt := "";
    for rung in ladder do
        try
            r := CubicSurfaceNoSplittingField(G, f : Universal := U, Resolvent27 := RD,
                     TmpDir := NOLTMP, Prec := rung, HeightBound := 6, SampleCap := 1500,
                     MaxSurfaces := 2, Print := false);
            surf := r`surface; osz := r`orbit_sizes;
            cts := LinesFrobeniusCycleTypesModQ(surf : MaxQ := 8, QStop := 500, Print := false);
            orders := { LCM(c[2]) : c in cts };
            consistent := orders subset { c[1] : c in ConjugacyClasses(r`galois_group) };
            ok := true; break;
        catch e
            errtxt := Sprint(e`Object);
        end try;
    end for;
    if ok then
        PrintFile(OUTFILE, Sprintf("%o:%o:%o:%o:%o:%o", label, polyOf[<dd,tt>], osz, orders, consistent, surf));
        printf "chunk %o: %o (d=%o ord=%o) -> OK digits=%o orders=%o consistent=%o\n",
            CHUNK, label, dd, Order(G), #Sprint(Max([Abs(c):c in Coefficients(surf)])), orders, consistent;
    else
        e1 := #errtxt gt 90 select errtxt[1..90] else errtxt;
        PrintFile(OUTFILE, Sprintf("# WALL %o (d=%o ord=%o) : %o", label, dd, Order(G), e1));
        printf "chunk %o: %o (d=%o ord=%o) -> WALL : %o\n", CHUNK, label, dd, Order(G), e1;
    end if;
    PrintFile(PROGFILE, label);
end for;
printf "chunk %o DONE3\n", CHUNK;
quit;
