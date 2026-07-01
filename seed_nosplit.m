/* seed_nosplit.m -- batch driver: certified cubic surfaces for the classes the
   field-based seeding skips, via cubic_surface_nosplit.m (no splitting field).

   One chunk of a round-robin parallel sweep.  Environment:
     CHUNK_ID, NUM_CHUNKS : this chunk's index and the total (round-robin split);
     NOLTMP               : a unique scratch dir for the PARI/GP call;
     OUTFILE, PROGFILE    : per-chunk output and progress files.

   Processes every W(E6) class NOT already in database_seed.txt (the small-surface
   seed) with group order in 13..2000 and minimal degree <= 10, using the
   WE6fields.txt polynomial for its (d,t).  Appends, per class,

     label:source_coeffs:orbit_sizes:frobenius_orders:consistent:cubic_form

   ("consistent" = the mod-q Frobenius orders are all element orders of the target
   group, i.e. the surface is not accidentally generic).  Resumable via PROGFILE.

   Launch (from repo root) with a shell loop, e.g.
     for k in $(seq 0 15); do NOLTMP=.../t$k OUTFILE=.../o$k PROGFILE=.../p$k \
       CHUNK_ID=$k NUM_CHUNKS=16 magma -b seed_nosplit.m & done
*/
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "cubic_surface_nosplit.m";
load "database_pipeline.m";
Q := RationalField(); Pq<t> := PolynomialRing(Q);

gets := func< n | GetEnvironmentValue(n) >;
geti := func< n, d | (s eq "" select d else StringToInteger(s)) where s := GetEnvironmentValue(n) >;
CHUNK  := geti("CHUNK_ID", 0);  NCH := geti("NUM_CHUNKS", 1);
NOLTMP := gets("NOLTMP");  OUTFILE := gets("OUTFILE");  PROGFILE := gets("PROGFILE");
if NOLTMP eq "" then NOLTMP := "."; end if;

U    := BuildUniversalCobleData(: Print := false);
RD   := BuildResolvent27Data(U : Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
polyOf := AssociativeArray();
for r in [ Split(x, ":") : x in Split(Read("WE6fields.txt")) ] do
    dt := <StringToInteger(r[1]), StringToInteger(r[2])>;
    if not IsDefined(polyOf, dt) then polyOf[dt] := r[3]; end if;
end for;

done := {};
for fn in ["database_seed.txt", PROGFILE] do
    try for l in Split(Read(fn)) do
        p := Split(l, ":");
        if #p ge 1 and #p[1] ge 1 and p[1][1] ne "#" then Include(~done, p[1]); end if;
    end for; catch e ; end try;
end for;

cand := [];
for si in [1..#subs] do s := subs[si];
    if s[1] in done or not IsDefined(polyOf, <s[3],s[4]>)
       or Order(s[2]) gt 2000 or s[3] gt 10 then continue; end if;
    Append(~cand, <Order(s[2]), si>);
end for;
Sort(~cand, func< a,b | a[1]-b[1] >);
myidx := [ cand[i][2] : i in [1..#cand] | (i-1) mod NCH eq CHUNK ];
printf "chunk %o/%o: %o classes to process\n", CHUNK, NCH, #myidx;

for si in myidx do
    s := subs[si]; label := s[1]; G := s[2]; dd := s[3]; tt := s[4];
    coeffs := eval polyOf[<dd,tt>];
    f := &+[ coeffs[k+1]*t^k : k in [0..#coeffs-1] ];
    ok := false; surf := 0; osz := []; orders := {}; consistent := false; prec := 1200;
    for attempt in [1..3] do
        try
            r := CubicSurfaceNoSplittingField(G, f : Universal := U, Resolvent27 := RD,
                     TmpDir := NOLTMP, Prec := prec, MaxSurfaces := 2, Print := false);
            surf := r`surface; osz := r`orbit_sizes;
            cts := LinesFrobeniusCycleTypesModQ(surf : MaxQ := 8, QStop := 500, Print := false);
            orders := { LCM(c[2]) : c in cts };
            consistent := orders subset { c[1] : c in ConjugacyClasses(r`galois_group) };
            ok := true; break;
        catch e
            prec *:= 2;
        end try;
    end for;
    if ok then
        PrintFile(OUTFILE, Sprintf("%o:%o:%o:%o:%o:%o", label, polyOf[<dd,tt>], osz, orders, consistent, surf));
    else
        PrintFile(OUTFILE, Sprintf("# FAIL %o (order %o)", label, Order(G)));
    end if;
    PrintFile(PROGFILE, label);
    printf "chunk %o: %o ord %o -> ok=%o digits=%o orders=%o consistent=%o\n",
        CHUNK, label, Order(G), ok, ok select #Sprint(Max([Abs(c):c in Coefficients(surf)])) else 0, orders, consistent;
end for;
printf "chunk %o DONE\n", CHUNK;
quit;
