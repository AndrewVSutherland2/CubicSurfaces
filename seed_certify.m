/* seed_certify.m -- re-certify the no-split database with the fuller cycle-type
   certificate (LinesGaloisCertificate): for each class, obtain the surface (read
   the stored one, or realize the previously-unreached ones with the current
   module) and compare its mod-q Frobenius cycle types on the 27 lines to the
   target group's.  One round-robin chunk.  Env: CHUNK_ID, NUM_CHUNKS, NOLTMP,
   OUTFILE, PROGFILE.  Output per class:
     label:source_coeffs:orbit_sizes:verdict:seen/target:primes:cubic
   verdict in {CERTIFIED, CONSISTENT, REJECTED}; "# WALL"/"# NOPOLY" for misses. */
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "cubic_surface_nosplit.m";
load "database_pipeline.m";
Q := RationalField(); Pq<t> := PolynomialRing(Q);
Pxyzw<x,y,z,w> := PolynomialRing(Q, 4);

gets := func< n | GetEnvironmentValue(n) >;
geti := func< n, dd | (s eq "" select dd else StringToInteger(s)) where s := GetEnvironmentValue(n) >;
CHUNK := geti("CHUNK_ID", 0); NCH := geti("NUM_CHUNKS", 1);
NOLTMP := gets("NOLTMP"); OUTFILE := gets("OUTFILE"); PROGFILE := gets("PROGFILE");
if NOLTMP eq "" then NOLTMP := "."; end if;

U := BuildUniversalCobleData(: Print := false);
RD := BuildResolvent27Data(U : Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
byLabel := AssociativeArray(); for s in subs do byLabel[s[1]] := s; end for;
phi27 := RD`phi27;
polyOf := AssociativeArray();
for row in [ Split(ln, ":") : ln in Split(Read("WE6fields.txt")) ] do
    dt := <StringToInteger(row[1]), StringToInteger(row[2])>;
    if not IsDefined(polyOf, dt) then polyOf[dt] := row[3]; end if;
end for;

storedSurf := AssociativeArray(); storedSrc := AssociativeArray(); storedOrb := AssociativeArray();
for ln in Split(Read("database_seed_nosplit.txt")) do
    if #ln lt 6 or ln[1] eq "#" then continue; end if;
    p := Split(ln, ":");
    storedSurf[p[1]] := p[6]; storedSrc[p[1]] := p[2]; storedOrb[p[1]] := p[3];
end for;

labels := [ k : k in Keys(storedSurf) ];
for ln in Split(Read("nosplit_unrealized.txt")) do
    if #ln ge 5 and ln[1..5] eq "51840" then
        sp := Index(ln, " "); Append(~labels, sp gt 0 select ln[1..sp-1] else ln);
    end if;
end for;
labels := Sort(labels);
mine := [ labels[i] : i in [1..#labels] | (i-1) mod NCH eq CHUNK ];

done := {};
try for ln in Split(Read(PROGFILE)) do
    pp := Split(ln, ":"); if #pp ge 1 and #pp[1] ge 1 and pp[1][1] ne "#" then Include(~done, pp[1]); end if;
end for; catch e; end try;

printf "chunk %o/%o: %o labels\n", CHUNK, NCH, #[ z0 : z0 in mine | z0 notin done ];
for label in mine do
    if label in done then continue; end if;
    s := byLabel[label]; G := s[2]; dd := s[3]; tt := s[4];
    surf := Pxyzw ! 0; src := ""; orb := "";
    haveSurf := false;
    if IsDefined(storedSurf, label) then
        surf := eval storedSurf[label]; src := storedSrc[label]; orb := storedOrb[label];
        haveSurf := true;
    elif IsDefined(polyOf, <dd,tt>) then
        coeffs := eval polyOf[<dd,tt>]; f := &+[ coeffs[k+1]*t^k : k in [0..#coeffs-1] ];
        src := polyOf[<dd,tt>];
        try
            r := CubicSurfaceNoSplittingField(G, f : Universal := U, Resolvent27 := RD, TmpDir := NOLTMP,
                     Prec := 2400, HeightBound := 6, SampleCap := 1500, MaxSurfaces := 2, Print := false);
            surf := r`surface; orb := Sprint(r`orbit_sizes); haveSurf := true;
        catch e
            PrintFile(OUTFILE, Sprintf("# WALL %o (d=%o ord=%o) : %o", label, dd, Order(G), e`Object));
        end try;
    else
        PrintFile(OUTFILE, Sprintf("# NOPOLY %o", label));
    end if;
    if haveSurf then
        cert := LinesGaloisCertificate(surf, G, phi27 : Print := false);
        PrintFile(OUTFILE, Sprintf("%o:%o:%o:%o:%o/%o:%o:%o", label, src, orb,
            cert`verdict, cert`n_seen, cert`n_target, cert`n_primes, surf));
        printf "chunk %o: %o -> %o %o/%o (%o primes)\n", CHUNK, label, cert`verdict,
            cert`n_seen, cert`n_target, cert`n_primes;
    end if;
    PrintFile(PROGFILE, label);
end for;
printf "chunk %o DONEC\n", CHUNK;
quit;
