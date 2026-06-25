/* gen_hard_examples.m -- collect cubic-surface models that do NOT minimize
   within a 120s timeout, for S3 (3T2) number fields.  Sweeps several S3 classes
   to span Delta_Cl, then keeps the SMALLEST-Delta_Cl non-minimizers (the most
   compact genuinely-hard models).  Writes unminimized_examples.txt. */
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "database_pipeline.m";
Q := RationalField(); Z := IntegerRing();
Pf<a> := PolynomialRing(Q);
R4<x,y,z,w> := PolynomialRing(Z, 4);
dir := "/tmp/claude-1000/-home-claude-CubicSurfaces/dfcfed9a-a863-46fb-80be-f21032adb5c7/scratchpad/race";

U := BuildUniversalCobleData(: Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
bm := U`bm;

keepchars := {@ "0","1","2","3","4","5","6","7","8","9","-","," @};
strip := function(s)
    t := ""; for i in [1..#s] do c := s[i]; if c in keepchars then t cat:= c; end if; end for; return t;
end function;
nflabel := AssociativeArray();
for ln in Split(Read("lmfdb_fields.csv")) do
    ci := Index(ln, ",");
    if ci gt 0 then lab := Substring(ln,1,ci-1);
        if lab ne "label" then nflabel[strip(Substring(ln,ci+1,#ln))] := lab; end if; end if;
end for;

s3 := [ Split(r,":") : r in Split(Read("WE6fields.txt")) | Split(r,":")[1] eq "3" and Split(r,":")[2] eq "2" ];
S3idx := [ i : i in [1..#subs] | subs[i][3] eq 3 and subs[i][4] eq 2 ];
NC := Min(4, #S3idx);    /* sweep this many S3 classes per field */
printf "%o S3 fields, sweeping %o of %o S3 classes\n", #s3, NC, #S3idx;

cand := [];   /* <gl, fieldpoly, lmfdb, classlabel, clebsch, DeltaCl, fieldidx> */
for fi in [1..#s3] do
    r := s3[fi]; coeffs := eval r[3];
    f := &+[ coeffs[k+1]*a^k : k in [0..#coeffs-1] ];
    lab := IsDefined(nflabel, strip(r[3])) select nflabel[strip(r[3])] else "?";
    try
        fz := IntegralMonicPolynomial(f);
        P, rts, GalData := GaloisGroup(fz);
        gfull := GaloisSubgroup(GalData, sub<P|>); Kf := NumberField(gfull); dK := Degree(gfull);
        autos := Automorphisms(Kf); rtsK := [rr[1] : rr in Roots(gfull, Kf)];
        Sd := SymmetricGroup(dK);
        perms := [ Sd ! [ Index(rtsK, au(rtsK[i])) : i in [1..dK] ] : au in autos ];
        PK := sub<Sd|perms>;
        for cj in [1..NC] do
            Gwe6 := subs[S3idx[cj]][2];
            rhos := TwistEmbeddings(U, Gwe6, PK);
            Brho := DescendedBasisConcrete(Kf, autos, perms, rhos[1], bm : Print := false);
            Drho := SourceDescendedBasis(Kf, autos, perms, rhos[1], U`phi6 : Print := false);
            m := GenerateModuliPoints(Kf, bm, Brho, Drho : Count := 2, MaxTries := 2500);
            for kk in [1..Min(1,#m)] do
                ii := m[kk][3]; ok := true; gl := 0;
                try gl := inv_2_gl(ii); catch e ok := false; end try;
                if ok and IsSmoothCubicSurface(gl) then
                    Append(~cand, <gl, f, lab, subs[S3idx[cj]][1], ii, DeltaClebsch(ii), fi>);
                end if;
            end for;
        end for;
        printf "  field %o (%o): %o candidates so far\n", fi, lab, #cand;
    catch e printf "  field %o (%o): ERROR %o\n", fi, lab, e`Object; end try;
end for;

printf "\n%o candidate models; minimizing with a 120s timeout...\n", #cand;
raws := [ c[1] : c in cand ];
fin := ParallelMinimizeReduce(raws, dir : Timeout := 120, Print := true);
finidx := { t[1] : t in fin };
printf "minimized %o / %o\n", #finidx, #raws;

/* non-minimizers, smallest Delta_Cl first, preferring distinct fields */
nz := [ i : i in [1..#cand] | i notin finidx ];
Sort(~nz, func< i,j | #Sprint(Abs(cand[i][6])) - #Sprint(Abs(cand[j][6])) >);
chosen := []; seenf := {};
for pass in [1,2] do
    for i in nz do
        if i in chosen then continue; end if;
        if pass eq 1 and cand[i][7] in seenf then continue; end if;
        Include(~seenf, cand[i][7]); Append(~chosen, i);
        if #chosen ge 10 then break; end if;
    end for;
    if #chosen ge 10 then break; end if;
end for;
printf "chosen %o hard examples\n", #chosen;

out := "unminimized_examples.txt";
PrintFile(out,
 "# Cubic surfaces from the Elsenhans-Jahnel pipeline (this repo) that do NOT minimize\n" cat
 "# within a 120-second MinimizeReduceCubicSurface timeout.  All are S3 (3T2) fields.\n" cat
 "# The intrinsic discriminant Delta_Cl is large enough that the reconstructed model is\n" cat
 "# badly non-minimal: its model discriminant carries huge unfactorable spurious primes,\n" cat
 "# so MinimizeReduceCubicSurface stalls.  ReduceCubicSurface (LLL reduction only) shrinks\n" cat
 "# the coefficients by <2% here, so the size below reflects genuine non-minimality.\n" cat
 "#\n" cat
 "# Per example: the number field (its defining polynomial) and its LMFDB label, the W(E6) class (Galois action\n" cat
 "# on the 27 lines), the Clebsch invariants [I8,I16,I24,I32,I40] (which determine the\n" cat
 "# surface), Delta_Cl, and an explicit cubic in P^3 (x,y,z,w), ReduceCubicSurface-reduced\n" cat
 "# but NOT minimized.\n"
 : Overwrite := true);
n := 0;
for i in chosen do
    c := cand[i]; n +:= 1;
    glr := c[1]; try glr := ReduceCubicSurface(c[1]); catch e ; end try;
    dc := c[6];
    PrintFile(out, Sprintf("\n=== Example %o ===", n));
    PrintFile(out, Sprintf("number field : %o", c[2]));
    PrintFile(out, Sprintf("LMFDB        : %o  (https://www.lmfdb.org/NumberField/%o , Galois group S3 = 3T2)", c[3], c[3]));
    PrintFile(out, Sprintf("W(E6) class  : %o", c[4]));
    PrintFile(out, Sprintf("Clebsch      : %o", c[5]));
    PrintFile(out, Sprintf("Delta_Cl     : %o  (%o digits)", dc, #Sprint(Abs(Numerator(dc))*Denominator(dc))));
    PrintFile(out, Sprintf("max |coef|   : %o digits", Max([ #Sprint(Abs(co)) : co in Coefficients(glr) ])));
    PrintFile(out, Sprintf("cubic        : %o", glr));
end for;
printf "wrote %o with %o examples\n", out, n;
quit;
