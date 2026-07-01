/* bootstrap_gaps.m -- extract source fields for the W(E6) classes that have
   no polynomial in WE6fields.txt, from the Elsenhans-Jahnel surfaces
   classified to those classes (classify_ej.m).

   Every gap class has a faithful line orbit of size exactly its designated
   degree (24 or 27), so the lines polynomial of the EJ surface in that class
   has an irreducible factor of that degree whose Galois closure realizes the
   class group with the correct 27-line embedding available.  The extracted
   field only needs Gal(f) isomorphic to the class group -- the no-split
   pipeline (which certifies its own output) does the rest, so a misclassified
   EJ surface is caught downstream, not propagated.

   Output: gap_fields.txt with rows  label:[coeffs]  (constant first). */
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "cubic_surface_nosplit.m";
Q := Rationals();
Pxyzw<x,y,z,w> := PolynomialRing(Q, 4);
load "ElsenhansJahnelDatabase.txt";
TMP := GetEnvironmentValue("BGTMP");
OUT := "gap_fields.txt";

/* degree-27 lines polynomial of a cubic surface (the elimination from
   LinesGaloisGroupOfCubicSurface, returning the polynomial itself) */
LinesPolynomial := function(gl : Tries := 12)
    R4 := PolynomialRing(Q, 4);
    glQ := R4 ! gl;
    PT<T> := PolynomialRing(Q);
    for attempt := 1 to Tries do
        if attempt eq 1 then
            g := glQ;
        else
            repeat M := Matrix(Q, 4, 4, [Random([-3..3]) : i in [1..16]]);
            until Determinant(M) ne 0;
            g := Evaluate(glQ, [&+[M[i,j]*R4.j : j in [1..4]] : i in [1..4]]);
        end if;
        A4 := PolynomialRing(Q, 4);
        Pxy<X,Y> := PolynomialRing(A4, 2);
        Fb := Evaluate(g, [X, Y, A4.1*X + A4.2*Y, A4.3*X + A4.4*Y]);
        cfs := [ MonomialCoefficient(Fb, X^k*Y^(3-k)) : k in [0..3] ];
        Sch := Scheme(AffineSpace(A4), cfs);
        if Dimension(Sch) eq 0 and Degree(Sch) eq 27 then
            A5 := PolynomialRing(Q, 5);
            emb := hom< A4 -> A5 | [A5.1, A5.2, A5.3, A5.4] >;
            lam := A5.1 + Random([2..4])*A5.2 + Random([5..9])*A5.3 + Random([10..16])*A5.4;
            Jel := EliminationIdeal(ideal< A5 | [emb(F) : F in cfs] cat [A5.5 - lam] >, 4);
            gens := [gg : gg in Basis(Jel) | gg ne 0];
            if #gens ge 1 then
                fa := SquarefreePart(Evaluate(gens[1], [0, 0, 0, 0, T]));
                if Degree(fa) eq 27 then return true, fa; end if;
            end if;
        end if;
    end for;
    return false, PT!0;
end function;

/* label -> <EJ index, faithful orbit degree> */
jobs := [
    <"51840.b.80.b1.a1",   346, 27>,
    <"51840.b.90.a1.a1",   324, 24>,
    <"51840.b.90.b1.a1",   322, 24>,
    <"51840.b.90.c1.a1",   323, 24>,
    <"51840.b.120.b1.a1",  342, 27>,
    <"51840.b.180.a1.a1",  320, 24>,
    <"51840.b.240.h1.a1",  339, 27>,
    <"51840.b.240.i1.a1",  337, 27>,
    <"51840.b.270.d1.a1",  319, 24>
];

for job in jobs do
    label := job[1]; idx := job[2]; dwant := job[3];
    printf "== %o from EJ[%o], want degree %o\n", label, idx, dwant;
    ok, fa := LinesPolynomial(ElsenhansJahnelDatabase[idx]);
    if not ok then printf "  lines polynomial FAILED\n"; continue; end if;
    facs := [ fc[1] : fc in Factorization(fa) | Degree(fc[1]) eq dwant ];
    if #facs eq 0 then
        printf "  no degree-%o factor (degrees %o)\n", dwant,
            [ Degree(fc[1]) : fc in Factorization(fa) ];
        continue;
    end if;
    red := PolredbestResolvents([ facs[1] ], TMP);
    g := red[1][1];
    printf "  field (%o digits max coeff): %o\n",
        Max([ #Sprint(Abs(Numerator(c))) : c in Coefficients(g) ]), g;
    PrintFile(OUT, Sprintf("%o:%o", label,
        [ Numerator(c) : c in Coefficients(g) ]));
end for;
printf "bootstrap done\n";
quit;
