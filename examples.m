/* examples.m

   Worked examples for cubic_surface_resolvent_twist.m (Elsenhans--Jahnel,
   Algorithm 5.1).  For each cyclic group C_n we build a smooth cubic surface
   over Q whose Galois action on the 27 lines is C_n, and certify it by
   computing that action directly from the 27 lines of the output surface.

   Each example calls CubicSurfaceFromPolynomial, which finds the W(E6)-subgroup
   classes isomorphic to Gal(f) and tries them until the rational-point search
   on the descended twist succeeds.  Run time is a few minutes per example.

      load "examples.m";
*/

load "cubic_surface_resolvent_twist.m";

Pq<t> := PolynomialRing(RationalField());

examples := [
    <"C4",  t^4 + t^3 + t^2 + t + 1>,                                       // Q(zeta_5)
    <"C5",  t^5 + t^4 - 4*t^3 - 3*t^2 + 3*t + 1>,                           // real Q(zeta_11)
    <"C9",  t^9 - 19*t^8 + 152*t^7 - 665*t^6 + 1729*t^5 - 2717*t^4
                + 2508*t^3 - 1254*t^2 + 285*t - 19>,                        // EJ subgroup no. 73
    <"C12", &+[t^i : i in [0..12]]>                                         // Q(zeta_13)
];

U := BuildUniversalCobleData(: Print := false);

for ex in examples do
    name := ex[1]; f := ex[2];
    printf "\n================ %o : Gal(f) = %o ================\n", name, name;
    tt := Cputime();
    R, ci := CubicSurfaceFromPolynomial(f : Universal := U, Print := false);
    printf "built in %o sec using subgroup class %o\n", Cputime(tt), ci;
    printf "Clebsch invariants : %o\n", R`invariants;
    printf "cubic surface      : %o\n", R`surface;

    P3<x0,x1,x2,x3> := ProjectiveSpace(RationalField(), 3);
    S := Surface(P3, Evaluate(PolynomialRing(RationalField(), 4)!R`surface, [x0,x1,x2,x3]));
    G := LinesGaloisGroupOfCubicSurface(R`surface : Print := false);
    printf "smooth over Q      : %o\n", IsNonsingular(S);
    printf "27-lines Galois grp: %o  (isomorphic to the prescribed group: %o)\n",
           GroupName(G), IsIsomorphic(G, R`subgroup);
end for;
