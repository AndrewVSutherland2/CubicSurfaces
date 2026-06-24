/* cubic_surface_resolvent_twist.m

   Implementation of Algorithm 5.1 of Elsenhans--Jahnel, "Moduli spaces and the
   inverse Galois problem for cubic surfaces" (arXiv:1209.5591): given a subgroup
   G of W(E6) and a polynomial f whose Galois group is isomorphic to G, it builds
   a smooth cubic surface over Q whose Galois action on the 27 lines realizes G.

   Main entry points:

      U := BuildUniversalCobleData();
      R := CubicSurfaceFromSubgroup(G, f : Universal := U);
      R := CubicSurfaceFromResolventPolynomial(classNo, f : Universal := U);

   Here f is a monic irreducible polynomial over Z or Q.  Magma's GaloisGroup is
   used to recognize Gal(f) and to produce a defining polynomial for the splitting
   field L; L is then realized as a concrete number field K, over which the two
   Galois descents (the twisted moduli space and the cubic equation) are carried
   out by exact arithmetic.  This is practical whenever [L:Q] = |G| is moderate.

   The optional Iso parameter is an explicit isomorphism Gal(L/Q) -> G, given as a
   homomorphism from the permutation realization of Aut(K) on the roots of the
   splitting-field polynomial.  If omitted, Magma's IsIsomorphic chooses one.

   Warning: the subgroup-class numbering used by CubicSurfaceFromResolventPolynomial
   is Magma's SubgroupClasses ordering for the 40-dimensional matrix realization
   constructed below.  It is not asserted to equal the GAP 4.4.12 numbering in
   Elsenhans--Jahnel.
*/

Q := RationalField();
Z := IntegerRing();

ring_z4<x,y,z,w> := PolynomialRing(IntegerRing(),4);

CobleDataRF := recformat<we6, bm, se, rel3, mon3, r10, mat_l>;
CubicTwistResultRF := recformat<
    surface,
    surface_unreduced,
    invariants,
    pentahedral_quintic,
    point,
    descended_cubics,
    basis_matrix,
    orbit_data,
    subgroup,
    galois_group,
    universal
>;

/* ------------------------------------------------------------------------- */
/* Equation problem: Clebsch invariants -> pentahedral quintic -> cubic.      */
/* ------------------------------------------------------------------------- */

penta_pol_2_gl := function(pol5)
    error if not (Degree(pol5) eq 5), "penta_pol_2_gl expects a degree-5 polynomial";
    error if not (Discriminant(pol5) ne 0), "The pentahedral quintic is not separable";

    alg := quo<Parent(pol5) | pol5>;

    nf_r4 := PolynomialRing(alg,4);
    bm0 := BasisMatrix(Kernel(Transpose(Matrix(
                         [[Trace(alg.1^i) : i in [0..4]]] ))));
    bas := [&+ [alg.1^(i-1) * bm0[j,i] : i in [1..5]] : j in [1..4]];
    lf := &+[bas[i] * nf_r4.i : i in [1..4]];
    kf := lf^3 * alg.1;
    gl := &+[Monomials(a)[1] * Trace(Coefficients(a)[1]) : a in Terms(kf)];
    gl := PolynomialRing(RationalField(),4)!gl;
    denom := LCM([Denominator(cc) : cc in Coefficients(gl)]);
    gl2 := ring_z4!(gl * denom);
    gl2 := gl2 div GCD(Coefficients(gl2));
    return gl2;
end function;

ClebschInvariantsToPentahedralPolynomial := function(inv2)
    error if not (#inv2 eq 5), "Need [A,B,C,D,E]";
    A := Q!inv2[1]; B := Q!inv2[2]; C := Q!inv2[3];
    D := Q!inv2[4]; E := Q!inv2[5];
    error if not (E ne 0), "E = 0: outside the proper pentahedral locus";

    /* Simplified form of Appendix A.8:
       [sigma_1,...,sigma_5] = [B, D, (C^2-AE)/4, CE, E^2]. */
    s1 := B;
    s2 := D;
    s3 := (C^2 - A*E)/4;
    s4 := C*E;
    s5 := E^2;

    P5<T> := PolynomialRing(Q);
    return T^5 - s1*T^4 + s2*T^3 - s3*T^2 + s4*T - s5;
end function;

inv_2_gl := function(inv2)
    pol5 := ClebschInvariantsToPentahedralPolynomial(inv2);
    return penta_pol_2_gl(pol5), pol5;
end function;

IsSmoothCubicSurface := function(gl)
    /* gl is a ternary... quaternary cubic form; the surface gl = 0 in P^3 is
       smooth iff its four partial derivatives have no common projective zero. */
    P3 := ProjectiveSpace(Q, 3);
    F := CoordinateRing(P3) ! gl;
    return IsEmpty(Scheme(P3, [Derivative(F, i) : i in [1..4]]));
end function;

/* ------------------------------------------------------------------------- */
/* Coble gamma coordinates and universal data.                               */
/* ------------------------------------------------------------------------- */

minor_ijk := function(mat,i,j,k)
    return Determinant(Submatrix(mat,[i,j,k],[1,2,3]));
end function;

gamma_3 := function(mat,i,j,k,l,m,n)
    return [minor_ijk(mat,i,k,l),  minor_ijk(mat,j,k,l),
            minor_ijk(mat,k,m,n),  minor_ijk(mat,l,m,n),
            minor_ijk(mat,m,i,j),  minor_ijk(mat,n,i,j)];
end function;

gamma_2 := function(mat, i,j,k,l,m,n)
    return [minor_ijk(mat,i,j,k), minor_ijk(mat,l,m,n),
            minor_ijk(mat,3,4,1)*minor_ijk(mat,5,6,1)*minor_ijk(mat,5,3,2)*minor_ijk(mat,4,6,2)
          - minor_ijk(mat,3,4,2)*minor_ijk(mat,5,6,2)*minor_ijk(mat,5,3,1)*minor_ijk(mat,4,6,1)];
end function;

all_gamma := function(mat)
    /* mat is the 6x3 matrix whose rows are the six blow-up points. */
    res1 := [&*gamma_3(mat, 1, i,j,k,l,m) : i,j,k,l,m in [2..6] |
             (#{i,j,k,l,m} eq 5) and (j lt k) and (l lt m)];
    res2 := [&*gamma_2(mat, 1, i,j,k,l,m) : i,j,k,l,m in [2..6] |
             (#{i,j,k,l,m} eq 5) and (i lt j) and (k lt l) and (l lt m)];
    return res1 cat res2;
end function;

RandomGeneralBlowupPointList := function()
    repeat
        a := [[1,0,0],[0,1,0],[0,0,1],[1,1,1],
              [Random([-100..100]),Random([-100..100]),Random([-100..100])],
              [Random([-100..100]),Random([-100..100]),Random([-100..100])]];
        if &or [&and [c eq 0 : c in a[i]] : i in [5,6]] then
            g := [Z|0];
        else
            g := all_gamma(Matrix(a));
        end if;
    until &and [u ne 0 : u in g];
    return a, g;
end function;

BuildUniversalCobleData := function(: InterpolationPoints := 300,
                                      RandomSeed := 1,
                                      Print := true)
    SetSeed(RandomSeed);

    if Print then printf "Evaluating Coble gamma coordinates for interpolation\n"; end if;
    bu_ll := [];
    gamma_ll := [];
    while #gamma_ll lt InterpolationPoints do
        a, g := RandomGeneralBlowupPointList();
        Append(~bu_ll, a);
        Append(~gamma_ll, g);
    end while;

    if Print then printf "Reconstructing the linear and cubic Coble relations\n"; end if;
    bm := BasisMatrix(RowSpace(Matrix(Q,gamma_ll)));
    error if not (Nrows(bm) eq 10), "Interpolation did not recover a 10-dimensional gamma span";

    r10 := PolynomialRing(Q,10);
    se := ElementToSequence(Vector([r10.i : i in [1..10]]) * ChangeRing(bm,r10));

    ul := [ElementToSequence(Solution(bm, Vector(Q,a))) : a in gamma_ll];
    mon3 := Monomials((&+[r10.i : i in [1..10]])^3);
    m3_mat := [[Evaluate(f,p) : f in mon3] : p in ul];
    ker := Kernel(Transpose(Matrix(m3_mat)));
    rel3_co := BasisMatrix(ker);
    rel3 := [&+[mon3[i] * rel3_co[k,i] : i in [1..#mon3]] : k in [1..Nrows(rel3_co)]];
    error if not (#rel3 eq 30), "Interpolation did not recover the 30 cubic relations";

    if Print then printf "Constructing the 40-dimensional signed W(E6) representation\n"; end if;
    bo :=  [[1,0,0],[0,1,0],[0,0,1],[1,1,1],[13,19,-31],[5,81,-31]];
    bo2 := [bo[i] : i in [1..4]] cat [[1/a : a in bo[i]] : i in [5,6]];
    bo3 := [bo[i] : i in [2,3,4,5,6,1]];
    bo4 := [bo[i] : i in [2,1,3,4,5,6]];

    g_ll_gen := [all_gamma(Matrix(a)) : a in [bo, bo2, bo3, bo4]];
    mul_l := [LCM([Denominator(b) : b in a]) : a in g_ll_gen];
    g_ll_gen := [[Z!(mul_l[i] * a) : a in g_ll_gen[i]] : i in [1..#g_ll_gen]];
    ggt_l := [GCD(a) : a in g_ll_gen];
    g_ll_gen := [[a div ggt_l[i] : a in g_ll_gen[i]] : i in [1..#g_ll_gen]];

    mat_l := [];
    for i := 2 to #g_ll_gen do
        mat := ZeroMatrix(Z,40,40);
        for j := 1 to 40 do
            k := Index(g_ll_gen[i], g_ll_gen[1][j]);
            if k eq 0 then k := Index(g_ll_gen[i], -g_ll_gen[1][j]); end if;
            error if not (k ne 0), "Could not recover a signed gamma permutation";
            mat[k,j] := g_ll_gen[1][j] div g_ll_gen[i][k];
        end for;
        error if not (IsRegular(mat)), "Recovered Weyl generator matrix is singular";
        Append(~mat_l, mat);
    end for;

    we6 := sub<GL(40,Z) | mat_l>;
    error if not (Order(we6) eq 51840), "The generated group is not W(E6) of order 51840";

    return rec<CobleDataRF | we6 := we6, bm := bm, se := se,
                            rel3 := rel3, mon3 := mon3, r10 := r10,
                            mat_l := mat_l>;
end function;

/* ------------------------------------------------------------------------- */
/* Small exact linear algebra helpers.                                       */
/* ------------------------------------------------------------------------- */

IntegralMonicPolynomial := function(f)
    Pq := PolynomialRing(Q);
    fq := Pq!f;
    error if not (IsMonic(fq)), "The input polynomial must be monic";
    error if not (&and [Denominator(c) eq 1 : c in Eltseq(fq)]), "The input polynomial must have integer coefficients";
    Pz := PolynomialRing(Z);
    fz := Pz!fq;
    error if not (IsIrreducible(fz)), "The input polynomial must be irreducible";
    return fz;
end function;

MatrixDenominator := function(M)
    den := Z!1;
    for i := 1 to Nrows(M) do
        for j := 1 to Ncols(M) do
            den := LCM(den, Denominator(Q!M[i,j]));
        end for;
    end for;
    return den;
end function;

LLLReduceRationalRowBasis := function(B)
    /* Return a short Z-basis of the rational row space of B.  Clearing
       denominators yields a sublattice that is in general imprimitive; saturating
       before LLL recovers the full integer lattice of the span, whose reduced
       basis is markedly smaller.  This matters a lot downstream: the descended
       cubics are formed in the coordinate system fixed by this basis, and the
       saturated basis shrinks their coefficients by orders of magnitude (e.g.
       from ~10^70 to ~10^31 for an A4 twist), which is what makes the subsequent
       rational-point search feasible. */
    den := MatrixDenominator(B);
    BZ := ChangeRing(den * ChangeRing(B, Q), Z);
    return ChangeRing(LLL(Saturation(BZ)), Q);
end function;

ExpCode := function(exp)
    /* Valid for total degree <= 3.  Base 4 prevents carries. */
    return &+[Z!exp[i] * 4^(i-1) : i in [1..#exp]];
end function;

MonomialCodeIndex := function(mon3)
    A := AssociativeArray(Integers());
    for i := 1 to #mon3 do
        A[ExpCode([Degree(mon3[i],j) : j in [1..10]])] := i;
    end for;
    return A;
end function;

/* ------------------------------------------------------------------------- */
/* Action of W(E6) on the ten Coble coordinates.                             */
/* ------------------------------------------------------------------------- */

WeylActionOnYCoords := function(Mg, bmQ, GramInv)
    /* Mg is a 40x40 signed permutation matrix realizing an element g of W(E6)
       on the forty gamma coordinates (acting on the gamma column by Mg).  The
       gamma linear span is g-stable, so there is a unique 10x10 rational matrix
       R_g with  bm^T R_g = Mg bm^T  on that span; equivalently it describes the
       action of g on the ten Coble coordinates y, where gamma = bm^T y.  We
       recover it by least squares against the (rational) Gram inverse of bm. */
    return GramInv * bmQ * ChangeRing(Mg, Q) * Transpose(bmQ);
end function;

WeylSubgroupByClassNumber := function(we6, classNo)
    cl := SubgroupClasses(we6);
    error if not (1 le classNo and classNo le #cl), "classNo outside Magma's SubgroupClasses list";
    return cl[classNo]`subgroup;
end function;

/* ------------------------------------------------------------------------- */
/* Twisted Galois descent over the concrete splitting field.                 */
/*                                                                            */
/* Let L be the splitting field of f, realized as a concrete number field K,  */
/* with rho : Gal(L/Q) -> G the prescribed isomorphism onto the chosen        */
/* W(E6)-subgroup.  A Q-rational point of the twisted moduli space is a vector */
/* y in the ten Coble coordinates over L satisfying                           */
/*                                                                            */
/*        R_{rho(sigma)} . sigma(y) = y       for all sigma in Gal(L/Q).      */
/*                                                                            */
/* Writing y_j = sum_b U_{bj} omega_b in a fixed Q-basis omega of L turns this */
/* into a rational linear system in the 10[L:Q] unknowns U; its solution space */
/* is ten dimensional and a basis of it is the descended basis B_rho.         */
/* ------------------------------------------------------------------------- */

DescendedBasisConcrete := function(K, autos, perms, psi, bm : LLLReduce := true, Print := true)
    d  := Degree(K);
    aa := K.1;
    bmQ     := ChangeRing(bm, Q);
    GramInv := (bmQ * Transpose(bmQ))^-1;
    N  := 10*d;
    Sd := Universe(perms);
    PK := sub<Sd | perms>;

    /* one semilinear constraint  R_g . sigma(y) = y  per generator suffices */
    gens := [ Sd | g : g in Generators(PK) | g ne Id(Sd) ];

    constraintRows := [];
    for g in gens do
        idx := Index(perms, g);
        error if not (idx ne 0), "Could not match a Galois generator to an automorphism";
        /* Magma labels an automorphism by the permutation it induces on the
           roots, and that labelling is an ANTI-homomorphism (since (i^h)^g =
           i^(hg)).  To make sigma |-> (its action) a genuine homomorphism -- so
           that the maps y |-> R_{rho(sigma)} sigma(y) form a group action and the
           fixed space is ten-dimensional -- we pair the generator with the
           inverse automorphism.  (For abelian G this is immaterial.) */
        au := Inverse(autos[idx]);
        Rg := WeylActionOnYCoords(psi(g), bmQ, GramInv);
        /* S^au : column b holds au(a^(b-1)) in the power basis 1,a,...,a^(d-1) */
        Sa := Transpose(Matrix(Q, d, d, [ Eltseq(au(aa)^(b-1)) : b in [1..d] ]));
        /* Operator A on vec(U), column-stacked as vec[(b-1)*10+j] = U_{bj}:
              (A.vec)[(c-1)*10+i] = sum_{b,j} Sa[c,b] Rg[i,j] U_{bj}.
           The constraint is A.vec = vec, i.e. (A - I).vec = 0. */
        AmI := ZeroMatrix(Q, N, N);
        for b := 1 to d do
            for j := 1 to 10 do
                p := (b-1)*10 + j;
                for c := 1 to d do
                    for i := 1 to 10 do
                        AmI[(c-1)*10+i, p] +:= Sa[c,b]*Rg[i,j];
                    end for;
                end for;
            end for;
        end for;
        AmI := AmI - IdentityMatrix(Q, N);
        constraintRows cat:= [AmI[r] : r in [1..N]];
    end for;

    fixed := #constraintRows eq 0 select VectorSpace(Q, N)
                                   else  Kernel(Transpose(Matrix(constraintRows)));
    error if not (Dimension(fixed) eq 10),
        Sprintf("Twisted fixed space has dimension %o, not 10", Dimension(fixed));

    BF := BasisMatrix(fixed);
    if LLLReduce then
        /* a reduced Z-basis of the fixed lattice gives small-height coordinates,
           which makes the rational-point search much easier (cf. EJ, Sec. 5). */
        BF := LLLReduceRationalRowBasis(BF);
    end if;

    /* each fixed vector becomes a y-vector y in K^10 */
    Brho := [];
    for r := 1 to 10 do
        Append(~Brho, [ K![ BF[r][(b-1)*10 + j] : b in [1..d] ] : j in [1..10] ]);
    end for;
    return Brho;
end function;

/* ------------------------------------------------------------------------- */
/* Descending the cubic relation space.                                       */
/*                                                                            */
/* The thirty universal cubics F in the ten Coble coordinates y are pulled    */
/* back along y = B_rho z, where B_rho has entries in the splitting field K.   */
/* Each resulting cubic in z has K-coefficients; splitting those coefficients  */
/* over the power basis of K and row-reducing yields the thirty rational       */
/* cubics that cut out the twist.  All of this is exact arithmetic over K.     */
/* ------------------------------------------------------------------------- */

DescendCubicsConcrete := function(rel3, Brho, K, bm, mon3 : Print := true)
    d := Degree(K);

    /* y_j as a K-linear form in the ten rational coordinates z */
    RK := PolynomialRing(K, 10);
    Ly := [ &+[ Brho[a][j]*RK.a : a in [1..10] ] : j in [1..10] ];

    codeToIdx := MonomialCodeIndex(mon3);
    nmon := #mon3;

    all_rows := [];
    for k := 1 to #rel3 do
        if Print then printf "  descending cubic relation %o/%o\n", k, #rel3; end if;
        Fk := Evaluate(rel3[k], Ly);
        rows := [ [Q!0 : i in [1..nmon]] : t in [1..d] ];
        mons := Monomials(Fk);
        cofs := Coefficients(Fk);
        for u := 1 to #mons do
            code := ExpCode([Degree(mons[u],j) : j in [1..10]]);
            idx := codeToIdx[code];
            cs := Eltseq(K!cofs[u]);
            for t := 1 to d do
                rows[t][idx] := cs[t];
            end for;
        end for;
        all_rows cat:= rows;
    end for;

    co_mat := Matrix(Q, all_rows);
    co_mat_2 := BasisMatrix(RowSpace(co_mat));
    error if not (Nrows(co_mat_2) eq 30),
        Sprintf("Descended cubic relation space has dimension %o, not 30", Nrows(co_mat_2));

    /* The substitution and row reduction leave the thirty cubics with enormous
       integer coefficients (growing fast with the field degree -- e.g. ~10^945
       at degree twelve), which cripples every later Groebner/point computation.
       Replacing the row space basis by an LLL-reduced basis of the saturated
       integer coefficient lattice recovers genuinely small representatives of
       the same thirty-dimensional space (~10^70 in that example) and speeds the
       rational-point search by orders of magnitude. */
    den := LCM([ Denominator(co_mat_2[i,j]) : i in [1..Nrows(co_mat_2)], j in [1..nmon] ]);
    red := LLL(Saturation(ChangeRing(den*co_mat_2, Z)));

    return [&+[red[i,j]*mon3[j] : j in [1..nmon]] : i in [1..Nrows(red)]];
end function;

/* ------------------------------------------------------------------------- */
/* Point search and invariant recovery.                                      */
/* ------------------------------------------------------------------------- */

SearchRationalPoint := function(kub_des : MaxSlices := 210, RandomSlices := 300,
                                          RandomBound := 2, RandomSeed := 1,
                                          Validate := false, Print := true)
    /* The descended twist is a rational fourfold in P^9, so a naive bounded
       search over P^9 is hopeless.  Instead we cut it down to dimension zero by
       intersecting with dim(X) hyperplanes and solve the resulting finite scheme
       exactly; its rational points are produced directly and tested against the
       validation predicate.  We use two complementary families of slices:
       coordinate hyperplanes, which tend to expose the special rational points
       sitting on coordinate subspaces, and random small-coefficient hyperplanes,
       which cut through the interior where the generic points live.  This is the
       heuristic step of Algorithm 5.1; for a stubborn twist one may also supply a
       point directly or pass a different generating polynomial. */
    hasValidate := Type(Validate) ne BoolElt;

    PP := ProjectiveSpace(Q, 9);
    R  := Parent(kub_des[1]);
    CR := CoordinateRing(PP);
    h  := hom< R -> CR | [PP.i : i in [1..10]] >;
    cub := [ h(F) : F in kub_des ];
    X := Scheme(PP, cub);
    dimX := Dimension(X);
    error if not (dimX ge 0), "Descended twist is empty";

    tryHyperplanes := function(hs)
        Y := Scheme(PP, cub cat hs);
        if Dimension(Y) ne 0 then return false, []; end if;
        for p in RationalPoints(Y) do
            sol := Eltseq(p);
            if (not hasValidate) or Validate(sol) then return true, sol; end if;
        end for;
        return false, [];
    end function;

    if Print then printf "Searching coordinate slices of the descended twist (dim %o)\n", dimX; end if;
    tried := 0;
    for S in Subsets({1..10}, dimX) do
        tried +:= 1;
        if tried gt MaxSlices then break; end if;
        ok, sol := tryHyperplanes([ CR.i : i in S ]);
        if ok then
            if Print then printf "  usable point on coordinate slice %o\n", S; end if;
            return true, sol;
        end if;
    end for;

    if Print then printf "Searching %o random hyperplane slices\n", RandomSlices; end if;
    SetSeed(RandomSeed);
    for iter := 1 to RandomSlices do
        hs := [ &+[ Random([-RandomBound..RandomBound])*CR.j : j in [1..10] ] : k in [1..dimX] ];
        ok, sol := tryHyperplanes(hs);
        if ok then
            if Print then printf "  usable point on random slice %o\n", iter; end if;
            return true, sol;
        end if;
    end for;

    if Print then printf "  no usable rational point found\n"; end if;
    return false, [];
end function;

ClebschInvariantsConcrete := function(sol, Brho, K, bm : Print := true)
    /* From a rational point z = sol on the descended twist, reconstruct the
       forty gamma coordinates over K, form the even power sums
       P_2k = sum_{i=1}^{40} gamma_i^{2k}, which are W(E6)-invariant and hence
       Galois-fixed (rational), and assemble Clebsch's invariant vector through
       Theorem 3.9 / formulas (invvsinv) of Elsenhans--Jahnel. */
    y   := [ &+[ (Q!sol[a]) * Brho[a][j] : a in [1..10] ] : j in [1..10] ];
    gam := [ &+[ y[t] * bm[t,jj] : t in [1..10] ] : jj in [1..40] ];

    powerSum := function(e)
        val := &+[ gam[jj]^e : jj in [1..40] ];
        ok, r := IsCoercible(Q, val);
        error if not (ok),
            Sprintf("Power sum P_%o is not Galois-fixed; the point is not on the twist", e);
        return r;
    end function;

    if Print then printf "Computing the even gamma power sums\n"; end if;
    P_2  := powerSum(2);
    P_4  := powerSum(4);
    P_6  := powerSum(6);
    P_8  := powerSum(8);
    P_10 := powerSum(10);

    i8 := (-6) * P_2;
    i16 := -24 * P_4  + (41/16) * P_2^2;
    i24 := 576/13 * P_6 -396/13 * P_4 * P_2 + 29/13 * P_2^3;
    i32 := -62208/1171 *P_8 + 54864/1171*P_6 * P_2
              + 203616/1171 * P_4^2 -61287/1171 * P_4 * P_2^2
              + 13393/4684 * P_2^4;

    i40 := 41472/155 * P_10 -4605984/36301 *P_8 * P_2
              -106272/403 * P_6 * P_4 +19990440/471913 * P_6 * P_2^2
           + 47719206/471913 * P_4^2 * P_2 -7468023/471913 * P_4 * P_2^3
           + 10108327/18876520 * P_2^5;

    return [i8,i16,i24,i32,i40];
end function;

/* ------------------------------------------------------------------------- */
/* Independent certification: the Galois action on the 27 lines.             */
/* ------------------------------------------------------------------------- */

LinesGaloisGroupOfCubicSurface := function(gl : Tries := 12, Print := true)
    /* For the smooth cubic surface gl = 0, compute the Galois group of the
       field of definition of its 27 lines, returned as a permutation group on
       the 27 lines.  A line that is a graph over the (x_0 : x_1) ruling, namely
       x_2 = a x_0 + b x_1, x_3 = c x_0 + d x_1, lies on the surface iff the
       binary cubic obtained by substitution vanishes identically; the resulting
       four equations cut out a zero-dimensional scheme of degree 27 whenever the
       chart captures every line.  We apply random linear coordinate changes
       until that happens, then read off the Galois group from the degree-27
       resolvent of a generic linear form in (a,b,c,d).  A single coordinate such
       as a often generates only a subfield of the field of definition of the
       lines, returning a proper quotient of the true Galois group; a generic
       linear form generates the whole field, so the resolvent's splitting field
       is the full field and GaloisGroup returns the correct group acting on the
       27 lines. */
    R4 := PolynomialRing(Q, 4);
    glQ := R4 ! gl;
    PT<T> := PolynomialRing(Q);

    for attempt := 1 to Tries do
        if attempt eq 1 then
            g := glQ;
        else
            repeat
                M := Matrix(Q, 4, 4, [Random([-3..3]) : i in [1..16]]);
            until Determinant(M) ne 0;
            g := Evaluate(glQ, [&+[M[i,j]*R4.j : j in [1..4]] : i in [1..4]]);
        end if;

        A4 := PolynomialRing(Q, 4);
        Pxy<X,Y> := PolynomialRing(A4, 2);
        Fb := Evaluate(g, [X, Y, A4.1*X + A4.2*Y, A4.3*X + A4.4*Y]);
        cfs := [ MonomialCoefficient(Fb, X^k*Y^(3-k)) : k in [0..3] ];
        Sch := Scheme(AffineSpace(A4), cfs);

        if Dimension(Sch) eq 0 and Degree(Sch) eq 27 then
            /* minimal polynomial of a generic linear form u = a + r2 b + r3 c + r4 d
               on the 27 line-points, by elimination of (a,b,c,d). */
            A5 := PolynomialRing(Q, 5);
            emb := hom< A4 -> A5 | [A5.1, A5.2, A5.3, A5.4] >;
            lam := A5.1 + Random([2..4])*A5.2 + Random([5..9])*A5.3 + Random([10..16])*A5.4;
            Jel := EliminationIdeal(ideal< A5 | [emb(F) : F in cfs] cat [A5.5 - lam] >, 4);
            gens := [gg : gg in Basis(Jel) | gg ne 0];
            if #gens ge 1 then
                fa := SquarefreePart(Evaluate(gens[1], [0, 0, 0, 0, T]));
                if Degree(fa) eq 27 then
                    if Print then printf "  all 27 lines captured; computing their Galois group\n"; end if;
                    return GaloisGroup(fa);
                end if;
            end if;
        end if;
    end for;

    error "Could not capture all 27 lines in a single chart; increase Tries";
end function;

/* ------------------------------------------------------------------------- */
/* Main functions.                                                           */
/* ------------------------------------------------------------------------- */

CubicSurfaceFromSubgroup := function(Gwe6, f : Universal := false,
                                             Iso := false,
                                             HeightBound := 2,
                                             Point := false,
                                             RandomSlices := 300,
                                             RandomBound := 2,
                                             RandomSeed := 1,
                                             ProveGalois := false,
                                             ProofLinearRelations := false,
                                             LLLReduceBasis := true,
                                             ReturnTwistOnly := false,
                                             DoMinimizeReduce := true,
                                             Print := true)
    Udata := Universal;
    if Type(Udata) eq BoolElt then
        Udata := BuildUniversalCobleData(: Print := Print);
    end if;

    fz := IntegralMonicPolynomial(f);

    if Print then printf "Computing Galois group of input polynomial\n"; end if;
    P, roots, GalData := GaloisGroup(fz);
    if ProveGalois then
        error if not (GaloisProof(fz, GalData)), "GaloisProof failed";
    end if;

    error if not (Order(P) eq Order(Gwe6)), Sprintf("Order(Gal(f)) = %o, but Order(G) = %o", Order(P), Order(Gwe6));

    /* Realise the splitting field concretely as K = Q[x]/(gfull), with gfull
       the resolvent of the trivial subgroup, and read off its automorphisms. */
    if Print then printf "Constructing the splitting field and its automorphisms\n"; end if;
    gfull := GaloisSubgroup(GalData, sub<P | >);
    K := NumberField(gfull);
    dK := Degree(gfull);
    if Print then printf "  splitting field has degree %o over Q\n", dK; end if;

    autos := Automorphisms(K);
    rtsK := [r[1] : r in Roots(gfull, K)];
    Sd := SymmetricGroup(dK);
    perms := [ Sd ! [ Index(rtsK, au(rtsK[i])) : i in [1..dK] ] : au in autos ];
    PK := sub<Sd | perms>;

    if Type(Iso) eq BoolElt then
        ok, psi := IsIsomorphic(PK, Gwe6);
        error if not (ok), "Gal(f) is not isomorphic to the selected W(E6)-subgroup";
    else
        psi := Iso;
    end if;

    if Print then printf "Solving the twisted Galois descent for the 10-dimensional fixed space\n"; end if;
    Brho := DescendedBasisConcrete(K, autos, perms, psi, Udata`bm
                                   : LLLReduce := LLLReduceBasis, Print := Print);

    if Print then printf "Descending the 30 cubic relations\n"; end if;
    kub_des := DescendCubicsConcrete(Udata`rel3, Brho, K, Udata`bm, Udata`mon3 : Print := Print);

    if ReturnTwistOnly then
        return rec<CubicTwistResultRF |
            surface := false,
            surface_unreduced := false,
            invariants := false,
            pentahedral_quintic := false,
            point := false,
            descended_cubics := kub_des,
            basis_matrix := Brho,
            orbit_data := false,
            subgroup := Gwe6,
            galois_group := P,
            universal := Udata>;
    end if;

    /* A point on the projective closure of the twist need not lie in the open
       locus of smooth marked surfaces with a proper pentahedron.  We reject any
       candidate whose Clebsch vector has E = 0 (improper pentahedron), whose
       pentahedral quintic is inseparable (an Eckardt point / nontrivial
       automorphism), or whose reconstructed cubic surface is singular (the
       discriminant may vanish even when E != 0), and keep searching, exactly as
       prescribed by Algorithm 5.1. */
    okClebsch := function(s)
        ii := ClebschInvariantsConcrete(s, Brho, K, Udata`bm : Print := false);
        if ii[5] eq 0 then return false; end if;
        pol5 := ClebschInvariantsToPentahedralPolynomial(ii);
        if Discriminant(pol5) eq 0 then return false; end if;
        return IsSmoothCubicSurface(penta_pol_2_gl(pol5));
    end function;

    if Type(Point) eq BoolElt then
        okpt, sol := SearchRationalPoint(kub_des : RandomSlices := RandomSlices,
                                                  RandomBound := RandomBound,
                                                  RandomSeed := RandomSeed,
                                                  Validate := okClebsch, Print := Print);
        error if not (okpt),
            "No usable rational point found by the heuristic search "
            * "(supply Point, raise RandomSlices, or try another polynomial)";
    else
        sol := [Q!a : a in Point];
        error if not (#sol eq 10), "Point must have 10 coordinates";
        error if not (&and [Evaluate(F, sol) eq 0 : F in kub_des]), "The supplied point is not on the descended twist";
    end if;

    if Print then printf "Using descended point %o\n", sol; end if;
    invs := ClebschInvariantsConcrete(sol, Brho, K, Udata`bm : Print := Print);

    if Print then printf "Clebsch invariants: %o\n", invs; end if;
    gl, pol5 := inv_2_gl(invs);
    gl2 := gl;

    if DoMinimizeReduce then
        if Print then printf "Trying MinimizeReduceCubicSurface\n"; end if;
        try
            gl2 := MinimizeReduceCubicSurface(gl);
        catch e
            if Print then printf "  MinimizeReduceCubicSurface failed; returning unreduced model\n"; end if;
            gl2 := gl;
        end try;
    end if;

    return rec<CubicTwistResultRF |
        surface := gl2,
        surface_unreduced := gl,
        invariants := invs,
        pentahedral_quintic := pol5,
        point := sol,
        descended_cubics := kub_des,
        basis_matrix := Brho,
        orbit_data := false,
        subgroup := Gwe6,
        galois_group := P,
        universal := Udata>;
end function;

CubicSurfaceFromResolventPolynomial := function(classNo, f : Universal := false,
                                                           Iso := false,
                                                           HeightBound := 2,
                                                           Point := false,
                                                           RandomSlices := 300,
                                                           RandomBound := 2,
                                                           RandomSeed := 1,
                                                           ProveGalois := false,
                                                           ProofLinearRelations := false,
                                                           LLLReduceBasis := true,
                                                           ReturnTwistOnly := false,
                                                           DoMinimizeReduce := true,
                                                           Print := true)
    Udata := Universal;
    if Type(Udata) eq BoolElt then
        Udata := BuildUniversalCobleData(: Print := Print);
    end if;
    Gwe6 := WeylSubgroupByClassNumber(Udata`we6, classNo);
    return CubicSurfaceFromSubgroup(Gwe6, f : Universal := Udata,
                                             Iso := Iso,
                                             HeightBound := HeightBound,
                                             Point := Point,
                                             RandomSlices := RandomSlices,
                                             RandomBound := RandomBound,
                                             RandomSeed := RandomSeed,
                                             ProveGalois := ProveGalois,
                                             ProofLinearRelations := ProofLinearRelations,
                                             LLLReduceBasis := LLLReduceBasis,
                                             ReturnTwistOnly := ReturnTwistOnly,
                                             DoMinimizeReduce := DoMinimizeReduce,
                                             Print := Print);
end function;

CubicSurfaceFromPolynomial := function(f : Universal := false,
                                           MaxClasses := 6,
                                           RandomSlices := 300,
                                           RandomBound := 2,
                                           Print := true)
    /* Convenience driver: realise the Galois group of f on the 27 lines without
       naming a W(E6)-subgroup.  A given abstract group can sit inside W(E6) in
       several non-conjugate ways, and these give genuinely different twists --
       different points of the inverse Galois problem -- whose rational-point
       searches differ widely in difficulty.  We therefore try the matching
       subgroup classes in turn until one yields a smooth surface. */
    Udata := Universal;
    if Type(Udata) eq BoolElt then
        Udata := BuildUniversalCobleData(: Print := Print);
    end if;
    fz := IntegralMonicPolynomial(f);
    Pgal := GaloisGroup(fz);
    cl := SubgroupClasses(Udata`we6);
    cand := [ c`subgroup : c in cl |
              Order(c`subgroup) eq Order(Pgal) and IsIsomorphic(c`subgroup, Pgal) ];
    error if not (#cand ge 1), "No subgroup of W(E6) is isomorphic to Gal(f)";
    if Print then printf "Gal(f) matches %o subgroup class(es) of W(E6)\n", #cand; end if;

    for ci := 1 to Min(MaxClasses, #cand) do
        if Print then printf "Trying subgroup class %o of %o\n", ci, #cand; end if;
        try
            return CubicSurfaceFromSubgroup(cand[ci], fz : Universal := Udata,
                                                          RandomSlices := RandomSlices,
                                                          RandomBound := RandomBound,
                                                          Print := Print), ci;
        catch e
            if Print then printf "  class %o gave no usable surface (%o)\n", ci, e`Object; end if;
        end try;
    end for;
    error Sprintf("No usable surface found among the first %o of %o subgroup classes",
                  Min(MaxClasses, #cand), #cand);
end function;

RunC9Example := function(: HeightBound := 2,
                           ReturnTwistOnly := false,
                           Print := true)
    U := BuildUniversalCobleData(: Print := Print);
    cl := SubgroupClasses(U`we6);
    c9s := [a`subgroup : a in cl | IsCyclic(a`subgroup) and Order(a`subgroup) eq 9];
    error if not (#c9s ge 1), "No cyclic subgroup of order 9 found";
    c9 := c9s[1];

    Pq<t> := PolynomialRing(Q);
    pol9 := t^9 - 19*t^8 + 152*t^7 - 665*t^6 + 1729*t^5 - 2717*t^4
              + 2508*t^3 - 1254*t^2 + 285*t - 19;

    return CubicSurfaceFromSubgroup(c9, pol9 : Universal := U,
                                             HeightBound := HeightBound,
                                             ReturnTwistOnly := ReturnTwistOnly,
                                             Print := Print);
end function;

/*
Example session (the cyclic-of-order-9 example, subgroup no. 73 of Elsenhans--Jahnel):

    load "cubic_surface_resolvent_twist.m";
    R := RunC9Example();
    R`surface;
        x^2*y - x^2*w - 3*x*y*w - x*z^2 + 3*x*z*w + x*w^2
            - y^2*z - y*z^2 + y*z*w + y*w^2 + z*w^2

    // independent certification: the 27 lines and their Galois group
    G := LinesGaloisGroupOfCubicSurface(R`surface);
    #G, IsCyclic(G);           // 9 true   -- the action on the 27 lines is C9

Without naming a subgroup, CubicSurfaceFromPolynomial just takes a polynomial and
realises its Galois group on the 27 lines, trying the matching W(E6)-subgroup
classes in turn (a given abstract group embeds in W(E6) in several inequivalent
ways, giving genuinely different twists):

    U  := BuildUniversalCobleData();
    Pq<t> := PolynomialRing(RationalField());
    R  := CubicSurfaceFromPolynomial(t^5 + t^4 - 4*t^3 - 3*t^2 + 3*t + 1 : Universal := U);
    LinesGaloisGroupOfCubicSurface(R`surface);   // C5

The following all run end to end and are certified to realise the named group on
the 27 lines (each picks the first workable subgroup class):

    C4 :  t^4 + t^3 + t^2 + t + 1
    C5 :  t^5 + t^4 - 4*t^3 - 3*t^2 + 3*t + 1
    C9 :  t^9 - 19*t^8 + ... + 285*t - 19          (RunC9Example, subgroup no. 73)
    C12:  t^12 + t^11 + ... + t + 1

To drive the descent for an explicit subgroup, call CubicSurfaceFromSubgroup(G, f)
with G a subgroup of U`we6, or CubicSurfaceFromResolventPolynomial(classNo, f) to
select it by Magma's SubgroupClasses index.

The final step -- searching for a rational point on the descended twist -- is the
heuristic part of the algorithm (cf. EJ, Section 5).  The search solves exact
coordinate- and random-hyperplane slices of the (rational, four-dimensional) twist.
The substituted cubics would otherwise carry enormous coefficients (growing with the
field degree), so the descent reduces them twice -- a saturated LLL basis for the
fixed space and an LLL-reduced integer basis for the thirty cubics -- which keeps the
slice solves cheap and is what lets the cyclic examples above run in seconds.

Difficulty nonetheless varies enormously between twists.  Cyclic groups behave well;
for several non-cyclic and very small groups (e.g. C2, C3, S3, D4) every accessible
rational point of the twist is degenerate (a singular surface), and no smooth point
was reached by slicing, bounded search, or line walking -- consistent with EJ, who
handled those "easy" conjugacy classes by direct search and reserved this descent for
the seven hard classes (no. 73 = C9 among them).  For a stubborn twist one may raise
RandomSlices, supply a known point through the Point parameter, try another subgroup
class, or (as the paper recommends) repeat with a different polynomial f.
*/
