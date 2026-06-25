/* database_pipeline.m

   A generation pipeline on top of cubic_surface_resolvent_twist.m for building a
   database of explicit cubic surfaces.  Given a polynomial f defining a number
   field with Galois group G, then for each conjugacy class of W(E6)-subgroups
   isomorphic to G and each isomorphism rho : G -> subgroup, it:

     1. generates many distinct moduli points via the A^6 cuspidal map (cheap;
        no reduction),
     2. ranks them by the Clebsch--Salmon discriminant  Delta_Cl  (the intrinsic
        discriminant, computed from the Clebsch invariants -- EJ Lemma 2.6),
     3. reconstructs the lowest-Delta_Cl candidates and minimizes them in parallel
        (spawned Magma jobs, each under a timeout, so an unfactorable model
        discriminant cannot stall the batch),
     4. keeps the smallest explicit cubics.

   Load AFTER the core file:
       load "cubic_surface_resolvent_twist.m";
       load "database_pipeline.m";
*/

/* ------------------------------------------------------------------------- */
/* The Clebsch--Salmon discriminant from the Clebsch invariants.             */
/* (A,B,C,D,E) = (I8,I16,I24,I32,I40);  E does not appear (EJ Lemma 2.6).     */
/* ------------------------------------------------------------------------- */

DeltaClebsch := function(ii)
    A := ii[1]; B := ii[2]; C := ii[3]; D := ii[4];
    return A^4 - 128*A^2*B + 4096*B^2 - 2048*A*C - 16384*D;
end function;

IntegerSize := function(x)
    /* number of decimal digits of the integer carrying all primes of rational x */
    return #Sprint(Abs(Numerator(x)) * Denominator(x));
end function;

/* ------------------------------------------------------------------------- */
/* Generate distinct moduli points for one twist (K, Brho, Drho given).       */
/* Returns a list of <DeltaClSize, z, clebsch, DeltaCl>, sorted by DeltaCl     */
/* size (smallest intrinsic discriminant first).                              */
/* ------------------------------------------------------------------------- */

GenerateModuliPoints := function(Kf, bm, Brho, Drho : Count := 150, HeightBound := 12, Print := false)
    cols := []; Msel := ZeroMatrix(Q, 10, 0);
    for j := 1 to 40 do
        Mtry := HorizontalJoin(Msel, Matrix(Q, 10, 1, [bm[r,j] : r in [1..10]]));
        if Rank(Mtry) eq Ncols(Mtry) then Msel := Mtry; Append(~cols, j); end if;
        if #cols eq 10 then break; end if;
    end for;
    bmIinv  := ChangeRing(Matrix(Q, [[bm[r,cols[c]] : c in [1..10]] : r in [1..10]])^-1, Kf);
    BmatInv := Matrix(Kf, [Brho[a] : a in [1..10]])^-1;

    a6 := function(u)
        tt := [ &+[ u[r]*Drho[r][j] : r in [1..6] ] : j in [1..6] ];
        for i := 1 to 6 do for j := i+1 to 6 do if tt[i]-tt[j] eq 0 then return false, []; end if; end for; end for;
        for i := 1 to 6 do for j := i+1 to 6 do for k := j+1 to 6 do if tt[i]+tt[j]+tt[k] eq 0 then return false, []; end if; end for; end for; end for;
        if &+tt eq 0 then return false, []; end if;
        gam := all_gamma(Matrix(Kf, 6, 3, [[tt[i], tt[i]^3, 1] : i in [1..6]]));
        if not &and[ gam[c] ne 0 : c in cols ] then return false, []; end if;
        y := Vector(Kf, [gam[c] : c in cols]) * bmIinv;
        w := y * BmatInv;
        jw := 0; for a := 1 to 10 do if w[a] ne 0 then jw := a; break; end if; end for;
        if jw eq 0 then return false, []; end if;
        zq := [];
        for a := 1 to 10 do
            okc, r := IsCoercible(Q, w[a]/w[jw]); if not okc then return false, []; end if; Append(~zq, r);
        end for;
        den := LCM([Denominator(x) : x in zq]); zz := [Z!(den*x) : x in zq]; gg := GCD(zz);
        return true, [Q!(x div gg) : x in zz];
    end function;

    seen := {}; out := [];
    for B := 1 to HeightBound do
        shell := [ [u[i] : i in [1..6]] : u in CartesianPower([-B..B], 6) ];
        shell := [ u : u in shell | Max([Abs(x) : x in u]) eq B and GCD(u) eq 1 ];
        Sort(~shell, func< a, b | &+[Abs(x):x in a] - &+[Abs(x):x in b] >);
        for u in shell do
            okp, zq := a6(u); if not okp then continue; end if;
            ii := ClebschInvariantsConcrete(zq, Brho, Kf, bm : Print := false);
            if ii[5] eq 0 or ii[1] eq 0 then continue; end if;
            key := <ii[2]/ii[1]^2, ii[3]/ii[1]^3, ii[4]/ii[1]^4, ii[5]/ii[1]^5>;
            if key in seen then continue; end if;
            Include(~seen, key);
            dc := DeltaClebsch(ii);
            Append(~out, <IntegerSize(dc), zq, ii, dc>);
            if #out ge Count then break B; end if;
        end for;
    end for;
    Sort(~out, func< a, b | a[1] - b[1] >);
    return out;
end function;

/* ------------------------------------------------------------------------- */
/* Minimize a batch of cubic surfaces in parallel, each under a timeout.      */
/* Returns [<orig_index, minimized_surface, maxcoef>] for those that finished. */
/* ------------------------------------------------------------------------- */

ParallelMinimizeReduce := function(surfaces, dir : Timeout := 60, Poll := 2, Print := false)
    n := #surfaces;
    _ := System(Sprintf("mkdir -p %o", dir));
    _ := System(Sprintf("rm -f %o/pm_*", dir));
    for i := 1 to n do
        PrintFile(Sprintf("%o/pm_surf_%o.m", dir, i),
            Sprintf("R4<x,y,z,w>:=PolynomialRing(IntegerRing(),4);\ngl:=%o;", surfaces[i]) : Overwrite := true);
        job := Sprintf("SetColumns(0);\nload \"%o/pm_surf_%o.m\";\ngl2:=MinimizeReduceCubicSurface(gl);\nSetOutputFile(\"%o/pm_res_%o.txt\":Overwrite:=true);\nprintf \"%%o\", gl2;\nUnsetOutputFile();\nquit;\n",
                       dir, i, dir, i);
        PrintFile(Sprintf("%o/pm_job_%o.m", dir, i), job : Overwrite := true);
        _ := System(Sprintf("(timeout %o magma -b %o/pm_job_%o.m >/dev/null 2>&1; touch %o/pm_done_%o) &",
                            Timeout, dir, i, dir, i));
    end for;
    waited := 0;
    repeat
        _ := System(Sprintf("sleep %o", Poll));
        waited +:= Poll;
        ndone := #[ i : i in [1..n] | System(Sprintf("test -f %o/pm_done_%o", dir, i)) eq 0 ];
    until ndone eq n or waited gt Timeout + 4*Poll + 5;
    R4<x,y,z,w> := PolynomialRing(IntegerRing(), 4);
    results := [];
    for i := 1 to n do
        if System(Sprintf("test -s %o/pm_res_%o.txt", dir, i)) eq 0 then
            try
                gl2 := eval Read(Sprintf("%o/pm_res_%o.txt", dir, i));
                Append(~results, <i, gl2, Max([Abs(c) : c in Coefficients(gl2)])>);
            catch e
                ;
            end try;
        end if;
    end for;
    if Print then printf "    parallel-minimized %o/%o within %os\n", #results, n, Timeout; end if;
    return results;
end function;

/* ------------------------------------------------------------------------- */
/* Best explicit cubics for one twist (K, autos, perms, psi=rho).             */
/* Returns up to PerClass records <surface, maxcoef, z, clebsch, DeltaCl>,     */
/* the smallest-discriminant ones that minimized within the timeout.          */
/* ------------------------------------------------------------------------- */

BestCubicsForTwist := function(Kf, autos, perms, psi, U, dir
        : PerClass := 5, GenCount := 150, MinimizeTop := 16, Timeout := 60, Print := false)
    bm := U`bm;
    Brho := DescendedBasisConcrete(Kf, autos, perms, psi, bm : Print := false);
    Drho := SourceDescendedBasis(Kf, autos, perms, psi, U`phi6 : Print := false);
    moduli := GenerateModuliPoints(Kf, bm, Brho, Drho : Count := GenCount, Print := Print);
    if Print then printf "    %o distinct moduli points; smallest Delta_Cl sizes %o\n",
        #moduli, [moduli[i][1] : i in [1..Min(8,#moduli)]]; end if;

    raws := []; provs := [];
    i := 1;
    while #raws lt MinimizeTop and i le #moduli do
        ii := moduli[i][3]; i +:= 1;
        ok := true; gl := 0;
        try gl := inv_2_gl(ii); catch e ok := false; end try;
        if ok and IsSmoothCubicSurface(gl) then
            Append(~raws, gl); Append(~provs, moduli[i-1]);
        end if;
    end while;

    fin := ParallelMinimizeReduce(raws, dir : Timeout := Timeout, Print := Print);
    /* fin : <index-into-raws, minsurf, maxcoef>; sort the finishers by Delta_Cl */
    pairs := [ <provs[t[1]][1], t> : t in fin ];   /* <DeltaClSize, <idx,surf,maxcoef>> */
    Sort(~pairs, func< a, b | a[1] - b[1] >);
    res := [];
    for k := 1 to Min(PerClass, #pairs) do
        t := pairs[k][2]; prov := provs[t[1]];
        Append(~res, <t[2], t[3], prov[2], prov[3], prov[4]>);  /* surface, maxcoef, z, clebsch, DeltaCl */
    end for;
    return res;
end function;

/* ------------------------------------------------------------------------- */
/* Distinct isomorphisms rho : G -> C (= Gwe6), up to conjugacy by the        */
/* normalizer of C in W(E6).  Conjugate embeddings give isomorphic surfaces,  */
/* so these are the genuinely different twists for the class.  (Any that slip  */
/* through still get removed by the downstream Clebsch deduplication.)         */
/* ------------------------------------------------------------------------- */

TwistEmbeddings := function(U, Gwe6, PK : MaxAut := 2000)
    ok, psi0 := IsIsomorphic(PK, Gwe6);
    error if not ok, "PK is not isomorphic to the chosen W(E6)-subgroup";
    np := Ngens(PK);
    single := [ hom< PK -> Gwe6 | [ PK.i -> psi0(PK.i) : i in [1..np] ] > ];
    bfs := function(idA, gens)   /* all elements of <gens> reachable from idA */
        S := {@ idA @}; fr := [idA];
        while #fr gt 0 do
            nf := [];
            for a in fr do for g in gens do
                b := a*g; if b notin S then Include(~S, b); Append(~nf, b); end if;
            end for; end for;
            fr := nf;
        end while;
        return S;
    end function;
    try
        A := AutomorphismGroup(Gwe6);
        if Order(A) gt MaxAut then return single; end if;     /* too many twists to enumerate */
        N := Normalizer(U`we6, Gwe6);
        ng := Ngens(Gwe6);
        idA := Id(A);
        innerN := [ A ! hom< Gwe6 -> Gwe6 | [ Gwe6.i -> n^-1*Gwe6.i*n : i in [1..ng] ] >
                    : n in Generators(N) ];
        ANset := bfs(idA, innerN);
        allA  := bfs(idA, [A.i : i in [1..Ngens(A)]]);
        /* left-coset reps of AN in A: rho ~ conj_n o rho twists alpha by AN.alpha */
        covered := {@ @}; reps := [];
        for a in allA do
            if a in covered then continue; end if;
            Append(~reps, a);
            for b in ANset do Include(~covered, b*a); end for;
        end for;
        return [ hom< PK -> Gwe6 | [ PK.i -> (al)(psi0(PK.i)) : i in [1..np] ] > : al in reps ];
    catch e
        return single;
    end try;
end function;

/* ------------------------------------------------------------------------- */
/* Best explicit cubics for one W(E6)-conjugacy class, pooling moduli points   */
/* over all distinct embeddings rho.  inv_2_gl reconstructs from the Clebsch    */
/* invariants alone, so the reconstruction/minimization is shared across rho.   */
/* Returns up to PerClass records <surface, maxcoef, z, clebsch, DeltaCl, rho>. */
/* ------------------------------------------------------------------------- */

BestCubicsForClass := function(Kf, autos, perms, PK, Gwe6, U, dir
        : PerClass := 5, GenCount := 150, MinimizeTop := 16, Timeout := 60, Print := false)
    bm := U`bm;
    rhos := TwistEmbeddings(U, Gwe6, PK);
    if Print then printf "    %o distinct embedding(s) rho\n", #rhos; end if;
    allmod := []; seen := {};
    for ri := 1 to #rhos do
        Brho := DescendedBasisConcrete(Kf, autos, perms, rhos[ri], bm : Print := false);
        Drho := SourceDescendedBasis(Kf, autos, perms, rhos[ri], U`phi6 : Print := false);
        m := GenerateModuliPoints(Kf, bm, Brho, Drho : Count := GenCount);
        for x in m do
            ii := x[3];
            key := <ii[2]/ii[1]^2, ii[3]/ii[1]^3, ii[4]/ii[1]^4, ii[5]/ii[1]^5>;
            if key in seen then continue; end if;
            Include(~seen, key);
            Append(~allmod, <x[1], x[2], x[3], x[4], ri>);
        end for;
    end for;
    Sort(~allmod, func< a, b | a[1] - b[1] >);
    if Print then printf "    %o distinct moduli points across embeddings; smallest Delta_Cl %o\n",
        #allmod, [allmod[i][1] : i in [1..Min(6,#allmod)]]; end if;

    raws := []; provs := [];
    i := 1;
    while #raws lt MinimizeTop and i le #allmod do
        ii := allmod[i][3]; i +:= 1;
        ok := true; gl := 0;
        try gl := inv_2_gl(ii); catch e ok := false; end try;
        if ok and IsSmoothCubicSurface(gl) then Append(~raws, gl); Append(~provs, allmod[i-1]); end if;
    end while;

    fin := ParallelMinimizeReduce(raws, dir : Timeout := Timeout, Print := Print);
    pairs := [ <provs[t[1]][1], t> : t in fin ];
    Sort(~pairs, func< a, b | a[1] - b[1] >);
    res := [];
    for k := 1 to Min(PerClass, #pairs) do
        t := pairs[k][2]; prov := provs[t[1]];
        Append(~res, <t[2], t[3], prov[2], prov[3], prov[4], prov[5]>);
    end for;
    return res, #rhos;
end function;

/* ------------------------------------------------------------------------- */
/* Top-level driver: a polynomial f -> the database rows for every W(E6)-      */
/* conjugacy class of subgroups isomorphic to Gal(f).  Each row is             */
/*   <class_index, surface, maxcoef, z, clebsch, DeltaCl, rho_index>.          */
/* ------------------------------------------------------------------------- */

GenerateDatabase := function(f, U, dir : PerClass := 5, GenCount := 150,
        MinimizeTop := 16, Timeout := 60, MaxClasses := 10000, Print := true)
    fz := IntegralMonicPolynomial(f);
    P, rts, GalData := GaloisGroup(fz);
    gfull := GaloisSubgroup(GalData, sub<P|>);
    Kf := NumberField(gfull); dK := Degree(gfull);
    autos := Automorphisms(Kf);
    rtsK := [r[1] : r in Roots(gfull, Kf)];
    Sd := SymmetricGroup(dK);
    perms := [ Sd ! [ Index(rtsK, au(rtsK[i])) : i in [1..dK] ] : au in autos ];
    PK := sub<Sd|perms>;
    cl := SubgroupClasses(U`we6);
    classes := [ c`subgroup : c in cl |
                 Order(c`subgroup) eq Order(P) and IsIsomorphic(c`subgroup, P) ];
    if Print then printf "Gal(f) order %o; %o matching W(E6) conjugacy class(es)\n", Order(P), #classes; end if;
    rows := [];
    for ci := 1 to Min(MaxClasses, #classes) do
        if Print then printf "  class %o/%o (order %o):\n", ci, #classes, Order(classes[ci]); end if;
        res, nrho := BestCubicsForClass(Kf, autos, perms, PK, classes[ci], U, dir
            : PerClass := PerClass, GenCount := GenCount, MinimizeTop := MinimizeTop,
              Timeout := Timeout, Print := Print);
        for r in res do
            Append(~rows, <ci, r[1], r[2], r[3], r[4], r[5], r[6]>);
        end for;
        if Print then printf "    -> %o cubic(s)\n", #res; end if;
    end for;
    return rows;
end function;

/* ------------------------------------------------------------------------- */
/* Canonical W(E6) subgroup labels (WE6subgroups.txt).                        */
/*                                                                            */
/* The file gives, per line "label:gens:d:t", a conjugacy class of subgroups  */
/* of W(E6) = 27T1161 by generators in S27 (acting on the 27 lines), with the  */
/* degree d and T-number t of the minimal transitive group isomorphic to it.   */
/* IsIsomorphic(27T1161, we6) transports each into the 40-dimensional gamma     */
/* representation used by the descent (this transport is the geometric one --   */
/* checked: it agrees with the explicit 27-line action of the generators on all */
/* 350 classes, so the W(E6) outer automorphism does not affect the labels).    */
/* ------------------------------------------------------------------------- */

LoadWE6Subgroups := function(U, filename)
    we6 := U`we6;
    S27 := SymmetricGroup(27);
    S := [ Split(r, ":") : r in Split(Read(filename)) ];
    G := [ sub<S27 | eval(r[2])> : r in S ];
    full := [ i : i in [1..#S] | StringToInteger(S[i][3]) eq 27 and StringToInteger(S[i][4]) eq 1161 ];
    error if not (#full ge 1), "Could not find the full W(E6) (27T1161) line in the file";
    ok, psi := IsIsomorphic(G[full[1]], we6);
    error if not ok, "27T1161 is not isomorphic to the 40-dimensional W(E6)";
    subs := [];
    for i in [1..#S] do
        Gwe6 := sub<we6 | [ psi(g) : g in Generators(G[i]) ]>;
        Append(~subs, < S[i][1], Gwe6, StringToInteger(S[i][3]), StringToInteger(S[i][4]) >);
    end for;
    return subs;
end function;

/* ------------------------------------------------------------------------- */
/* Realize the Galois group of f on the 27 lines as cubic surfaces, labelled   */
/* by the canonical W(E6) subgroup class.  Returns a list of <label, cubic     */
/* form> tuples, up to n per label, the labels ranging over the W(E6) subgroup  */
/* classes isomorphic to Gal(f).                                               */
/* ------------------------------------------------------------------------- */

RealizeCubicSurfaces := function(f, n, U, subs, dir
        : GenCount := 150, MinimizeTop := 16, Timeout := 60, Print := true)
    fz := IntegralMonicPolynomial(f);
    P, rts, GalData := GaloisGroup(fz);
    t, d := TransitiveGroupIdentification(P);
    T := [ i : i in [1..#subs] | subs[i][3] eq d and subs[i][4] eq t ];
    if Print then printf "Gal(f) = %oT%o; %o matching W(E6) subgroup label(s)\n", d, t, #T; end if;

    gfull := GaloisSubgroup(GalData, sub<P|>);
    Kf := NumberField(gfull); dK := Degree(gfull);
    autos := Automorphisms(Kf);
    rtsK := [r[1] : r in Roots(gfull, Kf)];
    Sd := SymmetricGroup(dK);
    perms := [ Sd ! [ Index(rtsK, au(rtsK[i])) : i in [1..dK] ] : au in autos ];
    PK := sub<Sd|perms>;

    result := [];
    for idx in T do
        label := subs[idx][1]; Gwe6 := subs[idx][2];
        if Print then printf "  label %o:\n", label; end if;
        res := BestCubicsForClass(Kf, autos, perms, PK, Gwe6, U, dir
            : PerClass := n, GenCount := GenCount, MinimizeTop := MinimizeTop,
              Timeout := Timeout, Print := Print);
        for r in res do Append(~result, < label, r[1] >); end for;
        if Print then printf "    -> %o cubic(s)\n", #res; end if;
    end for;
    return result;
end function;
