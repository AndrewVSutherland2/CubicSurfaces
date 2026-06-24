/* point_search_helpers.m

   Exploratory routines used while investigating the rational-point search on the
   descended twist -- in particular for the harder, non-cyclic groups whose
   accessible points tend to be degenerate.  These are research helpers, not part
   of the certified pipeline in cubic_surface_resolvent_twist.m; they are kept
   here for the record.  Load from the repository root:

      load "experiments/point_search_helpers.m";

   Provides:
     BuildTwist(U, G, f)            -- descent data (K, autos, perms, psi) for a pair
     BFofBrho / BrhoOfBF            -- convert a fixed-space basis to/from a flat matrix
     MinkowskiReduce(Brho, K)       -- Minkowski/T2-reduce the fixed lattice
     MakeValidator(Brho, K, bm)     -- smooth + proper-pentahedron + separable test
     FiberSearch(...)               -- search the fourfold by fixing four coordinates
     BoundedSearch(kub, B, validate)-- bounded integer-point search (the EJ style)
*/

load "cubic_surface_resolvent_twist.m";

// Build descent data for a (subgroup, polynomial) pair.
BuildTwist := function(U, Gwe6, f)
   fz := IntegralMonicPolynomial(f);
   P, roots, GalData := GaloisGroup(fz);
   gfull := GaloisSubgroup(GalData, sub<P|>);
   K := NumberField(gfull); dK := Degree(gfull);
   autos := Automorphisms(K); rtsK := [r[1] : r in Roots(gfull, K)];
   Sd := SymmetricGroup(dK);
   perms := [ Sd ! [ Index(rtsK, au(rtsK[i])) : i in [1..dK] ] : au in autos ];
   PK := sub<Sd | perms>; ok, psi := IsIsomorphic(PK, Gwe6);
   error if not ok, "not isomorphic";
   return K, dK, autos, perms, psi, gfull;
end function;

// Raw fixed-space basis (BF, rows = vec(U)) from a Brho.
BFofBrho := function(Brho, K)
   d := Degree(K);
   BF := ZeroMatrix(Q, 10, 10*d);
   for r in [1..10] do for j in [1..10] do
      cf := Eltseq(Brho[r][j]);
      for b in [1..d] do BF[r][(b-1)*10+j] := cf[b]; end for;
   end for; end for;
   return BF;
end function;

BrhoOfBF := function(BF, K)
   d := Degree(K);
   return [ [ K![ BF[r][(b-1)*10+j] : b in [1..d] ] : j in [1..10] ] : r in [1..10] ];
end function;

// Minkowski reduction of the fixed lattice using the T2 (Minkowski) Gram of K.
// For totally real K the trace form is exact; otherwise we use a high-precision
// real approximation, which still yields an exact integral transformation.
MinkowskiReduce := function(Brho, K)
   d := Degree(K); a := K.1;
   BF := BFofBrho(Brho, K);
   if IsTotallyReal(K) then
      G := Matrix(Q, d, d, [ Trace(a^(p+q-2)) : p,q in [1..d] ]);
      gden := LCM([Denominator(G[i,j]) : i,j in [1..d]]);
      Gi := Matrix(Z, d, d, [ Z!(gden*G[i,j]) : i,j in [1..d] ]);
   else
      prec := 200;
      cj := func< x | Conjugates(x : Precision := prec) >;
      RR := RealField(prec);
      pb := [a^(i-1) : i in [1..d]];
      C := [cj(b) : b in pb];
      // T2 inner product <x,y> = sum_k Re(sigma_k(x) conj(sigma_k(y)))
      Gr := Matrix(RR, d, d, [ &+[ Real(C[p][k]*ComplexConjugate(C[q][k])) : k in [1..d] ] : p,q in [1..d] ]);
      Gi := Matrix(Z, d, d, [ Round(10^20*Gr[i,j]) : i,j in [1..d] ]);
   end if;
   N := 10*d;
   BigG := ZeroMatrix(Z, N, N);
   for j in [1..10] do for b in [1..d] do for bp in [1..d] do
      BigG[(b-1)*10+j, (bp-1)*10+j] := Gi[b,bp];
   end for; end for; end for;
   bfden := MatrixDenominator(BF);
   BFint := ChangeRing(bfden*BF, Z);
   Lat := Lattice(BFint, BigG);
   BFr := ChangeRing(BasisMatrix(LLL(Lat)), Q);
   return BrhoOfBF(BFr, K);
end function;

// Validation: smooth, proper pentahedron, separable.
MakeValidator := function(Brho, K, bm)
   return function(s)
      ii := ClebschInvariantsConcrete(s, Brho, K, bm : Print := false);
      if ii[5] eq 0 then return false; end if;
      pol5 := ClebschInvariantsToPentahedralPolynomial(ii);
      if Discriminant(pol5) eq 0 then return false; end if;
      return IsSmoothCubicSurface(penta_pol_2_gl(pol5));
   end function;
end function;

// dim-4 fibre search: fix 4 coordinate-ratios (over coord "nrm"), search a box,
// substitute, and solve the residual 0-dimensional scheme in P^5.
FiberSearch := function(kub, Brho, K, bm, B, charts, validate : Print := false)
   R := Parent(kub[1]);
   for ch in charts do
      S := ch[1]; nrm := ch[2];      // S = 4 coords to fix as ratios of coord nrm
      rest := [j : j in [1..10] | j notin S];   // 6 remaining coords (includes nrm)
      P5 := ProjectiveSpace(Q, 5);
      for tup in CartesianPower([-B..B], 4) do
         c := [tup[1],tup[2],tup[3],tup[4]];
         if GCD([Z|x : x in c] cat [1]) ne 1 then continue; end if;
         // substitution: z_{S[t]} = c[t]*z_nrm ; z_{rest} -> P5 vars
         img := [R| 0 : j in [1..10]];
         posn := Index(rest, nrm);
         for t in [1..4] do img[S[t]] := c[t]*P5.posn; end for;
         for t in [1..6] do img[rest[t]] := P5.t; end for;
         hm := hom< R -> CoordinateRing(P5) | img >;
         eqs := [hm(F) : F in kub];
         Y := Scheme(P5, eqs);
         if Dimension(Y) ne 0 then continue; end if;
         for p in RationalPoints(Y) do
            pc := Eltseq(p);
            sol := [Q| 0 : j in [1..10]];
            for t in [1..6] do sol[rest[t]] := pc[t]; end for;
            for t in [1..4] do sol[S[t]] := c[t]*pc[posn]; end for;
            if validate(sol) then return true, sol; end if;
         end for;
      end for;
   end for;
   return false, [];
end function;

// Bounded integer point search (the EJ style): enumerate primitive projective
// integer vectors of height <= B and short-circuit the 30 cubics.
BoundedSearch := function(kub, B, validate)
   f1 := kub[1]; rest := [kub[i] : i in [2..#kub]];
   vals := [-B..B];
   for v in CartesianPower(vals, 10) do
      sol := [v[i] : i in [1..10]];
      if Evaluate(f1, sol) ne 0 then continue; end if;
      // primitive + canonical sign
      nz := [i : i in [1..10] | sol[i] ne 0];
      if #nz eq 0 then continue; end if;
      if sol[nz[1]] lt 0 then continue; end if;
      if GCD([Z|x : x in sol]) ne 1 then continue; end if;
      if forall{F : F in rest | Evaluate(F, sol) eq 0} and validate(sol) then
         return true, sol;
      end if;
   end for;
   return false, [];
end function;
