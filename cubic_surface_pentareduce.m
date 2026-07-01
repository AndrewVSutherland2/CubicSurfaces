/* cubic_surface_pentareduce.m -- small models from Clebsch invariants via the
   (K5, lambda) pentahedral coordinates, without factoring model discriminants.

   The reconstruction penta_pol_2_gl builds  Tr_{A/Q}(lambda * l(x)^3)  on the
   trace-zero hyperplane of A = Q[T]/(pol5) from the raw power basis, with all
   of the gauge freedom  y -> c*y  (lambda -> lambda*c^3)  unused; its models
   carry ~3x the quintic's height and their discriminants are unfactorable.
   In the (K5, lambda) coordinates the only integers that ever need factoring
   are E = I40 (N(lambda) = E^2) and Delta_Cl -- heights of the *moduli point*,
   not the model.  This module factors those with bounded effort, hands the
   quintic + prime list to pentareduce.gp (polredbest + ideal cube-stripping +
   weighted LLL of the trace-zero lattice), certifies the returned model in
   Magma (same weighted Clebsch point, nonzero discriminant), and polishes with
   MinimizeCubicSurface at the now-known primes followed by ReduceCubicSurface.

   Load after cubic_surface_resolvent_twist.m.  Needs `gp` on the PATH.

   Main entry points:
     PentahedralReduce(inv : ...)        -> ok, cubic form, metadata record
     PentaReduceCubicSurface(F : ...)    -> ok, cubic form, metadata record
*/

PentaReduceRF := recformat<
    surface,            // the reduced cubic form (over Q[x1..x4])
    digits,             // max |coefficient| digit count
    primes,             // known prime support used (SetEnum)
    leftover,           // unfactored composite parts of E / Delta_Cl (SetEnum)
    verified,           // invariant-proportionality certificate passed
    diagnostics         // gp diagnostics string
>;

/* ---- bounded-effort integer factorization ------------------------------- */
/* Returns <primes, leftover-composites>.  Trial division to TrialLimit, then
   perfect-power reduction, ECM rounds (ECMRounds = [<B1, curves>, ...]), and
   a full Factorization (MPQS) only for composites of at most MPQSDigitCap
   digits.  Probable primes are accepted as primes (fine here: they feed
   polredbest/nfinit prime lists and MinimizeCubicSurface; a pseudoprime would
   at worst weaken the reduction, not break correctness).                     */
BoundedFactorInteger := function(n : TrialLimit := 10^6,
                                     ECMRounds := [<5*10^4, 40>, <25*10^4, 30>, <10^6, 15>],
                                     MPQSDigitCap := 68)
    n := Abs(n);
    if n le 1 then return {Integers()|}, {Integers()|}; end if;
    fs, comps := TrialDivision(n, TrialLimit);
    primes := { f[1] : f in fs };
    todo := [ c : c in comps ];
    leftover := {Integers()|};
    while #todo gt 0 do
        c := todo[#todo]; Prune(~todo);
        if c eq 1 then continue; end if;
        if IsProbablePrime(c) then Include(~primes, c); continue; end if;
        ispow, root := IsPower(c);
        if ispow then Append(~todo, root); Append(~todo, c div root); continue; end if;
        split := false;
        for round in ECMRounds do
            for curve in [1..round[2]] do
                g := ECM(c, round[1]);
                if g ne 0 and g ne 1 and g ne c then
                    Append(~todo, g); Append(~todo, c div g);
                    split := true; break;
                end if;
            end for;
            if split then break; end if;
        end for;
        if not split and #Sprint(c) le MPQSDigitCap then
            try
                for f in Factorization(c : Proof := false) do
                    Include(~primes, f[1]);
                end for;
                split := true;
            catch e ; end try;
        end if;
        if not split then Include(~leftover, c); end if;
    end while;
    return primes, leftover;
end function;

/* ---- the intrinsic discriminant from the Clebsch invariants -------------- */
DeltaClOfInvariants := function(inv)
    A := inv[1]; B := inv[2]; C := inv[3]; D := inv[4];
    return A^4 - 128*A^2*B + 4096*B^2 - 2048*A*C - 16384*D;
end function;

/* ---- weighted-projective equality of two invariant vectors --------------- */
/* I_k has weight 8k, so vectors agree on the same moduli point iff the
   zero-patterns match and (Ik'/Ik)^j = (Ij'/Ij)^k for all nonzero pairs.     */
SameCleschPoint := function(inv1, inv2)
    for k in [1..5] do
        if (inv1[k] eq 0) ne (inv2[k] eq 0) then return false; end if;
    end for;
    nz := [ k : k in [1..5] | inv1[k] ne 0 ];
    for a in [1..#nz-1] do
        j := nz[a]; k := nz[a+1];
        if (inv2[k]/inv1[k])^j ne (inv2[j]/inv1[j])^k then return false; end if;
    end for;
    return true;
end function;

/* ---- core: Clebsch invariants -> reduced model --------------------------- */
PentahedralReduce := function(inv : TmpDir := ".", GPProg := "",
                                    TrialLimit := 10^6,
                                    ECMRounds := [<5*10^4, 40>, <25*10^4, 30>, <10^6, 15>],
                                    MPQSDigitCap := 68, Polish := true,
                                    PolishDigitCap := 400, Print := false)
    err := rec< PentaReduceRF | verified := false >;
    if inv[5] eq 0 then return false, err, "E = 0: outside the pentahedral locus"; end if;
    A := Rationals()!inv[1]; B := Rationals()!inv[2]; C := Rationals()!inv[3];
    D := Rationals()!inv[4]; E := Rationals()!inv[5];
    DCl := DeltaClOfInvariants([A,B,C,D,E]);
    if DCl eq 0 then return false, err, "Delta_Cl = 0: singular surface"; end if;

    /* pentahedral quintic, made integral-monic by T -> T/m (lambda -> m*lambda) */
    s := [B, D, (C^2 - A*E)/4, C*E, E^2];
    m := LCM([ Denominator(x) : x in s ]);
    sc := [ Integers() | m^k * s[k] : k in [1..5] ];
    PC := [ -sc[5], sc[4], -sc[3], sc[2], -sc[1], 1 ];

    /* known prime support: E, Delta_Cl, the scaling m */
    primes := {Integers()|}; leftover := {Integers()|};
    for t in [ Numerator(E), Denominator(E), Numerator(DCl), Denominator(DCl), m ] do
        p1, l1 := BoundedFactorInteger(t : TrialLimit := TrialLimit,
            ECMRounds := ECMRounds, MPQSDigitCap := MPQSDigitCap);
        primes join:= p1; leftover join:= l1;
    end for;
    if Print then
        printf "pentareduce: %o known primes, %o unfactored composite(s)\n",
            #primes, #leftover;
    end if;

    /* hand off to pentareduce.gp */
    gpprog := GPProg eq "" select GetCurrentDirectory() cat "/pentareduce.gp" else GPProg;
    inlines := Sprintf("PC = %o;\nSP = %o;\n",
        PC, Sort(SetToSequence(primes)));
    PrintFile(TmpDir cat "/pentared_in.gp", inlines : Overwrite := true);
    System("cd " cat TmpDir cat " && rm -f pentared_out.txt && gp -q " cat
           gpprog cat " > pentared_gp.log 2>&1");
    outstr := "";
    try outstr := Read(TmpDir cat "/pentared_out.txt"); catch e ; end try;
    if outstr eq "" then return false, err, "gp produced no output file"; end if;
    out := Split(outstr, "\n");
    if #out lt 2 or #out[1] lt 2 or out[1][1..2] ne "OK" then
        return false, err, "gp: " cat (#out ge 1 select out[1] else "empty output");
    end if;
    cfz := [ StringToInteger(x) : x in Split(out[2], " ") | #x gt 0 ];
    if #cfz ne 20 then return false, err, "expected 20 coefficients from gp"; end if;

    /* rebuild the form; monomial order x_i x_j x_k, i <= j <= k */
    R4 := PolynomialRing(Rationals(), 4);
    AssignNames(~R4, ["x", "y", "z", "w"]);
    idx := 0; F := R4!0;
    for i in [1..4] do for j in [i..4] do for k in [j..4] do
        idx +:= 1; F +:= cfz[idx] * R4.i * R4.j * R4.k;
    end for; end for; end for;

    /* certify: same weighted Clebsch point, nonzero discriminant */
    invN, disc := ClebschSalmonInvariants(F);
    if disc eq 0 or not SameCleschPoint([A,B,C,D,E], invN) then
        return false, err, "verification failed: model does not match the moduli point";
    end if;

    /* polish: minimize at the known primes (cheap local linear algebra --
       no factoring), then LLL-reduce; keep whatever is smallest */
    digitsOf := func< f | Max([ #Sprint(Abs(Numerator(c))) : c in Coefficients(f) ]) >;
    best := F;
    if Polish and digitsOf(F) le PolishDigitCap then
        Fm := F;
        for p in Sort(SetToSequence(primes)) do
            try Fm := MinimizeCubicSurface(Fm, p); catch e ; end try;
        end for;
        try Fm := ReduceCubicSurface(PolynomialRing(Integers(),4)!Fm); catch e ; end try;
        Fm := R4!Fm;
        invM, discM := ClebschSalmonInvariants(Fm);
        if discM ne 0 and SameCleschPoint([A,B,C,D,E], invM)
           and digitsOf(Fm) lt digitsOf(best) then best := Fm; end if;
    end if;

    res := rec< PentaReduceRF |
        surface := best, digits := digitsOf(best), primes := primes,
        leftover := leftover, verified := true,
        diagnostics := #out ge 3 select out[3] else "" >;
    return true, res, _;
end function;

/* ---- convenience wrapper: reduce an explicit cubic surface --------------- */
/* Computes the invariants from the given model, so only sensible when the
   input model is small enough for ClebschSalmonInvariants (its invariants
   inherit the model's scaling; the weighted-projective gauge is recovered
   through the prime list, so a badly scaled huge input gives a weak prime
   list -- prefer PentahedralReduce on reconstruction-normalized invariants). */
PentaReduceCubicSurface := function(F : TmpDir := ".", GPProg := "",
                                        TrialLimit := 10^6, Polish := true,
                                        Print := false)
    inv, disc := ClebschSalmonInvariants(F);
    if disc eq 0 then
        return false, rec< PentaReduceRF | verified := false >, "singular surface";
    end if;
    return PentahedralReduce(inv : TmpDir := TmpDir, GPProg := GPProg,
        TrialLimit := TrialLimit, Polish := Polish, Print := Print);
end function;
