/* seed_database.m -- initial seeding of the cubic-surface database.

   Runs every polynomial in WE6fields.txt (rows "d:t:[f0,...,fd]") through
   RealizeCubicSurfaces(f, n) and appends, per realized class, lines

       label:source_coeffs:cubic_form

   to database_seed.txt.  Resumable: seed_progress.txt records the last completed
   row, so re-running continues where it stopped.  Run from the repo root:

       magma -b seed_database.m

   Optional: set N (per-label count), GENCOUNT, MINTOP, TMO via the environment.
*/
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "database_pipeline.m";
Q := RationalField(); Z := IntegerRing();
Pq<t> := PolynomialRing(Q);

geti := function(name, def)
    s := GetEnvironmentValue(name);
    return s eq "" select def else StringToInteger(s);
end function;
N        := geti("N", 3);
GENCOUNT := geti("GENCOUNT", 50);
MINTOP   := geti("MINTOP", 10);
TMO      := geti("TMO", 30);

dir      := "/tmp/claude-1000/-home-claude-CubicSurfaces/dfcfed9a-a863-46fb-80be-f21032adb5c7/scratchpad/race";
outfile  := "database_seed.txt";
progfile := "/tmp/claude-1000/-home-claude-CubicSurfaces/dfcfed9a-a863-46fb-80be-f21032adb5c7/scratchpad/seed_progress.txt";

U    := BuildUniversalCobleData(: Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
rows := [ Split(r, ":") : r in Split(Read("WE6fields.txt")) ];
printf "loaded %o subgroup labels and %o polynomials (N=%o, GenCount=%o, MinTop=%o, Timeout=%o)\n",
       #subs, #rows, N, GENCOUNT, MINTOP, TMO;

start := 0;
try start := StringToInteger(Read(progfile)); catch e start := 0; end try;
if start gt 0 then printf "resuming after row %o\n", start; end if;

tStart := Cputime();
totcubics := 0;
for i := start+1 to #rows do
    r := rows[i];
    coeffs := eval r[3];
    f := &+[ coeffs[k+1]*t^k : k in [0..#coeffs-1] ];
    t0 := Cputime();
    nc := 0;
    try
        res := RealizeCubicSurfaces(f, N, U, subs, dir
            : GenCount := GENCOUNT, MinimizeTop := MINTOP, Timeout := TMO, Print := false);
        for tup in res do
            PrintFile(outfile, Sprintf("%o:%o:%o", tup[1], r[3], tup[2]));
            nc +:= 1;
        end for;
    catch e
        PrintFile(outfile, Sprintf("# ERROR row %o coeffs %o: %o", i, r[3], e`Object));
    end try;
    totcubics +:= nc;
    PrintFile(progfile, Sprintf("%o", i) : Overwrite := true);
    printf "[%o/%o] %oT%o %o -> %o cubic(s)  %.1os  (total %o, %.0os elapsed)\n",
           i, #rows, r[1], r[2], r[3], nc, Cputime(t0), totcubics, Cputime(tStart);
end for;
printf "DONE: %o cubics from %o polynomials in %.0os\n", totcubics, #rows, Cputime(tStart);
quit;
