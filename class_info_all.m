/* class_info_all.m -- emit  label:order:groupname:orbit_sizes  for every W(E6)
   class present in either database_seed.txt or database_seed_nosplit.txt, using
   the S27 generators in WE6subgroups.txt.  Feeds the combined HTML/TSV. */
SetColumns(0);
S27 := SymmetricGroup(27);
labels := {};
for fn in ["database_seed.txt", "database_seed_nosplit.txt"] do
    try for l in Split(Read(fn)) do
        if #l ge 5 and l[1..5] eq "51840" then Include(~labels, Split(l, ":")[1]); end if;
    end for; catch e; end try;
end for;
gensOf := AssociativeArray();
for l in Split(Read("WE6subgroups.txt")) do
    p := Split(l, ":"); if #p ge 2 then gensOf[p[1]] := p[2]; end if;
end for;
OUT := "class_info_all.txt"; System("rm -f " cat OUT);
n := 0;
for lab in Sort(SetToSequence(labels)) do
    if not IsDefined(gensOf, lab) then continue; end if;
    H := eval("sub<S27 | " cat gensOf[lab] cat ">");
    orb := Sort([ #o : o in Orbits(H) ]);
    name := "?"; try name := GroupName(H); catch e; end try;
    PrintFile(OUT, Sprintf("%o:%o:%o:%o", lab, Order(H), name, orb));
    n +:= 1;
end for;
printf "wrote class_info_all.txt for %o labels\n", n;
quit;
