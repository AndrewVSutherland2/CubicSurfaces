/* class_info_nosplit.m -- emit label:d:t:order:groupname for every class in
   database_seed_nosplit.txt (and nosplit_unrealized.txt), for the HTML browser. */
SetColumns(0);
load "cubic_surface_resolvent_twist.m";
load "database_pipeline.m";
U    := BuildUniversalCobleData(: Print := false);
subs := LoadWE6Subgroups(U, "WE6subgroups.txt");
byLabel := AssociativeArray(); for s in subs do byLabel[s[1]] := s; end for;

labels := [];
for fn in ["database_seed_nosplit.txt", "nosplit_unrealized.txt"] do
    for l in Split(Read(fn)) do
        p := Split(l, ":");
        if #p ge 1 and #p[1] ge 6 and p[1][1..5] eq "51840" then Append(~labels, p[1]); end if;
    end for;
end for;

OUT := "class_info_nosplit.txt";
System("rm -f " cat OUT);
for lab in labels do
    if not IsDefined(byLabel, lab) then continue; end if;
    s := byLabel[lab]; G := s[2];
    name := "?";
    try name := GroupName(G); catch e try name := GroupName(G : TeX := false); catch ee ; end try; end try;
    PrintFile(OUT, Sprintf("%o:%o:%o:%o:%o", lab, s[3], s[4], Order(G), name));
end for;
printf "wrote class_info_nosplit.txt for %o labels\n", #labels;
quit;
