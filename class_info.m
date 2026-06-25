/* class_info.m -- emit class_info.txt with one line "label:d:t:order:groupname"
   for every W(E6) class that appears in database_seed.txt.  Used by
   make_seed_html.py to label the browser.  Run from the repo root:
       magma -b class_info.m
*/
SetColumns(0);
S27 := SymmetricGroup(27);
info := AssociativeArray();
for r in [ Split(x, ":") : x in Split(Read("WE6subgroups.txt")) ] do
    info[r[1]] := <r[2], StringToInteger(r[3]), StringToInteger(r[4])>;
end for;

labels := {};
for ln in Split(Read("database_seed.txt")) do
    if #ln gt 0 and ln[1] ne "#" then Include(~labels, Split(ln, ":")[1]); end if;
end for;

namecache := AssociativeArray();
lines := [];
for lbl in labels do
    rec := info[lbl];
    key := <rec[2], rec[3]>;
    if not IsDefined(namecache, key) then
        G := sub<S27 | eval(rec[1])>;
        namecache[key] := <#G, GroupName(G)>;
    end if;
    nm := namecache[key];
    Append(~lines, Sprintf("%o:%o:%o:%o:%o", lbl, rec[2], rec[3], nm[1], nm[2]));
end for;
Sort(~lines);
PrintFile("class_info.txt", Join(lines, "\n") : Overwrite := true);
printf "wrote class_info.txt: %o classes\n", #lines;
quit;
