#!/usr/bin/env python3
"""merge_databases.py -- fold the re-emission and new-class sweeps into the
database files, then regenerate the tsv/html.

Inputs (paths given as arguments):
  --reemit DIR   chunk outputs o* from reemit_nosplit.m
                 (label:src:orb:verdict:digits:cubic, plus "# KEEP label ..." rows)
  --seed2 DIR    chunk outputs o* from seed_nosplit2.m (same format,
                 plus "# WALL/#NOPOLY" rows)
  --trivial FILE optional "label:src:cubic" row for the split (trivial) class,
                 appended to database_seed.txt if the label is absent

Updates database_seed_nosplit.txt in place (4-field rows label:src:orb:cubic,
smallest model wins), writes nosplit_verdicts.txt (label:verdict:digits), and
prints a summary.  Run class_info_all.m and make_seed_html.py afterwards.
"""
import argparse, glob, os, re, sys

def coefdigits(cubic):
    return max(len(m) for m in re.findall(r"\d+", cubic))

def read_chunks(d):
    rows, keeps = {}, []
    for fn in sorted(glob.glob(os.path.join(d, "o*"))):
        for ln in open(fn):
            ln = ln.rstrip("\n")
            if not ln:
                continue
            if ln.startswith("#"):
                keeps.append(ln)
                continue
            p = ln.split(":")
            if len(p) < 6:
                continue
            label, src, orb, verdict, digits = p[0], p[1], p[2], p[3], p[4]
            cubic = ":".join(p[5:])
            if verdict not in ("CERTIFIED", "CONSISTENT"):
                print(f"DROPPED (verdict {verdict!r}): {label}")
                continue
            rows[label] = dict(src=src, orb=orb, verdict=verdict,
                               digits=int(digits), cubic=cubic)
    return rows, keeps

ap = argparse.ArgumentParser()
ap.add_argument("--reemit")
ap.add_argument("--seed2")
ap.add_argument("--trivial")
args = ap.parse_args()

old = {}
order = []
for ln in open("database_seed_nosplit.txt"):
    ln = ln.rstrip("\n")
    if not ln or ln.startswith("#"):
        continue
    p = ln.split(":")
    label = p[0]
    old[label] = dict(src=p[1], orb=p[2], cubic=":".join(p[3:]))
    order.append(label)

verdicts = {}
improved = same = added = 0

if args.reemit:
    rows, keeps = read_chunks(args.reemit)
    for label, r in rows.items():
        if label not in old:
            print(f"WARNING: reemit produced unknown label {label}")
            continue
        oldd = coefdigits(old[label]["cubic"])
        verdicts[label] = (r["verdict"], r["digits"])
        if r["digits"] < oldd:
            old[label] = dict(src=r["src"], orb=r["orb"], cubic=r["cubic"])
            improved += 1
        else:
            same += 1
    print(f"reemit: {improved} improved, {same} kept, {len(keeps)} comment rows")

if args.seed2:
    rows, keeps = read_chunks(args.seed2)
    for label, r in rows.items():
        if label in old:
            oldd = coefdigits(old[label]["cubic"])
            if r["digits"] < oldd:
                old[label] = dict(src=r["src"], orb=r["orb"], cubic=r["cubic"])
                improved += 1
        else:
            old[label] = dict(src=r["src"], orb=r["orb"], cubic=r["cubic"])
            order.append(label)
            added += 1
        verdicts[label] = (r["verdict"], r["digits"])
    walls = [k for k in keeps if "WALL" in k or "NOPOLY" in k]
    print(f"seed2: {added} new classes, {len(walls)} walls")
    for w in walls:
        print("   ", w)

with open("database_seed_nosplit.txt", "w") as f:
    for label in order:
        r = old[label]
        f.write(f"{label}:{r['src']}:{r['orb']}:{r['cubic']}\n")

with open("nosplit_verdicts.txt", "w") as f:
    f.write("# label:verdict:max_coeff_digits  (LinesGaloisCertificate on the stored model)\n")
    for label in sorted(verdicts):
        v = verdicts[label]
        f.write(f"{label}:{v[0]}:{v[1]}\n")

if args.trivial and os.path.exists(args.trivial):
    row = open(args.trivial).read().strip()
    label = row.split(":")[0]
    seed = open("database_seed.txt").read()
    if label not in seed:
        with open("database_seed.txt", "a") as f:
            f.write(row + "\n")
        print(f"trivial class appended: {label}")

digits = sorted(coefdigits(r["cubic"]) for r in old.values())
print(f"nosplit rows: {len(order)}; digit quartiles "
      f"{digits[len(digits)//4]}/{digits[len(digits)//2]}/{digits[3*len(digits)//4]}, max {digits[-1]}")
