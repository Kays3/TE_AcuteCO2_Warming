"""A1.9: for Unknown families with FDR < 0.05, the share of their genomic copies (TE GTF records) that overlap
a RefSeq exon by at least 1 bp. Usage: unknown_exon_share.py FAMILIES.txt TE.gtf REFSEQ.gtf OUT.tsv"""
import sys, collections, bisect
fams = set(l.strip() for l in open(sys.argv[1]) if l.strip())
ex = collections.defaultdict(list)
for line in open(sys.argv[3]):
    if line.startswith("#"): continue
    p = line.split("\t", 8)
    if len(p) > 2 and p[2] == "exon": ex[p[0]].append((int(p[3]), int(p[4])))
merged = {}
for c, iv in ex.items():
    iv.sort(); m = []
    for s, e in iv:
        if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    merged[c] = ([s for s, e in m], [e for s, e in m])
n = collections.Counter(); hit = collections.Counter()
for line in open(sys.argv[2]):
    p = line.split("\t")
    if len(p) < 9: continue
    g = p[8].split('gene_id "')[1].split('"')[0]
    if g not in fams: continue
    c, s, e = p[0], int(p[3]), int(p[4]); n[g] += 1
    if c in merged:
        st, en = merged[c]; i = bisect.bisect_right(st, e) - 1
        if i >= 0 and en[i] >= s: hit[g] += 1
with open(sys.argv[4], "w") as o:
    o.write("family\tn_copies\tn_exon_overlap\tshare\n")
    for g in sorted(fams): o.write(f"{g}\t{n[g]}\t{hit[g]}\t{hit[g]/n[g] if n[g] else float('nan'):.3f}\n")
