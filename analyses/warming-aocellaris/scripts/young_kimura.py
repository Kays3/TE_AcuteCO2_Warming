"""5.3.3 / A1.12: median CpG-adjusted Kimura distance per family from the RepeatMasker .align (9d0105e5).
Each alignment record header line is followed, a few lines later, by 'Kimura (with divCpGMod) = x'.
The family is the repeat name in the header (the field 'name#class/superfamily'); output: family, class, n, median.
Usage: young_kimura.py genome.fa.align OUT.tsv"""
import sys, re, statistics, collections
kim = collections.defaultdict(list); cls = {}
fam = None
head = re.compile(r"^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+\S+\s+\d+\s+\d+\s+\(\d+\)\s+(C\s+)?(\S+)#(\S+)")
with open(sys.argv[1]) as f:
    for line in f:
        m = head.match(line)
        if m:
            fam = m.group(2); cls[fam] = m.group(3); continue
        if line.startswith("Kimura (with divCpGMod)") and fam:
            kim[fam].append(float(line.split("=")[1])); fam = None
with open(sys.argv[2], "w") as o:
    o.write("family\tclass\tn_records\tmedian_kimura_cpg\n")
    for k in sorted(kim):
        o.write(f"{k}\t{cls[k]}\t{len(kim[k])}\t{statistics.median(kim[k]):.3f}\n")
print("families", len(kim), "records", sum(map(len, kim.values())))
