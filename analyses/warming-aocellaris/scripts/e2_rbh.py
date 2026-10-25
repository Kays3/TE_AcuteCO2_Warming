"""E2 (A1.6): one-to-one reciprocal best hits passing 80-80-80 on the best HSP.
Best hit per query = highest bitscore. A pair passes if, on that HSP, identity >= 80 %, alignment length >= 80 bp,
and the HSP covers >= 80 % of both the query and the subject length.
Usage: e2_rbh.py AO_VS_AP.tsv AP_VS_AO.tsv OUT.tsv"""
import sys
import pandas as pd
cols = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()
def best(f):
    d = pd.read_csv(f, sep="\t", names=cols)
    d = d.sort_values(["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid")
    d["qcov"] = (d.qend - d.qstart + 1) / d.qlen
    d["scov"] = ((d.send - d.sstart).abs() + 1) / d.slen
    return d
a, b = best(sys.argv[1]), best(sys.argv[2])
m = a.merge(b[["qseqid", "sseqid"]].rename(columns={"qseqid": "sseqid", "sseqid": "back"}), on="sseqid")
m = m[m.back == m.qseqid]
ok = m[(m.pident >= 80) & (m.length >= 80) & (m.qcov >= 0.8) & (m.scov >= 0.8)]
ok[["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcov", "scov", "bitscore"]].rename(
    columns={"qseqid": "aocel_family", "sseqid": "apoly_family"}).to_csv(sys.argv[3], sep="\t", index=False)
print(f"queries with hits {len(a)}; reciprocal best {len(m)}; passing 80-80-80 {len(ok)}")
