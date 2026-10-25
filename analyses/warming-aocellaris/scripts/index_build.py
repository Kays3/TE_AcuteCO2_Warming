#!/usr/bin/env python3
"""5.2.2 global TE index for PRJNA970523 (warming). Copied from external-co2-stress-validation-20260907/
scripts/build_index_and_matrix.py; only the input glob and output paths changed. is_strict_te needs no change:
the A. ocellaris TEcount ids carry the class as their last field, as in OA.
Index definition copied verbatim in logic from paper/results/te-uniqueread-sensitivity-20260903/
scripts/e4_index.py: sum(strict-TE counts)/sum(all counts)*100 per library."""
import os, glob
CNT = os.path.expanduser("~/jim_te_warm/out")
STRICT = {"DNA", "LTR", "LINE", "SINE", "RC", "PLE", "Retroposon"}

def is_strict_te(fid):
    parts = fid.split(":")
    return len(parts) >= 3 and parts[-1] in STRICT

files = sorted(glob.glob(f"{CNT}/*/*.multi.cntTable"))
samples, idx, mat = [], [], {}
for fp in files:
    sid = os.path.basename(fp).replace(".multi.cntTable", "")
    samples.append(sid)
    st = 0.0; tot = 0.0
    with open(fp) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2: continue
            try: c = float(p[1])
            except ValueError: continue
            tot += c
            if is_strict_te(p[0]):
                st += c
                mat.setdefault(p[0], {})[sid] = int(c)
    idx.append((sid, st, tot, 100.0 * st / tot if tot else 0.0))

with open(os.path.expanduser("~/jim_te_warm/analysis/te_index_per_sample.tsv"), "w") as o:
    o.write("sample_id\tstrict_te_sum\ttotal_counts\tte_index_pct\n")
    for sid, s, t, fr in idx:
        o.write(f"{sid}\t{s:.1f}\t{t:.1f}\t{fr:.6f}\n")

with open(os.path.expanduser("~/jim_te_warm/analysis/te_feature_counts.tsv"), "w") as o:
    o.write("feature\t" + "\t".join(samples) + "\n")
    for feat in sorted(mat):
        row = mat[feat]
        o.write(feat + "\t" + "\t".join(str(row.get(s, 0)) for s in samples) + "\n")
print(f"samples={len(samples)} te_features={len(mat)}")
