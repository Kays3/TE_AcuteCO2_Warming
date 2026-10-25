#!/usr/bin/env python3
"""Stage 3: detection, self-driven / passenger classes, family outcomes, shuffle p, and the family sums for Stage 4.

Rules: PREREGISTRATION sections 4, 5, 7, 8 with amendments A2 (points 7-12) and A3. Reads the featureCounts tables of
Stage 2. Unique fragments drive every locus call; 1/NH counts enter family sums only.
"""
import argparse, glob, gzip, os, sys
import numpy as np
import pandas as pd

MIN_FRAG, MIN_LIBS = 5, 3          # detection: >= 5 unique fragments in >= 3 libraries of a tissue
SELF, PASS = 3.0, 2.0              # section 5
GAP_FRAC, GAP_MIN_BP = 0.25, 200   # A2.8
TISSUES = ("brain", "digestivetract", "liver", "muscle")  # warming E1 (A1.5): the only change from the mobilisation script


def read_fc(path):
    """featureCounts table -> (ids, lengths, counts). Rows come in SAF order; duplicates impossible (unique ids)."""
    df = pd.read_csv(path, sep="\t", comment="#", usecols=[0, 5, 6], dtype={0: str})
    df.columns = ["id", "len", "n"]
    return df


def fpkb(n, length):
    return n / (length / 1000.0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--stage2", required=True)
    p.add_argument("--gene-counts", required=True, help="folder of OA *.cntTable")
    p.add_argument("--out", required=True)
    p.add_argument("--samples", help="test only: a subset list; the real run uses stage2/samples.txt")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)

    samples = [l.strip() for l in open(a.samples or f"{a.stage2}/samples.txt") if l.strip()]
    tissue = np.array([s.split("_")[0] for s in samples])
    assert a.samples or set(tissue) == set(TISSUES), set(tissue)
    ctx = pd.read_csv(f"{a.stage2}/locus_context.tsv.gz", sep="\t", dtype={"rm_id": str})
    ctx["lid"] = "L:" + ctx.rm_id
    nlib = len(samples)

    # ---- unique pass: loci, flanks, gaps, intron bins, shuffles 0-98 ----
    first = read_fc(f"{a.stage2}/fc/{samples[0]}.uniq.txt.gz")
    ids, lens = first.id.values, first.len.values.astype(float)
    kind = pd.Series(ids).str.split(":").str[0].values
    is_s = pd.Series(kind).str.match(r"^S\d+$").values
    keep = ~is_s
    U = np.zeros((keep.sum(), nlib), dtype=np.int32)          # non-shuffle regions x libraries
    shuf_ids = ids[is_s]
    S_ge = {t: np.zeros(is_s.sum(), dtype=np.int16) for t in TISSUES}  # libraries with >= MIN_FRAG, per tissue
    for j, s in enumerate(samples):
        df = first if j == 0 else read_fc(f"{a.stage2}/fc/{s}.uniq.txt.gz")
        assert (df.id.values == ids).all(), s
        n = df.n.values
        U[:, j] = n[keep]
        S_ge[tissue[j]] += (n[is_s] >= MIN_FRAG)
    rid = ids[keep]; rlen = lens[keep]; rkind = kind[keep]
    row = {r: i for i, r in enumerate(rid)}

    # ---- 1/NH pass (loci only) ----
    M = np.zeros((len(ctx), nlib), dtype=np.float64)
    for j, s in enumerate(samples):
        df = read_fc(f"{a.stage2}/fc/{s}.multi.txt.gz").set_index("id")
        M[:, j] = df.loc[ctx.lid, "n"].values

    # ---- A3 extra shuffles (unique pass) ----
    extra = sorted(glob.glob(f"{a.stage2}/fc_extra/*.extra.txt.gz"))
    have_extra = len(extra) == nlib
    if have_extra:
        e0 = read_fc(f"{a.stage2}/fc_extra/{samples[0]}.extra.txt.gz")
        ex_ids = e0.id.values
        E_ge = {t: np.zeros(len(ex_ids), dtype=np.int16) for t in TISSUES}
        for j, s in enumerate(samples):
            df = e0 if j == 0 else read_fc(f"{a.stage2}/fc_extra/{s}.extra.txt.gz")
            assert (df.id.values == ex_ids).all(), s
            E_ge[tissue[j]] += (df.n.values >= MIN_FRAG)

    # ---- gene expression (A2.8a) ----
    G = {}
    for j, s in enumerate(samples):
        g = pd.read_csv(f"{a.gene_counts}/{s}.cntTable", sep="\t", header=0, names=["gene", "n"])
        g = g[~g.gene.str.contains(":")]
        G[s] = g.set_index(g.gene.str.strip('"'))["n"]
    G = pd.DataFrame(G).fillna(0)
    gene_expr = {t: ((G.loc[:, tissue == t] >= MIN_FRAG).sum(axis=1) >= MIN_LIBS) for t in TISSUES}

    # ---- per-locus calls ----
    Li = np.array([row[l] for l in ctx.lid])
    Lc = U[Li]                                              # unique fragments, loci x libraries
    llen = rlen[Li]
    def side(prefix):
        idx = np.array([row.get(f"{prefix}:{r}", -1) for r in ctx.rm_id])
        return idx
    fl, fr, gp = side("FL"), side("FR"), side("GAP")
    # intron bins grouped by locus
    ib_rows = np.where(rkind == "IB")[0]
    ib_locus = pd.Series(rid[ib_rows]).str.split(":").str[1].values
    ib_groups = pd.Series(ib_rows).groupby(ib_locus).apply(lambda x: x.values).to_dict()

    out = ctx[["rm_id", "chrom", "start", "end", "tier", "family", "context", "host_gene", "n_intron_bins", "gap_len"]].copy()
    out["unique_total"] = Lc.sum(axis=1)
    out["multi_total"] = M.sum(axis=1)
    out["not_resolvable"] = out.unique_total == 0
    for t in TISSUES:
        cols = np.where(tissue == t)[0]
        lib_ok = Lc[:, cols] >= MIN_FRAG
        nl = lib_ok.sum(axis=1)
        det = nl >= MIN_LIBS
        out[f"{t}_n_libs"] = nl
        out[f"{t}_detected"] = det
        # pooled over the libraries in which the locus reaches MIN_FRAG (section 5)
        pool = lambda idx: np.where(idx >= 0, (U[np.maximum(idx, 0)][:, cols] * lib_ok).sum(axis=1), np.nan)
        e_n = (Lc[:, cols] * lib_ok).sum(axis=1)
        e = fpkb(e_n, llen)
        f_l = fpkb(pool(fl), np.where(fl >= 0, rlen[np.maximum(fl, 0)], np.nan))
        f_r = fpkb(pool(fr), np.where(fr >= 0, rlen[np.maximum(fr, 0)], np.nan))
        g_n = pool(gp)
        g_c = fpkb(g_n, np.where(gp >= 0, rlen[np.maximum(gp, 0)], np.nan))
        g_libs = np.where(gp >= 0, ((U[np.maximum(gp, 0)][:, cols] >= 1) & lib_ok).sum(axis=1), 0)
        host = np.full(len(ctx), np.nan)
        for k in np.where(det & (ctx.context.values == "intronic"))[0]:
            rows_k = ib_groups.get(ctx.rm_id.values[k])
            if rows_k is not None and len(rows_k):
                b = fpkb((U[rows_k][:, cols] * lib_ok[k]).sum(axis=1), rlen[rows_k])
                host[k] = np.median(b)
        r_l, r_r = (e + 1) / (f_l + 1), (e + 1) / (f_r + 1)
        r_big = (e + 1) / (np.fmax(f_l, f_r) + 1)
        r_host = (e + 1) / (host + 1)
        expr = ctx.host_gene.map(lambda g: bool(gene_expr[t].get(g, False)) if isinstance(g, str) else False).values
        down = ctx.context.values == "downstream5kb"
        gap_len = ctx.gap_len.fillna(0).values
        continuous = np.where(gap_len < GAP_MIN_BP, True, (g_c >= GAP_FRAC * e) & (g_libs >= MIN_LIBS))
        passenger = ((ctx.context.values == "exonic")
                     | (r_big < PASS)
                     | ((ctx.context.values == "intronic") & (r_host < PASS))
                     | (down & expr & continuous))
        for thr, tag in ((SELF, ""), (2.0, "_2x"), (5.0, "_5x")):
            selfd = (np.isin(ctx.context.values, ["intergenic", "intronic"])
                     & ~np.isnan(f_l) & ~np.isnan(f_r) & (r_l >= thr) & (r_r >= thr)
                     & ((ctx.context.values != "intronic") | (~np.isnan(host) & (r_host >= thr)))
                     & ~passenger)
            cls = np.where(selfd, "self_driven", np.where(passenger, "passenger", "undetermined"))
            out[f"{t}_class{tag}"] = np.where(det, cls, "")
        out[f"{t}_elem_fpkb"] = np.where(det, e, np.nan)
        out[f"{t}_flank_l_fpkb"] = np.where(det, f_l, np.nan)
        out[f"{t}_flank_r_fpkb"] = np.where(det, f_r, np.nan)
        out[f"{t}_host_intron_fpkb"] = np.where(det, host, np.nan)
        out[f"{t}_gap_fpkb"] = np.where(det, g_c, np.nan)
        out[f"{t}_host_gene_expressed"] = expr
    out.to_csv(f"{a.out}/locus_calls.tsv.gz", sep="\t", index=False)

    # ---- family outcomes (section 8; registered on T1; T1u and SINE-T reported the same way, descriptively) ----
    fam_rows = []
    nt = out[out.tier != "C1"]
    for intact, label in ((["T1"], "T1"), (["T1u"], "T1u"), (["SINE-T"], "SINE-T")):
        for fam, d in nt.groupby("family"):
            it = d[d.tier.isin(intact)]
            t2 = d[d.tier == "T2"]
            if not len(it) and label != "T1":
                continue
            rec = dict(family=fam, intact_tier=label, n_intact=len(it), n_T2=len(t2))
            n_self = 0
            for t in TISSUES:
                s_t = (it[f"{t}_class"] == "self_driven").any()
                n_self += s_t
                rec[f"{t}_self"] = s_t
            det_any = any(pd.concat([it, t2])[f"{t}_detected"].any() for t in TISSUES)
            if not len(it):
                o = "NO_INTACT_COPIES"
            elif n_self:
                o = "SELF_TRANSCRIBED_INTACT"
            elif det_any:
                o = "TRANSCRIBED_PASSENGER_ONLY"
            elif pd.concat([it, t2]).unique_total.sum() == 0:
                o = "NOT_RESOLVABLE"
            else:
                o = "INTACT_NOT_DETECTED"
            rec.update(outcome=o, n_tissues_self=n_self)
            fam_rows.append(rec)
    pd.DataFrame(fam_rows).to_csv(f"{a.out}/family_outcomes.tsv", sep="\t", index=False)

    # ---- detection rates by tier and context (control C1) ----
    rates = []
    for t in TISSUES:
        for (tier, c), d in out.groupby(["tier", "context"]):
            rates.append(dict(tissue=t, tier=tier, context=c, n=len(d), detected=int(d[f"{t}_detected"].sum()),
                              not_resolvable=int(d.not_resolvable.sum())))
    pd.DataFrame(rates).to_csv(f"{a.out}/detection_rates.tsv", sep="\t", index=False)

    # ---- shuffle p (section 7.3, A2.11, A3) ----
    sh = pd.DataFrame({"sid": shuf_ids})
    sh["k"] = sh.sid.str.split(":").str[0].str[1:].astype(int)
    sh["rm_id"] = sh.sid.str.split(":").str[1]
    for t in TISSUES:
        sh[t] = S_ge[t] >= MIN_LIBS
    if have_extra:
        ex = pd.DataFrame({"sid": ex_ids})
        ex["k"] = ex.sid.str.split(":").str[0].str[1:].astype(int)
        ex["rm_id"] = ex.sid.str.split(":").str[1]
        for t in TISSUES:
            ex[t] = E_ge[t] >= MIN_LIBS
        sh = pd.concat([sh, ex], ignore_index=True)
    tier_of = out.set_index("rm_id").tier
    sh["tier"] = sh.rm_id.map(tier_of)
    real = out[out.rm_id.isin(set(sh.rm_id))]
    prow = []
    for t in TISSUES:
        for tier in ("T1", "SINE-T", "T1u", "pooled"):
            rl = real if tier == "pooled" else real[real.tier == tier]
            if not len(rl):
                continue
            rate = rl[f"{t}_detected"].mean()
            ss = sh if tier == "pooled" else sh[sh.tier == tier]
            per_k = ss.groupby("k")[t].mean()
            for K, lab in ((99, "registered_99"), (999, "a3_999")):
                pk = per_k[per_k.index < K]
                if len(pk) < K:
                    continue
                prow.append(dict(tissue=t, tier=tier, n_loci=len(rl), real_rate=rate, shuffles=K,
                                 shuffle_mean=pk.mean(), shuffle_max=pk.max(),
                                 n_at_or_above=int((pk >= rate).sum()), p=(1 + (pk >= rate).sum()) / (K + 1)))
    pr = pd.DataFrame(prow)
    pr["bonferroni_6"] = np.where(pr.tier.isin(["T1", "SINE-T"]), np.minimum(pr.p * 6, 1.0), np.nan)
    pr.to_csv(f"{a.out}/shuffle_p.tsv", sep="\t", index=False)

    # ---- family sums for Stage 4 (T1, T1u, SINE-T, T2 loci; unique and 1/NH) ----
    sel = (ctx.tier != "C1").values
    fam = ctx.family.values[sel]
    pd.DataFrame(Lc[sel], columns=samples).groupby(fam).sum().to_csv(f"{a.out}/family_sums_unique.tsv.gz", sep="\t")
    pd.DataFrame(M[sel], columns=samples).groupby(fam).sum().to_csv(f"{a.out}/family_sums_multi.tsv.gz", sep="\t")
    selc = (ctx.tier == "C1").values
    pd.DataFrame(M[selc], columns=samples).groupby(ctx.family.values[selc]).sum().to_csv(
        f"{a.out}/family_sums_multi_C1.tsv.gz", sep="\t")
    pd.DataFrame(Lc[sel], columns=samples, index=ctx.rm_id.values[sel]).to_csv(
        f"{a.out}/locus_unique_noC1.tsv.gz", sep="\t")
    print("extra shuffles present:", have_extra)
    print(pd.DataFrame(fam_rows).groupby(["intact_tier", "outcome"]).size().to_string())
    print(pr.to_string())


if __name__ == "__main__":
    sys.exit(main())
