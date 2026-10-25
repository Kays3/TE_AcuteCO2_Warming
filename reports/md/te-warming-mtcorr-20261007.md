# Warming and TE expression in *A. ocellaris* after removing mitochondrial rRNA: post-hoc results

**Author:** Kaisar Dauyey, Laboratory of Mathematical Biology, Hokkaido University, Japan  
**Request:** the human, 6 October 2026, about 22:50 JST (13:50Z), passed on by the orchestrator (message e69a82): "for warming check muscle again for positive control and relax the criteria, reevaluate for all tissues based on correction of mt 16S RNA of gene counts".  
**Status:** complete. The work is post-hoc and follows Amendment 5 of the warming pre-registration (`PREREGISTRATION.md` sha256 `ecab474dc713`, committed at 5b09082, passed by Andy and Phyllis before any corrected number existed). Code and tables are on TE_OA branch `analysis/te-warming-mtcorr-20261006`, head `cf83ef6`. For review by Phyllis (science) and Andy (numbers).  
**Date:** 7 October 2026.

## Summary of findings

We asked whether the warming result for transposable elements (TEs) changes when the mitochondrial ribosomal RNA (mt-rRNA) genes are removed from the gene counts, and whether muscle passes a relaxed gene positive control. Under M1, removing the rRNA genes does not change the TE result in any tissue. Muscle passes the relaxed control by a small margin, and its null result can now be read, though as weak evidence.

- **The TE-family calls are the same after correction.** Under the post-hoc primary correction (M1, the 95 rRNA genes dropped), brain has 8 classified TE families at false discovery rate (FDR) 0.05 (permutation p = 0.0022, the smallest value the 924-labeling null can give, so brain's p is at the floor and is not ordered against the other tissues'). Digestive tract has 22 (p = 0.0065), liver 24 (p = 0.0043) and muscle 0 (p = 1.0). These are the same families as the registered fits, with log2 fold changes that rank identically (Spearman ρ = 1.00 in each tissue). Every permutation p is from 924 labelings of two tanks, so each detected shift is a difference between two tanks. (Figure 4; Table 1)
- **Readings.** Brain, digestive tract and liver read `POSTHOC_TANK_TEMPERATURE_SHIFT_DETECTED`. This is exploratory, and the temperature effect cannot be separated from the tank. Muscle reads `POSTHOC_NO_DETECTABLE_SHIFT`. (Figure 2)
- **The muscle control passes by 7 genes under M1, and by 9 on the registered fits.** Muscle has 101 genes at FDR 0.05. The rank-832 value (90th percentile) of its M1 null is 94, a margin of 7, and 90 of 924 labelings reach 101 or more. On the registered fits the observed count is also 101 and the rank-832 value is 92, a margin of 9 (descriptive only); Andy had disclosed that 92 before the run. At the registered 95th percentile (402), muscle still fails. Under the relaxed rule, a tissue's control passes by chance about 10 % of the time instead of 5 %, so the muscle null is weaker evidence than a null read through the registered control. The reading refers to family-level differential expression only. The descriptive global TE index is lower at 31 °C in muscle by −0.879 percentage points [−1.386, −0.372], the largest drop of any tissue, and like every comparison here it is a difference between two tanks. A shift spread thinly across many families can fail the family-by-family FDR test with six libraries per group while the aggregate index picks it up. The muscle reading therefore does not mean that muscle TE expression is unchanged. (Figures 5 and 6)
- **mt-rRNA is a small share of gene counts in these libraries.** The 12S and 16S genes together take 0.03–0.30 % of gene counts per library. Brain is highest (0.23–0.30 %) and the other tissues are at 0.11 % or less. That is too little for removing them to move the size factors noticeably. (Figure 3a)
- **The covariate sensitivity (M2) is confounded with tank in the digestive tract.** In the digestive tract the mt-rRNA covariate separates the two tanks (point-biserial r = 0.86). With it in the model, the TE count falls to 0 and the gene count to 25, below its own rank-832 value of 74. M2 gives no reading, by design. In the digestive tract it cannot separate composition from tank. (Figure 3b; Table 1)
- **The global TE index barely moves.** With rRNA removed from the denominator, the 31 minus 28 °C differences are within 0.005 percentage points of the registered ones. Liver moves from −0.144 to −0.147 points [−0.187, −0.108] and muscle from −0.875 to −0.879 [−1.386, −0.372], n = 12 per tissue. (Figure 6)
- **What the data cannot show.** No analysis here can tell a temperature effect from a tank effect. These results do not raise the registered readings to confirmatory status, or change them.

## What was asked

The human asked for two things. The first was to correct the gene counts for mt 16S rRNA and re-evaluate every tissue. The second was to check muscle against the gene positive control again, with a relaxed criterion. In the registered analysis, muscle had no TE family at FDR 0.05. Its gene control also failed (101 genes, against a registered 95th percentile of 402), so that null could not be read.

Amendment 5 sets out how this was done, before any corrected number was computed:

- **A5.1.** The rows to remove are the 95 genes with `gene_biotype "rRNA"` in the RefSeq annotation: the mitochondrial 16S (`KEG47_r01`) and 12S (`KEG47_r02`) genes, and 93 nuclear rRNA genes.
- **A5.3 (M1).** M1 drops these rows and reruns the registered procedure unchanged. It is the post-hoc primary.
- **A5.4 (M2).** M2 adds the per-library mt-rRNA log ratio as a covariate. It is a sensitivity that yields no reading.
- **A5.5 (relaxed control).** The control is applied per tissue, at the 90th percentile of the tissue's own 924-labeling gene null, where the registered threshold was the 95th.
- **A5.6 (readings).** The reading rules are the registered ones with a `POSTHOC_` prefix.
- **A5.7 (other tests).** The superfamily test (E3) and the global TE index are redone; E1, E2, 5.4 and the 5.2.4 tissue model are not.
- **Ruled out.** The request ruled out new alignment. Everything runs from the existing TEcount tables, under the study's 250 CPU-hour (CPU-h) ceiling.

## Methods

**Data.** The data are 48 poly(A) RNA-seq libraries from PRJNA970523. The design is four tissues (brain, digestive tract, liver, muscle) × two temperatures (28 °C and 31 °C) × 6 fish. Each temperature group is one tank. Counts are the existing TEcount 2.2.4 multi and unique tables from the registered run. That run used STAR 2.7.11b, reverse-stranded, against GCF_022539595.1 with the RefSeq annotation and the RepeatMasker TE annotation (`RUN_RECORD.md`).

```{.dot caption="Figure 1. Processing for Amendment 5. Only the steps in the dashed box are new; the counts are the registered run's. Tool versions and counts from TE_OA:paper/results/te-warming-aocellaris-20261003/RUN_RECORD.md at cf83ef6." width="95%"}
digraph G {
  rankdir=LR; nodesep=0.25; ranksep=0.35;
  node [shape=box, style=rounded, fontsize=10, color="#222222"]; edge [color="#222222", fontsize=9];
  R [label="48 TEcount tables\n(registered run; STAR 2.7.11b,\nTEcount 2.2.4)"];
  G [label="RefSeq GTF\ngene_biotype rRNA\n95 genes (2 mt, 93 nuclear)"];
  subgraph cluster_new { label="Amendment 5 (new, 15.98 CPU-h)"; fontsize=9; color="#8C8C8C"; style=dashed;
    S [label="A5.2 shares\nper library"];
    M1 [label="M1: drop 95 rows\nfilter ≥10 in ≥6, poscounts\nDESeq2 1.46.0 ~ temperature"];
    M2 [label="M2: M1 + mt-rRNA covariate\n~ mtrRNA + temperature"];
    E3 [label="E3 under M1\nsuperfamily level"];
    N [label="924 labelings each\nBH over classified families\ngene count at FDR 0.05"];
    I [label="global TE index\nrRNA out of denominator"]; }
  O [label="readings (A5.6)\nrelaxed control (A5.5)"];
  R -> S; G -> S; R -> M1; G -> M1; M1 -> M2; M1 -> E3; M1 -> N; M2 -> N; E3 -> N; N -> O; R -> I; G -> I;
}
```

**Statistics.** The model is DESeq2 1.46.0 on R 4.4.3, design `~ temperature` (M1) or `~ mtrRNA + temperature` (M2), with `poscounts` size factors and independent filtering off. Benjamini–Hochberg (BH) FDR is computed over the classified TE families (DNA, LINE, SINE, LTR, RC, PLE, Retroposon).

The null is exact: all 924 ways of splitting each tissue's 12 libraries into two groups of six, each fitted with the same function as the observed labeling. The permutation p is the share of the 924 labelings whose family count at FDR 0.05 is at least the observed count, and a shift reads as detected at p < 0.0125.

The relaxed gene control passes when the observed gene count at FDR 0.05 is strictly greater than the value at rank ⌈0.9 × 924⌉ = 832 of that tissue's ordered null counts. The M2 covariate is log((KEG47_r01 + KEG47_r02 + 0.5) / (non-rRNA gene counts + 0.5)), centred within tissue. The global TE index is strict-TE counts as a percentage of all counts. It is compared between temperatures with a Welch interval, Hedges' g, and Welch and Mann–Whitney tests.

```{.dot caption="Figure 2. Reading rules of Amendment 5 (A5.6) with the path each tissue took in bold vermillion. Brain, digestive tract and liver leave at the first test; muscle passes through the relaxed control. Rules from TE_OA:paper/results/te-warming-aocellaris-20261003/PREREGISTRATION.md at 5b09082; values on the edges from mtcorr/tables/a5_results.tsv at ea719be." width="92%"}
digraph G {
  rankdir=LR; nodesep=0.3; ranksep=0.5;
  node [shape=box, style=rounded, fontsize=10, color="#222222"]; edge [color="#8C8C8C", fontsize=9];
  A [label="M1 fit, per tissue\n924-labeling null"];
  Q1 [label="≥ 1 classified family at FDR 0.05\nand p < 0.0125?", shape=diamond, fontsize=9];
  Q2 [label="genes at FDR 0.05 >\nrank-832 value of own null?", shape=diamond, fontsize=9];
  D [label="POSTHOC_TANK_TEMPERATURE_\nSHIFT_DETECTED\n(exploratory; tank clause)"];
  N [label="POSTHOC_NO_DETECTABLE_SHIFT\n(weaker evidence: 10 % chance pass)"];
  U [label="POSTHOC_UNINFORMATIVE"];
  A -> Q1 [color="#D55E00", penwidth=2.2];
  Q1 -> D [label="brain 8, p .0022\nDT 22, p .0065\nliver 24, p .0043", color="#D55E00", penwidth=2.2, fontcolor="#D55E00"];
  Q1 -> Q2 [label="muscle 0, p 1.0", color="#D55E00", penwidth=2.2, fontcolor="#D55E00"];
  Q2 -> N [label="muscle 101 > 94", color="#D55E00", penwidth=2.2, fontcolor="#D55E00"];
  Q2 -> U [label="no"];
}
```

**Quality control.** No library was excluded. The registered exclusion check (all 48 above a 50 % uniquely mapped rate) still applies, because the alignments are unchanged. The completeness of the new nulls was checked: every tissue-by-fit set has all labelings 0–924, and none of the 148 job logs shows an error or a guard stop. The index reimplementation reproduces the registered numerator and denominator exactly (maximum difference 0) before the rRNA rows are removed.

## Results

### How much of the gene counts is mt-rRNA

![Figure 3. (a) Mitochondrial 12S + 16S rRNA counts as a share of gene counts, per library; (b) the M2 covariate (centred log ratio of mt-rRNA to non-rRNA gene counts), with the point-biserial correlation r between covariate and tank above each tissue. Blue 28 °C tank, vermillion 31 °C tank; bars are group means. n = 48 libraries (6 per group). Units: % of gene counts (a); natural-log ratio, centred within tissue (b). Descriptive. Source: `TE_OA:paper/results/te-warming-aocellaris-20261003/mtcorr/tables/a52_shares.tsv` at `ea719be`; drawn by `hive/reports/py/te-warming-mtcorr-20261007.py`.](../png/te-warming-mtcorr-20261007/fig1_shares.png){width=78%}

The mt-rRNA genes are a small part of the gene counts in this poly(A) data (Figure 3a). Brain is highest, with group means of 0.256 % (28 °C) and 0.263 % (31 °C). The other tissues sit between 0.039 % and 0.081 %. The 16S gene alone is 0.036–0.193 % by group, and all 95 rRNA genes together are 0.130–0.337 %. With shares this small, and `poscounts` size factors that resist a few dominant rows, we expected M1 to change little (A5.4 says so in advance).

The covariate behaves differently across tissues (Figure 3b). In the digestive tract it separates the tanks almost completely, with centred means −0.164 at 28 °C and +0.164 at 31 °C, r = 0.86. Liver has a weaker, opposite association (r = −0.48), and brain and muscle show little (r = 0.15 and 0.18). Where the covariate tracks the tank, M2 removes tank variation along with any compositional effect. These shares cannot tell us which of the two it removed.

### TE families under the correction

![Figure 4. Permutation nulls for the M1 fits: each grey point is the number of classified TE families at FDR 0.05 for one of the 924 labelings; vermillion line is the observed labeling, with the count of labelings at or above it. n = 924 labelings per tissue. Units: families, symmetric-log axis. Test: exact permutation over 924 labelings, threshold p < 0.0125. Source: `TE_OA:paper/results/te-warming-aocellaris-20261003/mtcorr/tables/null_mtcorr.tsv.gz` at `ea719be`; drawn by `hive/reports/py/te-warming-mtcorr-20261007.py`.](../png/te-warming-mtcorr-20261007/fig2_te_null_m1.png){width=88%}

Under M1 the observed family counts sit in the upper tail of the null in brain (2 of 924 labelings at or above 8), digestive tract (6 of 924 at or above 22) and liver (4 of 924 at or above 24) (Figure 4). In muscle, all 924 labelings reach the observed 0 (p = 1.0; Figure 4d). In this 6-vs-6 exact null each labeling and its complement give the same family count (checked for all 462 pairs in every tissue), so the smallest attainable p is 2/924 = 0.0022, and brain sits exactly at it.

These are the registered counts and permutation p values, and the same families. No family was gained or lost in any tissue, and the log2 fold changes of the tested families rank identically (Spearman ρ = 1.00 to two decimals). Most detected families are lower at 31 °C: 7 of 8 in brain, 18 of 22 in digestive tract and 20 of 24 in liver (`mtcorr/tables/m1_families_<tissue>.tsv`).

The unique-mode sensitivity under M1 gives 7, 24, 23 and 3 families, as in the registered unique fits. The superfamily test (E3) under M1 also matches its registered result:

- digestive tract: 5 superfamilies, p = 0.0043;
- liver: 2, p = 0.039;
- brain and muscle: 0.

| Fit | Brain | Digestive tract | Liver | Muscle |
|---|---|---|---|---|
| Registered (families; p) | 8; .0022 | 22; .0065 | 24; .0043 | 0; 1.0 |
| M1, post-hoc primary | 8; .0022 | 22; .0065 | 24; .0043 | 0; 1.0 |
| M2, covariate (no reading) | 5; .0108 | 0; 1.0 | 17; .0130 | 0; 1.0 |
| M1 genes at FDR 0.05 / rank-832 value | 918 / 57 | 1,013 / 43 | 1,435 / 73 | 101 / 94 |
| M2 genes at FDR 0.05 / rank-832 value | 1,035 / 59 | 25 / 74 | 760 / 80 | 84 / 105 |

Table 1. Classified TE families at FDR 0.05 with permutation p (924 labelings), and gene counts against the rank-832 value of each fit's own null. The exact null is one-sided over 462 distinct partitions, each labeling and its complement giving the same count, so .0022 is the floor: a p at the floor is not ordered against the others, and only the family counts separate those tissues. Source: `TE_OA:paper/results/te-warming-aocellaris-20261003/mtcorr/tables/a5_results.tsv` at `ea719be`.

M2 moves the counts most in the digestive tract, where the covariate tracks the tank (Table 1). The TE count there falls from 22 to 0. The gene count falls from 1,013 to 25, and fails the relaxed control at that tissue's rank-832 value of 74. A model that loses almost all of a tissue's gene signal when one tank-correlated covariate enters is not evidence that the TE shift was compositional; the same loss would follow if the covariate absorbed most of the tank difference. In liver, M2 keeps 17 families at p = 0.0130, just above the 0.0125 threshold. In brain it keeps 5 at p = 0.0108. M2 yields no reading under Amendment 5, and these numbers do not change the M1 readings.

### Gene positive control, relaxed

![Figure 5. Gene positive control per tissue: each grey point is the number of genes at FDR 0.05 for one of 924 labelings; vermillion line is the observed labeling; dashed line is the rank-832 value (90th percentile), the A5.5 threshold. Top row (a–d): registered fits and their stored nulls, descriptive only. Bottom row (e–h): M1 fits, which carry the A5.6 reading. Each panel prints the observed count, the number of labelings at or above it, and the rank-832 value. n = 924 labelings per panel. Units: genes, symmetric-log axis. Test: observed count strictly greater than the rank-832 value. Source: `TE_OA:paper/results/te-warming-aocellaris-20261003/tables/null_distributions.tsv.gz` at `266db41` (top) and `mtcorr/tables/null_mtcorr.tsv.gz` at `ea719be` (bottom); drawn by `hive/reports/py/te-warming-mtcorr-20261007.py`.](../png/te-warming-mtcorr-20261007/fig3_gene_control.png){width=100%}

Brain, digestive tract and liver pass the relaxed control by wide margins: under M1, 918 genes against 57, 1,013 against 43 and 1,435 against 73 (Figure 5). They also passed the registered 95th-percentile control.

Muscle is the case the request was about. Under M1, the observed muscle gene count of 101 is above the rank-832 value of 94, a margin of 7 genes. Ninety of the 924 labelings give 101 or more (45 distinct partitions, since the gene counts are also complement-exact), so the observed labeling sits just inside the top 10 % of its null. On the registered fits the observed count is also 101 and the rank-832 value is 92, a margin of 9 (Figure 5d, descriptive only). Andy had told me of the value 92 during the amendment review, before the run; the rule was fixed before that message (PREREGISTRATION, "Seen during review").

Muscle therefore reads `POSTHOC_NO_DETECTABLE_SHIFT`. Three qualifications go with that reading:

- **Weaker evidence.** The control was relaxed from the 95th to the 90th percentile. Under null labels a tissue's control passes by chance about 10 % of the time instead of 5 %, so this null is weaker evidence than one read through the registered control.
- **Small margin.** The margin is small in absolute terms. At the registered 95th percentile (402 genes), muscle still fails.
- **Family level only.** The reading concerns differential expression of single TE families. The descriptive global TE index is lower at 31 °C in muscle by −0.879 percentage points [−1.386, −0.372] (Hedges' g −2.06, Welch p 0.0032; Figure 6), the largest drop of any tissue and about a fifth of a 4.6 % index. Like every comparison here, it is a difference between two tanks. The family test is not designed to detect a shift shared thinly across many families: with six libraries per group, such a shift can miss FDR 0.05 in every family while the aggregate index detects it.

The muscle gene signal is modest in either form. Taken with the index, the muscle reading says that no single TE family shifts detectably. It does not say that muscle TE expression is unchanged.

The heat-shock comparison with the source paper, repeated under M1, is unchanged. The registered heat-shock and DnaJ genes at FDR 0.05 are 8 in brain, 11 in digestive tract, 11 in liver and 4 in muscle, the same counts as before. The largest change in any log2 fold change is 0.00016 (`mtcorr/tables/heatshock_m1_summary.tsv` at `cf83ef6`).

### Global TE index

![Figure 6. Global strict-TE index, 31 °C minus 28 °C, with Welch 95 % confidence intervals: registered (filled) and with the 95 rRNA gene rows removed from the denominator (open). n = 12 libraries per tissue (6 per tank). Units: percentage points of all counts. Test: Welch interval; descriptive. Source: `TE_OA:paper/results/te-warming-aocellaris-20261003/tables/index_5_2_2.tsv` at `93444f3` (registered) and `mtcorr/tables/a57_index.tsv` at `ea719be`; drawn by `hive/reports/py/te-warming-mtcorr-20261007.py`.](../png/te-warming-mtcorr-20261007/fig4_index.png){width=80%}

Removing rRNA from the denominator raises every index slightly, by up to 0.025 percentage points in the muscle group means, and leaves the temperature differences almost unchanged (Figure 6):

- **brain:** −0.039 points [−0.102, 0.024] (registered −0.039);
- **digestive tract:** −0.068 [−0.169, 0.033] (registered −0.069);
- **liver:** −0.147 [−0.187, −0.108] (registered −0.144), Hedges' g −4.41, Welch p 8.7 × 10⁻⁶;
- **muscle:** −0.879 [−1.386, −0.372] (registered −0.875), g −2.06, Welch p 0.0032.

The muscle Mann–Whitney p moves from 0.0043 to 0.0087. One pair of libraries changed order: a 31 °C library at 4.464 % and a 28 °C library at 4.462 % swapped ranks, a difference of 0.002 points. That shift reflects a rank-based test on six-per-group data and says nothing about the size of the effect. The index is descriptive in the registered design and remains so here.

## Limits

- **Tank and temperature cannot be separated.** Each temperature group is one tank of six fish. Every reading, registered or post-hoc, describes a difference between two tanks, and no correction for composition can change that.
- **Post-hoc status.** Amendment 5 was written after the registered results were known, at the human's request. Its readings are exploratory and do not replace the registered ones. The relaxed control in particular was introduced because muscle failed the registered one. The rule was fixed before the muscle null below the 95th percentile was seen (A5.5), but the motive for relaxing it was known.
- **Brain's permutation p is saturated and the three tissues' p values are not ordered.** The exact null is one-sided over 462 distinct partitions, each labeling and its complement giving the same family count, so the smallest attainable p is 2/924 = 0.0022 and brain sits exactly at it (line 96; Figure 4a). Digestive tract at 0.0065 = 6/924 and liver at 0.0043 = 4/924 are two steps and one step above the floor, so where the summary and Table 1 print the three together, brain's value carries no information beyond "no labeling but its own mirror reached the observed count". The family counts, 8 against 22 and 24, are what separate those tissues. The floor sits below the 0.0125 threshold, so this is not a test that could not have reached significance at this n.
- **Power.** With six libraries per group the gene null is wide (the muscle M1 null has a 95th percentile of 402 genes), and a TE shift of the size seen in the other tissues could be missed in muscle. These data do not let us estimate the smallest effect the muscle test could detect.
- **What M1 can and cannot show.** M1 shows that the registered calls do not depend on the rRNA rows. Because these rows are a small share of the counts, it does not test whether library composition more broadly drives the TE shifts. M2 was meant to probe that. In the digestive tract its covariate is confounded with tank, so it cannot answer the question there.
- **Not done.** E1, E2, 5.4 and the 5.2.4 tissue model were not re-evaluated, as Amendment 5 stated in advance. The M2 superfamily fit was not run. No new alignment or counting was done.

## Deviations, compute and records

**Deviations and additions.** Two changes were made relative to Amendment 5:

- **Script fix.** The summary script's glob for the registered nulls first picked up the registered superfamily files as well. It was tightened to the family chunk pattern before any reading was taken (`RUN_RECORD.md`).
- **Where it ran.** The summary ran on a lab workstation (C) instead of the compute workstation (S), because the S Python environment lacks scipy. The inputs were copied unchanged.

Nothing was added after seeing the data beyond what Amendment 5 lists. No registered quantity or reading changes.

![Figure 7. Cumulative compute against the projection and the 250 CPU-h ceiling: the registered study (grey) and the three Amendment 5 jobs, each including its observed fits. n = 148 timed processes for Amendment 5. Units: CPU-h, user + system time from `/usr/bin/time`. Descriptive. Source: `TE_OA:paper/results/te-warming-aocellaris-20261003/RUN_RECORD.md` at `cf83ef6`; drawn by `hive/reports/py/te-warming-mtcorr-20261007.py`.](../png/te-warming-mtcorr-20261007/fig5_compute.png){width=70%}

**Compute.** Amendment 5 used 15.98 CPU-h, below the 17–22 projected, bringing the study to 236.02 of 250 (Figure 7). The jobs were:

- M1 family: 5.911 CPU-h;
- M2 family: 4.534;
- M1 superfamily: 5.526;
- M1 unique observed fits: 0.012.

Every process ran pinned to cores 0–20 with the memory guard on. The M1 observed fits finished at 23:04 JST on 6 October (14:04Z), and the last null at 00:14 JST on 7 October (15:14Z).

**Records.**

- **Repository.** TE_OA, branch `analysis/te-warming-mtcorr-20261006`, head `cf83ef6`. The design is `paper/results/te-warming-aocellaris-20261003/PREREGISTRATION.md` (Amendment 5) at `5b09082`, and the run is in `RUN_RECORD.md` at `cf83ef6`.
- **Scripts** (in `mtcorr/scripts/`): `make_rrna_list.sh`, `deseq_fit_mtcorr.R`, `run_mtcorr.sh`, `summarize_mtcorr.py` and `heatshock_m1.py`.
- **Tables** (in `mtcorr/tables/`): `a5_results.tsv`, `a52_shares.tsv`, `a52_shares_by_group.tsv`, `null_mtcorr.tsv.gz`, `m1_families_<tissue>.tsv`, `a57_index.tsv`, `a57_index_per_library.tsv` and `heatshock_m1_{genes,summary}.tsv`.
- **Registered report.** `hive/reports/md/te-warming-aocellaris-20261003.md`, unchanged.
- **Rebuild.** The figures rebuild with `hive/tools/.venv-report/bin/python hive/reports/py/te-warming-mtcorr-20261007.py`, which copies its tables to `hive/reports/tsv/te-warming-mtcorr-20261007/`. The PDF rebuilds with `hive/tools/md2report.sh hive/reports/md/te-warming-mtcorr-20261007.md`.
