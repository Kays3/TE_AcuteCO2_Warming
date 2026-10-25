# Run record: acute CO2 TE check (PRJNA1102381)

Frozen design: `PREREGISTRATION.md` at 77270d9 (see `REGISTRATION.md`). Machines are named by workstation letter; S is the server that holds the warming study.

## Implementation choices fixed before any download (2026-10-05)

- **Brain gene control null.** The 200 assignments are drawn with `numpy.random.default_rng(20261005)` without replacement from the 599 within-cohort brain assignments that differ from the observed one. The observed assignment is excluded from the null. The threshold is the 190th of the 200 ordered null counts; the control passes if the observed gene count at FDR 0.05 is strictly greater. With 200 draws this threshold is somewhat noisier than one from all 600; that is accepted for a gate.
- **Order if the cap nears** (as the budget table): pipeline, primary TE tests and their permutations, brain gene control, cluster comparison, genic/intergenic split, then the section 5 sensitivities.

## Download (2026-10-05)

26 of 26 runs fetched with `prefetch` 3.4.1 on S, 14:58:48Z to 15:55Z. All are in SRA normalized format with full base qualities. Each passed `vdb-validate`, its spot count equals the NCBI runinfo, and its sha256 and md5 are recorded. The total is 82.3 GB, which is under the 100 GB download cap. Download CPU: 0.31 CPU-h.

## Reference (Tier B stage B2, shared)

The STAR index was built on S at 19:06Z and charged to Tier B. The archive keeper pasted the archive `reference.sha256` and the pipeline command lines, so the archive was not opened. Six files match the archive sha256:

- the packed Ensembl 112 genome and gene GTF, and UCSC danRer11 `rmsk.txt.gz`;
- the unpacked `genome.fa` and `genes.gtf`;
- the TE GTF rebuilt with the archive's `make_te_gtf.py` (3,562,104 loci, 11,255 subfamilies).

## Pipeline deviations and human steps

- **Genome loading.** Each STAR process loads the genome itself, where the archive loads it once into shared memory (`LoadAndKeep`). Loading mode does not change alignments. The orchestrator accepted this at 19:41Z.
- **Temporary files.** Temporary FASTQ files go to three fixed slot folders (`work/slot1` to `work/slot3`), and each new library overwrites its slot. The scripts delete nothing. **Human step:** delete the slot folders and, if wanted, the SRA files after the BAMs are verified.

## Rate check after the first 3 libraries (2026-10-05, 20:23Z)

| Run | CPU-h | Wall | STAR unique | Multi | Assigned (TEcount) |
|---|---|---|---|---|---|
| SRR28763293 | 1.80 | 42:55 | 86.8 % | 7.6 % | 24.9 M |
| SRR28763294 | 1.87 | 43:28 | 90.7 % | 3.9 % | 24.0 M |
| SRR28763295 | 1.79 | 39:28 | 89.3 % | 4.7 % | 24.4 M |

**Steps per library (CPU-h):**
- fasterq-dump 0.04;
- fastp 0.21 to 0.23;
- STAR 1.07 to 1.12;
- TEcount 0.41 to 0.48.

All three libraries pass the QC exclusion rule.

**Rate.** The rate is 0.060 CPU-h per million pairs (2.39 per 40 M), against the design's 1.46.

**Projection.**

| Item | CPU-h |
|---|---|
| Pipeline, 834.9 M pairs | 49.9 |
| Downstream analysis | 2.4 |
| Download | 0.31 |
| **Total** | **about 52.6** |

The stop rule is min(35, 1.5 × 30.5) = 35 CPU-h. The projection passes it, so the runner was paused before the fourth library, and the cap was put to the orchestrator and then to the human. The stated options are a raise to 55 CPU-h at full depth, or a stopped run with no reading.

**Disk.** Projected peak disk is about 252 GB:

| Item | GB |
|---|---|
| SRA files | 82.3 |
| Three slot folders (about 32 GB each: 28 GB uncompressed FASTQ and 4 GB trimmed) | 96.5 |
| 26 BAMs (2.6 to 2.9 GB each) | about 73 |

The preregistration's "about 180 GB" was an underestimate: it assumed 4 × SRA size per working library, and the real figure is about 10 ×. The 100 GB cap applies to downloads and is met. S keeps about 1.45 TB free at the peak, far above the 15 % floor.

## Cap change and resume (appended 2026-10-06)

The human ruled on 06 Oct 2026 (about 11:55 JST, relayed by the orchestrator) to raise the acute cap from 35 to 75 CPU-h and continue at full depth under the frozen design. Nothing in the preregistration changes: same libraries, pipeline, tests and readings. The new stop rule is 75 CPU-h; if the projection passes it, the run stops and the orchestrator is told.

The runner was resumed at about 2026-10-06T03:06:05Z on S. The remaining 23 libraries run three at a time on cores 0 to 20 with the same memory and disk guards. Rates are reported after the next three libraries and at completion. Temporary slot folders are kept, as the human asked.

## Pipeline complete (2026-10-06)

All 26 libraries finished at 08:37:46Z, with no failures and no guard stops. Measured pipeline CPU was 40.47 CPU-h for 834.9 M pairs, or 0.048 CPU-h per million pairs. With the download, the total is 40.8 CPU-h against the 75 cap.

**QC exclusions: none.** Across the 26 libraries:
- STAR uniquely mapped reads were 68.8 to 91.6 %; the exclusion rule is below 30 %.
- TEcount assigned 24.6 to 30.9 M fragments; the rule is below 5 M.
- One gill library, SRR28763310 (EG7, CO2, cohort 1), has 29.0 % multi-mapping reads, against 3.6 to 13.4 % for the other 25. It passes the rule and stays in. The report flags it.

Per-library values are in `tables/qc.tsv`.

## Downstream implementation choices (fixed 2026-10-06, before any downstream fit was run or any count matrix opened)

1. **TE set.** The TEcount TE rows (11,257) carry the rmsk class as the last name field. The tested TE set is every subfamily of class DNA, DNA?, LINE, LTR, SINE or RC, which is 1,011 rows before filtering.
   - Simple_repeat (10,134), Low_complexity, Satellite, the RNA-gene classes (tRNA, snRNA, rRNA, scRNA), ARTEFACT and Unknown are not tested. They also stay out of the size-factor matrix.
   - Unknown is left out as in the damselfish classified set.
   - The budget table's "2,356 TE rows" came from a timing matrix, not from these data.
2. **Rows and size factors.** These are fixed once per tissue and used unchanged in every fit, observed and permuted (one procedure for both).
   - The 08 filter keeps rows with at least 10 counts in at least as many libraries as the smaller group (brain 7, gill 6).
   - Size factors are DESeq2 median-of-ratios on the filtered gene plus TE-set rows. The TE statistic is fitted on the filtered TE-set rows, the brain gene control on the filtered gene rows.
3. **DESeq2.** DESeq2 1.46.0 uses its defaults: parametric dispersion fit, Cook's distance handling and the Wald test. The one exception is `independentFiltering = FALSE`, so that Benjamini–Hochberg runs over every tested row of the set. Effect sizes come from apeglm 1.28.0 on the observed fit.
   - The frozen environment on S has DESeq2 1.46.0 but no apeglm, so a separate environment with the same DESeq2 1.46.0, apeglm 1.28.0 and R 4.4.3 is used for all fits.
4. **Enumeration.** Within a cohort, libraries are sorted by run accession. CO2 positions come from `itertools.combinations` and cohorts are combined with `itertools.product`, cohort 1 varying slowest. Assignment 0 is the observed one.
   - This gives 600 brain and 210 gill assignments (`scripts/acute_assign.py`).
   - The 200 brain gene-control draws are `default_rng(20261005).choice(599, 200, replace=False)` over assignments 1 to 599 in that order.
5. **Sensitivities.** Each sensitivity is fitted on the observed labels only.
   - Its count of TE subfamilies at FDR 0.05 is placed against the primary null as the share of the 600 (brain) or 210 (gill) primary values that are at least as large.
   - The unique-count sensitivity uses TEcount `--mode uniq`, with rows and size factors rebuilt from the unique matrix by the same rule.
   - The label-exclusion sensitivities keep the primary size factors of the remaining libraries and re-apply the 08 filter with the subset's smaller group (6 brain, 5 gill).
   - The no-cohort sensitivity uses `~ treatment` on the primary rows.
6. **TE share.** The TE share is TE-set counts divided by all TEcount-assigned counts (genes and every repeat row) per library. It is fitted with `lm(share ~ cohort + treatment)` per tissue, and the treatment coefficient is reported in percentage points with a 95 % t interval.
7. **Cluster comparison (exploratory).** The zebrafish side is featureCounts 2.0.8 on the Tier A R1 cluster GTF (`r1_shared573.gtf`, sha256 e4831c84), with 1/NH fractions, `-O` and paired fragments.
   - A cluster is covered if it has at least 10 counts in at least 4 libraries of the tissue.
   - DESeq2 uses `~ cohort + treatment`, the tissue's primary size factors and apeglm.
   - The damselfish side follows the Tier A Q3 rule: full-cohort OA results with `strict_te`, `tested`, and the family value as the mean log2 fold change over the two parental phenotypes. A cluster takes the unweighted mean over its tested families. The comparisons are gill developmental, gill intergenerational and brain intergenerational.
   - Spearman ρ comes with a 95 % percentile interval from 5,000 bootstrap resamples of clusters (seed 20261005) and the same-sign share.
8. **Genic and intergenic reads** (only if a subfamily is called at FDR 0.05). The loci of each called subfamily come from the TE GTF. A locus is genic if it overlaps any Ensembl 112 `gene` record (exons or introns) and intergenic otherwise.
   - Counting is featureCounts 1/NH with the same settings as the cluster pass.
   - The genic share per group is the mean over its libraries of genic counts divided by genic plus intergenic counts.
9. **Report name.** The hive report is `zf-acute-20261006`, following the orchestrator's instruction (preregistration section 7 wrote `zf-acute-co2-te-<date>`). The content is unchanged.

## Downstream results (2026-10-06)

The extra count pass finished at 11:56:26Z on 26 of 26 BAMs, with no failures. Fits ran on S with the choices fixed at af6c2bd. Tables are under `tables/primary`, `tables/descriptive` and `tables/compare`.

**Primary tests (section 4).**

| Tissue | TE rows tested | At FDR 0.05 | Assignments | Permutation p | Holm p | Reading |
|---|---|---|---|---|---|---|
| Brain | 950 | 0 | 600 | 1.0 | 1.0 | UNINFORMATIVE |
| Gill | 870 | 0 | 210 | 1.0 | 1.0 | NO_DETECTABLE_TE_SHIFT |

- **Brain gene positive control: fails.** The observed fit gives 0 genes at FDR 0.05 out of 19,337. The 190th of the 200 ordered null counts is 65 (null median 1, maximum 2,019). The brain TE reading is therefore UNINFORMATIVE, as section 4 fixes. The brain TE count would also have been 0 against the null.
- **Gill genes:** 3 genes at FDR 0.05 out of 15,098. Gill has no control, as registered.
- **Null distributions:** both TE nulls are mostly 0, with a long tail (brain maximum 625, gill 142). No single library is over-represented among the assignments that give 50 or more.

**Sensitivities (section 5).** All six re-fits give 0 TE subfamilies at FDR 0.05, so each sits at the primary null's share of 1.0:
- unique-count TEcount, per tissue;
- excluding EB19 and EB21;
- excluding EG19 and EG22;
- no cohort term, per tissue.

**TE share (descriptive).**
- **Brain:** control 15.6 %, CO2 14.7 %; adjusted effect −0.78 percentage points [−3.27, 1.71].
- **Gill:** 4.6 % and 3.9 %; adjusted effect −0.88 [−2.61, 0.85].

**Genic and intergenic split:** not run, because no subfamily was called.

**Cluster comparison (exploratory, no reading).** 349 brain and 246 gill clusters are covered, and none is at FDR 0.05. The damselfish side has 433 clusters per comparison.

| Comparison | Matched | ρ (apeglm) [95 % CI] | Same sign | ρ (MLE) [95 % CI] |
|---|---|---|---|---|
| Gill, developmental | 70 | −0.05 [−0.29, 0.19] | 0.41 | −0.09 [−0.33, 0.15] |
| Gill, intergenerational | 70 | −0.17 [−0.39, 0.05] | 0.40 | −0.14 [−0.37, 0.10] |
| Brain, intergenerational | 98 | −0.22 [−0.40, −0.01] | 0.44 | −0.12 [−0.33, 0.10] |

**Note added after the fits (not a change to any registered quantity).** apeglm shrank every zebrafish effect to near zero: the largest absolute value is 3 × 10⁻⁵ for TE subfamilies and 1 × 10⁻⁵ for clusters. The MLE values reach 2.3. The signs of the shrunken values match the MLE signs, but their order rests on residues of order 10⁻⁶ (rank correlation with the MLE: 0.73 brain, 0.90 gill).

The registered apeglm ρ is reported as fixed. The MLE ρ, with a bootstrap from the same seed, is added next to it in `scripts/acute_compare.py`, and the figures plot MLE values. The brain interval that just excludes zero under apeglm does not do so under MLE. As registered, no reading attaches to either.

**Library composition (descriptive diagnostic, added after the fits).** The Ensembl mitochondrial 16S rRNA gene (ENSDARG00000080337, `Mt_rRNA`) takes a large share of gene counts in every library:
- **Brain:** 8 to 69 %. The first gene PCA axis (60 % of variance) tracks this share with Spearman 0.98.
- **Gill:** 44 to 83 %.
- **SRR28763296 (EB11, brain, CO2, cohort 2):** the highest brain share and the smallest size factor (0.20).
- **EG19 (gill):** the lowest gill share (44 %). It stands apart on gill PC1. It is one of the two gill libraries the paper's Suppl. Fig. 4 does not show.

This axis is a likely reason why the brain gene control fails, but the design cannot test that. The QC rules fixed before the data do not cover it, so no library was excluded.

**Downstream compute.**

| Step | CPU-h |
|---|---|
| featureCounts cluster pass | 0.63 |
| TEcount unique pass | 11.11 |
| DESeq2 fits (observed, 810 TE permutation fits, 200 brain gene fits, sensitivities, clusters) | 2.00 |
| Python and the PCA | under 0.01 |

On top of these, creating the R environment took a few CPU-minutes, and the cluster comparison on C took 4 CPU-s. The total for the check is about 54.5 CPU-h: pipeline 40.47, download 0.31 and downstream 13.74, against the 75 cap.

## Corrections after review (2026-10-06, append-only)

These corrections follow the numbers review of 3e73adc. No result or reading changes.

1. **Assigned fragments.** "Pipeline complete" gives TEcount assigned fragments as "24.6 to 30.9 M". The minimum in `tables/qc.tsv` is 23.99 M (SRR28763294), so the range is **24.0 to 30.9 M**.
2. **Fit count.** "Downstream results" lists the DESeq2 fits as "observed, 810 TE permutation fits, 200 brain gene fits". The itemised count is 1,023:
   - 4 observed fits (TE and gene per tissue, with apeglm on TE);
   - 810 TE permutation fits (600 brain and 210 gill), which include assignment 0, a refit of the observed labels;
   - 201 brain gene-control fits (observed labels plus the 200 draws);
   - 6 sensitivity fits;
   - 2 cluster fits.

   In `tables/primary/brain_gene_null.tsv`, data row 1 (file line 2) is assignment 0 and the 200 draws are file lines 3 to 202.
3. **apeglm bound.** The largest absolute apeglm value is 3.01 × 10⁻⁵ (gill TE), so "at most about 3 × 10⁻⁵" is the correct wording. For clusters it is 9.1 × 10⁻⁶ (brain) and 1.1 × 10⁻⁵ (gill).
