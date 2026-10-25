# Acute CO2 and TE expression in adult zebrafish brain and gill: pre-registration

**Author:** Kaisar Dauyey, Laboratory of Mathematical Biology, Hokkaido University, Japan
**Request:** the human ("green light for tier b and acidification zebrafish as acute TE activity check"), passed on by the orchestrator in message dc3322 at 2026-10-05T14:39:19Z.
**Status:** frozen before any download, at the commit that both reviewers pass. `REGISTRATION.md`, added in the next commit, states that commit and its time as the stated registration time. Nothing about these data beyond public metadata, the paper's abstract and its supplementary sample tables has been read.
**Standing:** exploratory relative to the damselfish OA and warming studies. It is a separate, small question, not part of their registered tests, and none of its outcomes changes theirs.

## 1. Question

Does five days at elevated CO2 shift transposable element (TE) expression in the brain or gill of adult zebrafish? This is an acute check: a short exposure in fresh water, one generation, no intergenerational arm. It asks whether a TE expression response to CO2 is visible at all in a fish with a well-annotated TE complement, which the damselfish studies cannot show on their own. A shift in TE RNA is not evidence of transposition, and the design cannot detect transposition.

## 2. Data

PRJNA1102381 (Sourisse et al. 2025, *Sci Total Environ* 987:179822, PMID 40472796). Adult AB zebrafish, 4 months old, held at about 500 µatm (control) or about 1,000 µatm CO2 for 5 days in three cohorts, each with both treatments. Brain 14 libraries (7 control, 7 CO2), gill 12 (6 control, 6 CO2); polyA, paired end, Illumina NovaSeq 6000, 28 to 36 M pairs each, 835 M pairs in total. Sex is not recorded.

Treatment labels are not in the public sample records. They come from the paper's supplementary tables (Suppl. Table 2, fish ID and treatment; Suppl. Table 5, sample IDs and read counts), mapped by Stanley: `samples.tsv` (copy of the hive file, sha256 `4e7d714e0c1cbdfd…`). The mapping rule is that sample EB*x* or EG*x* is fish *x*. It is supported by collection date (26 of 26), read count (26 of 26) and, for 10 gill samples, the treatment colouring of Suppl. Fig. 4 (10 of 10).

**Decided before the data:**
- **EB19 and EB21** (brain): the supplement flags brain RNA for fish 20 and 22, not 19 and 21. Fish 19 and 20 are both CO2 and fish 21 and 22 both control, so the treatment label is the same under either reading. Both are included with these labels. A sensitivity excludes both.
- **EG19 and EG22** (gill): present in the supplement's sample table but absent from Suppl. Fig. 4, with no stated reason. Both are included. A sensitivity excludes both.
- **Brain labels** rest on the naming rule plus the date and read-count checks only; no per-sample brain figure confirms them. This is stated as a limit.
- **NOT_TESTABLE** applies to a tissue if, after QC exclusions, either group has fewer than 4 libraries, or if any label is found to contradict the supplement before analysis.

## 3. Processing

The archive pipeline, as in Tier B of the zebrafish-background study: SRA download with `prefetch` (about 78 GB), fastp, STAR 2.7.11b against the GRCz11 Ensembl 112 index built at Tier B stage B2 (multi-mapping kept, as in the archive); this check waits for that index, and if Tier B stops at B0 the index is still built and charged to the Tier B budget, TEcount 2.2.4 in multi mode with the archive's gene GTF and danRer11 `rmsk` TE GTF. Counts are TE subfamilies plus genes. A second count of the same BAMs with featureCounts 2.0.8 (1/NH, `-O`) against the shared-library cluster GTF of Tier A R1 gives cluster counts for the comparison in section 5. Machine S, one third of cores free, MemAvailable ≥ 20 %, disk ≥ 15 % free.

**Genic and intergenic reads.** For every TE subfamily called at FDR 0.05, the same BAMs are counted with featureCounts against that subfamily's `rmsk` loci, split into copies overlapping an Ensembl gene body (exons or introns) and intergenic copies, and the report gives the genic share of its reads per group.

**QC exclusion, fixed now:** a library is excluded if STAR uniquely mapped reads are below 30 % or fewer than 5 M fragments are assigned by TEcount. Exclusions are reported by accession.

**Data size.** NCBI SRA lists the 26 runs at 78,456 MB in total (largest 3,388 MB), 2 × 151 bases. Peak disk on S is about 180 GB: all SRA files (78 GB), working fastq for three concurrent libraries (about 4 times one SRA file each, at most 41 GB) and the 26 sorted BAMs kept until the genic count is done (about 57 GB at the archive average of 2.2 GB). After the Tier B download S keeps about 1,890 GB free, far above the 15 % floor.

**Depth and budget.** All reads are used, with no subsampling. The orchestrator raised this check's cap from 25 to 35 CPU-h at 14:52Z on my projection of about 30 CPU-h for the pipeline. The full budget:

| Item | CPU-h |
|---|---|
| fastp, STAR, TEcount for 834.9 M pairs at 1.46 CPU-h per 40 M | 30.5 |
| featureCounts, cluster GTF, 26 BAMs at 0.033 per BAM (Tier A R1 rate) | 0.9 |
| TE-statistic permutations: 810 fits on TE rows, about 1 s each (timed on S, 12 libraries, 2,356 TE rows) | 0.3 |
| brain gene control: 200 of the 600 brain assignments drawn with seed 20261005, about 13 s each (timed, 20,095 gene rows) | 0.8 |
| genic/intergenic featureCounts for the called subfamilies | 0.4 |
| **Total** | **32.9** |

Size factors are computed once per tissue, because they do not depend on the labels. The gene control uses 200 drawn assignments rather than all 600 to stay inside the cap; its 95th percentile is then estimated from 200 values. The sensitivities in section 5 report the observed count of each re-fit against the primary null and are not permuted. After the first three libraries the projection is re-checked. The run stops at 1.5 times the pipeline estimate or at 35 CPU-h, whichever comes first, and the analysis then runs in the order of this table, primary first; past the cap it stops and asks.

## 4. Primary tests (two, one per tissue)

DESeq2 1.46.0 within each tissue, design `~ cohort + treatment`. Size factors from the full matrix (genes and TE subfamilies) after the 08 filter (≥ 10 counts in at least as many libraries as the smaller group); TE rows with the same filter; apeglm shrinkage for effect sizes; Benjamini–Hochberg within the TE rows of the tissue.

- **Statistic:** the number of TE subfamilies at FDR 0.05.
- **Null:** treatment labels permuted within cohort, enumerating every distinct assignment that keeps each cohort's group sizes (600 for brain, 210 for gill), with the statistic recomputed each time. Gill cohort 2 holds one library per group, so it contributes a single paired difference and 2 of the 210 assignments.
- **p:** the share of assignments, the observed one included, with a statistic at least as large. Holm over the two tissues.
- **Gene positive control (brain):** the number of genes at FDR 0.05 must exceed the 95th percentile of its own within-cohort null. The paper reports brain gene changes, so a failure here means the processing or labels cannot reproduce the published signal, and the brain TE reading becomes UNINFORMATIVE. It also checks the naming-inferred brain labels, since scrambled labels would most likely fail it. The paper reports no gill gene changes, so gill has no positive control. A gill TE shift without a gene shift is reported as such.

**Readings.**
- **CO2_TE_SHIFT_DETECTED:** Holm p < 0.05 (and, for brain, the gene control passes). The report gives the direction (up and down counts) and the class make-up.
- **NO_DETECTABLE_TE_SHIFT:** otherwise.
- **UNINFORMATIVE:** brain gene control fails.
- **NOT_TESTABLE:** as in section 2.

With 7 against 7 and 6 against 6 and a doubling of CO2 for five days, the test has power only for broad or strong shifts. A null reading would not show that TEs are unresponsive.

## 5. Secondary and descriptive

- **TE share of expression:** TE counts as a share of all assigned counts, per library; linear model share ~ cohort + treatment per tissue; effect with 95 % CI. Descriptive.
- **Comparison with the damselfish OA results, exploratory:** across shared clusters covered in zebrafish (≥ 10 counts in ≥ 4 libraries of the tissue), Spearman ρ between the zebrafish CO2 log2 fold change (apeglm, cluster counts) and the OA damselfish log2 fold change of the same tissue (gill: developmental and intergenerational; brain: intergenerational only, because the OA brain developmental contrast is confounded with sequencing run). ρ with a 95 % bootstrap interval and the same-sign share; no p-value, no reading. The species differ by 200 to 250 million years, in salinity and in exposure length, so agreement would be suggestive at most.
- **Sensitivities:** unique-count TEcount; excluding EB19 and EB21; excluding EG19 and EG22; dropping the cohort term.

## 6. Limits

- Two tissues, one strain, one laboratory, 4-month-old fish, five days of exposure, sexes not recorded.
- Fresh water: zebrafish acid-base regulation differs from that of reef fish in seawater, and 1,000 µatm in soft fresh water is a different chemical exposure from the damselfish treatments.
- Brain labels rest on the sample naming rule, checked indirectly by the brain gene control.
- TEcount multi mode assigns reads from TE copies inside introns and UTRs of expressed genes, so a subfamily count can move passively with its host genes. Brain gene changes are expected, so a brain CO2_TE_SHIFT_DETECTED cannot by itself separate autonomous TE transcription from host-gene read-through; the genic and intergenic shares of section 3 bound this, and the report says so where they cannot.
- TEcount multi mode distributes multi-mapping reads by an EM model; subfamily totals, not loci, are tested.
- No reading here bears on mobilisation; nothing in this design can show transposition.

## 7. Outputs

`RUN_RECORD.md` (times, CPU, exclusions, any deviation in a dated append-only note), tables, figures, and one hive report `zf-acute-co2-te-<date>`, reviewed by Phyllis (science) and Andy (numbers) before a TE_Ocean_Acidification PR. New files name machines by workstation letter.
