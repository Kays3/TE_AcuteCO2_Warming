# Pre-registration: does ocean warming change TE expression in clownfish tissues? (PRJNA970523)

Author: Kaisar Dauyey, Laboratory of Mathematical Biology, Hokkaido University, Japan. The human requested this on 3 October 2026, about 07:45 UTC. Michael approved it at 07:47:50 UTC (message a8f694, human rulings at 17:35 JST). Status: for review (Phyllis, science; Andy, numbers). Nothing has been downloaded or run. Only SRA run metadata (`runinfo`), the BioProject record, the RefSeq assembly listing and the paper's title were looked at.

**This is not an ocean-acidification result.** The stressor is temperature, in another species. It is reported as a separate study and never as evidence about CO2.

## 1. Question and what this design can and cannot show

Does rearing post-larval clownfish at 31 °C instead of 28 °C change the expression of transposable-element (TE) families in brain, liver, muscle or digestive tract? The data are bulk poly(A) RNA-seq. They can show changes in TE-derived transcript abundance at family level. They cannot show transposition, insertion or protein. All larvae come from one clutch, so the result cannot be generalised beyond that family, and family is not a confounder between temperatures.

**One tank per temperature: temperature is confounded with tank.** The paper's methods (section 2.1, read from the PDF the human supplied on 3 October 2026) state that all larvae came from a single breeding pair and that, the morning after hatching, "half of the newly hatched larvae were transferred to a separate larval tank, with one tank assigned a temperature of 28 °C and the other 31 °C". There were no replicate tanks. The 6 larvae per group are replicates within one tank, not replicate tanks. Every result is therefore phrased as a difference between the 28 °C tank and the 31 °C tank, which cannot separate temperature from any other tank effect, and the permutation p is read as a test of that tank difference.

## 2. Data

- **BioProject PRJNA970523**, *Amphiprion ocellaris*. Paper: "Tissue-specific transcriptional response of post-larval clownfish to ocean warming", Sci Total Environ (PMID 37923256).
- **Included:** the 48 runs with `LibrarySelection = PolyA`. These are 4 tissues × 2 temperatures × 6 libraries, all from clutch C5, 20 days post-hatch. **The 48 libraries come from 12 larvae, 6 per temperature, each giving all four tissues.** The paper (section 2.2) states that "six individuals were randomly sampled from each of the two treatments" and that "the brain, liver, muscle, and entire digestive tract of each larvae" were extracted, giving 48 tissues. The BioSample `isolate` attribute (Larva1–Larva48, numbered in tissue blocks) contradicts this and is treated as a submission artefact. The `L1`–`L10` labels in the sample names do not recover the fish either: at 28 °C nine different labels occur, and no label set gives six fish each present in all four tissues. Libraries therefore cannot be mapped to fish. Within a tissue, the 12 libraries are 12 different larvae, so the per-tissue tests are unaffected. Across tissues, the four tests use the same 12 larvae and are not independent; the Bonferroni split (0.05 / 4) remains valid under that dependence.
- **Excluded:** the 12 runs with `LibrarySelection = Inverse rRNA` (digestive-tract metatranscriptomes). They use a different library chemistry, and their source is metatranscriptomic.
- **Runs:** NovaSeq 6000, paired 2 × 150, 27.6–52.4 M pairs per run, 162.5 GiB in SRA (48 included runs; all 60 runs total 211.0 GB). The full list is `PRJNA970523_runinfo.csv`, committed with this file.
- **Reference:** GCF_022539595.1 (ASM2253959v1), with the RefSeq annotation GTF (release 2023_03). The sha256 of each file is recorded before use.

## 3. TE annotation (the main methodological choice)

No repeat annotation is distributed for this assembly, and Dfam has no curated *Amphiprion* families.

- **Library.** The A. polyacanthus de novo library `apoly-families.fa` (sha256 152a2a7e…, the OA library, same family Pomacentridae), plus Dfam curated families for Actinopterygii. The Dfam families are extracted with `famdb.py` from the `dfam/tetools:2.00` container on pSS (image a6179a355ba7). The container's RepeatMasker and Dfam versions are recorded before use. Family names are prefixed `apoly:` or `dfam:`.
- **Masking.** RepeatMasker with the combined library, `-s -a -nolow -no_is -gff`, on the soft-unmasked genome, 16 threads.
- **TE GTF.** The output is converted to TEcount format by a committed script:
  - one record per copy;
  - `gene_id` = family, `transcript_id` = copy, plus `family_id` and `class_id`;
  - same layout as the OA TE GTF (sha256 71bb6589…);
  - simple repeats, low complexity, rRNA, tRNA, snRNA and satellites are excluded, as in OA.
- **Adequacy check, fixed now.** The interspersed-repeat masked fraction of the genome is reported. If it is below 20 percent, I stop before any alignment and ask Michael. RepeatModeler2 needs a fresh sign-off.
- **Stated bias.** Young, lineage-specific *Amphiprion* families will be missed or merged into their closest pomacentrid or Dfam relative. The report states this next to every family-level result. Because the shared `apoly:` family IDs are used in both studies, the comparison with OA (section 5.4) is direct for those families.

## 4. Processing pipeline (frozen, as in the OA cohort)

- **Download:** SRA Toolkit `prefetch` and `fasterq-dump`, with `vdb-validate` on each `.sra`. Two runs at a time. Inbox re-checked before the download starts.
- **Trimming:** fastp 1.3.6 `--detect_adapter_for_pe`.
- **Alignment:** STAR 2.7.11b, `--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD`. The index is built with `--sjdbOverhang 149`, matched to the 150 bp reads. That is the one intended difference from OA, which used 100 for 101 bp reads.
- **Counting:** TEcount 2.2.4, `--mode multi`, with the RefSeq gene GTF and the section 3 TE GTF. samtools 1.23.1.
- **Strandedness: tested once, then fixed.** The metadata do not state strandedness. On one brain library, a 2 M-pair subset (seqtk, seed 20261003) is aligned with `--quantMode GeneCounts`, and the column 3 / column 4 ratio of `ReadsPerGene.out.tab` sets one `--stranded` value for all 48 libraries:
  - ratio > 0.8 or < 0.2 (column 3 / (column 3 + column 4)): `forward` or `reverse`;
  - otherwise `no`.
  This is decided before any full library is counted and recorded in the run record.
- **Environment:** a new micromamba environment `jim_te_warm` on pSS (STAR, fastp, samtools, sra-tools, TEtranscripts, seqtk, DESeq2 1.46.0 on R 4.4.3). The full export goes to the run record before the first computation. If TEcount 2.2.4 cannot be installed, the version used is fixed by a dated amendment before any count is made.
- **Exclusions, fixed now:** a library is excluded only if its STAR uniquely-mapped rate is below 50 percent or TEcount fails. No other exclusion is allowed.

## 5. Statistics

### 5.1 Primary: TE family differential expression, per tissue

- **Model:** DESeq2 1.46.0 on the joint gene + TE-family count table of each tissue, with design `~ temperature` (31 vs 28, 6 vs 6).
  - Size factors: `poscounts`.
  - Independent filtering is disabled.
  - A feature is kept if it has ≥ 10 counts in ≥ 6 libraries of that tissue.
- **FDR:** Benjamini–Hochberg over TE families only, within a tissue.
- **Outcome names (per tissue):**
  - `WARMING_SHIFT_DETECTED`: at least one TE family has FDR < 0.05. Each such family is reported with direction, log2 fold change and its 95 percent CI.
  - `NO_DETECTABLE_SHIFT`: no family passes. The CI of every family is still reported. Equivalence is not tested, and "no change" is never written.
- **Exact permutation null.** Within each tissue, all C(12, 6) = 924 assignments of the temperature label to the 12 larvae of that tissue are enumerated; the same 12 larvae appear in every tissue, but each tissue is tested on its own. For each labeling, the number of TE families at FDR < 0.05 is recorded. The permutation p is the rank of the observed count among the 924 labelings, the observed labeling included: p = (number of labelings with a count ≥ the observed) / 924. No seed is needed because the enumeration is complete. The labeling and its mirror image (28 ↔ 31) give the same count, so the null has 462 distinct values at most; this is stated with the result.

**Null after an exclusion.** If section 4 excludes libraries, the null enumerates all labelings of the remaining n larvae with the observed group sizes n28 and n31: p = (labelings with count ≥ observed) / C(n, n28). The minimum p is 2/C(n, n28) when the groups are equal (mirror tie) and 1/C(n, n28) otherwise. A tissue whose minimum reachable p is ≥ 0.0125 is reported as `NOT_TESTABLE` at the primary threshold, with its DESeq2 results shown descriptively. Minima: 6 v 5, 1/462 = 0.0022; 5 v 5, 2/252 = 0.0079; 6 v 4, 1/210 = 0.0048; 5 v 4, 1/126 = 0.0079; 6 v 3, 1/84 = 0.0119 (reachable). Untestable: 4 v 4 (2/70 = 0.029), 5 v 3 (1/56 = 0.018), 6 v 2 (1/28 = 0.036).

**One procedure for observed and permuted fits.** Every labeling, the observed one included, is fitted identically: the full joint gene + TE table of the tissue, the same feature filter, and the same size factors. The filter (≥ 10 counts in ≥ 6 libraries) and the `poscounts` size factors do not depend on the labels, so both are computed once per tissue on all 12 libraries and held fixed for all 924 fits. Dispersions are re-estimated in each fit, as they are in the observed fit. No reduced feature set is used anywhere, so the gene positive control (5.3.1) and the matched-decile statistic (5.2.3) are taken from the same 924 fits.

**Four tissues = four primary tests.** A tissue counts as positive only if both FDR < 0.05 for ≥ 1 family and permutation p < 0.05 / 4 = 0.0125. Because the observed labeling's mirror is among the 924 and always ties it, the smallest reachable p is 2/924 ≈ 0.0022, so the threshold is reachable.

### 5.2 Secondary

1. **Direction among passing families:** a two-sided binomial sign test of the up/down split. If fewer than 10 families pass, it is reported with an explicit underpowered flag and nothing further is built on it.
2. **Global TE index:** strict-TE counts as a percent of all counts per library. This uses the implementation of `paper/results/external-co2-stress-validation-20260907/scripts/build_index_and_matrix.py`, unchanged, with its `is_strict_te` rule adapted only for the `apoly:`/`dfam:` prefixes. Per tissue, the difference in means with 95 percent CI, Hedges' g, and Welch's t (two-sided) are reported, with Mann–Whitney U alongside. Absolute index values are not compared with OA or other cohorts.
3. **TE versus gene shift:** the median absolute log2 fold change of TE families against genes of matched base mean (10 baseMean deciles). This is descriptive and is reported with the permutation distribution of the same statistic.
4. **Tissue contrasts:** `~ tissue + temperature` over all 48 libraries. The fish term that the design calls for cannot be fitted, because libraries cannot be mapped to fish, so the model treats 48 libraries from 12 larvae as independent. Its p-values are therefore too small and are not reported; only effect sizes are shown. This is descriptive and is used only to show that TE expression is measurable and tissue-structured in this dataset.

### 5.3 Controls

1. **Pipeline positive control.** The number of genes at FDR < 0.05 per tissue is compared with its own permutation distribution.
   - If no tissue exceeds the 95th percentile of that distribution, TE nulls are reported as uninformative (low power or a weak warming response), not as evidence of no TE effect.
   - The gene result is also compared descriptively with the paper's DEG counts (brain 102, liver 1,785, digestive tract 380, muscle 447; liver strongest, brain weakest) and its 23 heat-shock protein genes up in all tissues. The pipelines differ, so only the ranking and the HSP direction are compared.
2. **Unique-only sensitivity:** TEcount `--mode uniq`, same models. A family that passes only in multi mode is flagged.
3. **Young families:** families whose median Kimura distance (CpG-adjusted, from the RepeatMasker `.align`) is below 5 percent are reported separately. This is descriptive, because they are also the families the cross-species library serves worst.

### 5.4 Comparison with OA (descriptive)

For `apoly:` families tested in both studies, the liver warming log2 fold changes (this study) are set against the OA liver CO2 fold changes (developmental and intergenerational, merged OA results). The comparison reports sign concordance and Spearman correlation. It is descriptive only; species, stressor and age differ.

## 6. Compute, host and stop rules

- **Host:** pSS (32 cores, x86_64), kaisar account, no sudo. Every process is pinned with `taskset -c 0-20`, so at least a third of the cores stay free. Every process runs under `/usr/bin/time -v`.
- **Estimate:**
  - RepeatMasker: 20–40 CPU-h.
  - STAR index: 2–3 CPU-h.
  - Alignment: 40–55 CPU-h.
  - TEcount: 10–20 CPU-h.
  - DESeq2 and permutations: 924 labelings × 4 tissues = 3,696 full fits on the filtered joint table (section 5.1), about 10–25 CPU-h at 10–25 s per fit. The per-fit time is measured on the observed fits before the permutations start.
  - **Total about 85–145 CPU-h.** Michael approved 75–120 CPU-h, set before the full-fit permutation was fixed. The stop rule below therefore covers the whole projection, DESeq2 included.
- **Pilot:** two libraries are processed end to end first. If their timing, plus the measured DESeq2 per-fit time once the observed fits exist, projects the total past 120 CPU-h, I stop and ask Michael before going further.
- **Disk:** about 1.3 TB free on pSS `/`. FASTQs are deleted after each library is counted. BAMs are kept until the report is merged.

## 7. Outputs

Code, this pre-registration, `RUN_RECORD.md` (hashes, versions, CPU-h, deviations) and compact tables go in `paper/results/te-warming-aocellaris-20261003/` on branch `analysis/te-warming-aocellaris-20261003` of the private TE_Ocean_Acidification repository. The report goes to `hive/reports/{md,html,pdf}/te-warming-aocellaris-20261003.*`, with figures rather than large tables. A PR follows: Phyllis reviews the science, Andy the numbers, and Michael merges. Any change after review is a dated, append-only amendment below, made before the stage it affects.

## 8. Amendments

### Amendment 1 (3 October 2026, ~09:50 UTC): species TE reference, environment, outcome names, and three extensions

Registered at Michael's request (message ad5d36, human at 18:35 JST: "let Jim conduct more work on warming"), after Phyllis's science PASS at da4d4a8 and before any download or computation. It supersedes the base text where they conflict. Nothing has been downloaded or run.

**A1.1 TE reference: the species' own library replaces the proxy (supersedes section 3, Library and Masking).** The Mobilome paper's comparative panel (TE_Mobilome_Paper, decision register D9/D10) already ran RepeatModeler 2.0.9 and RepeatMasker 4.2.4 (`dfam/tetools:2.00`, pSS) on this exact assembly, GCF_022539595.1. The library has 1,849 families (`families.fa`, sha256 16f27731…, stored under the run's internal database name `apoly-families.fa`); the RepeatMasker `genome.fa.out` sha256 is 9f51c6e7…, which matches the committed provenance file in the Mobilome repository. RepeatMasker masked 35.21 % of 863.5 Mb. The interspersed classes (DNA 12.79 %, Unknown 10.97 %, LINE 6.75 %, LTR 2.90 %, RC 0.19 %, SINE 0.19 %, PLE 0.11 %, Retroposon 0.03 %; class-wise merged sums) total about 33.9 %, above the 20 % stop line, so the adequacy check passes on existing output. The published value for this assembly is 44.7 % (Ryu et al. 2022, G3), with a different pipeline; the difference is stated, not corrected. The TE GTF is built from this `.out` by the committed converter in section 3 with the same exclusions; `Unknown` families stay in the GTF as features but are not strict TE for the index, as in OA. Consequences: the RepeatMasker step (20–40 CPU-h) is dropped; the young-*Amphiprion* bias in section 3 no longer applies in that form; family IDs are `rnd-*_family-*` of this library and are no longer shared with the OA library, so section 5.4 uses the mapping in A1.6.

**A1.2 Environment (supersedes section 4, Environment).** The frozen OA environment already on pSS (`envs/te-ocean-acidification`) carries the exact pinned versions: STAR 2.7.11b, fastp 1.3.6, samtools 1.23.1, TEcount 2.2.4, DESeq2 1.46.0 on R 4.4.3, TEsorter. It is used read-only by absolute path. A small new environment `jim_te_warm` adds only sra-tools, seqtk and Subread 2.0.8 (featureCounts, as in the mobilisation study). Both exports go to the run record before the first computation.

**A1.3 Outcome names (Phyllis, item 1).** `WARMING_SHIFT_DETECTED` becomes `TANK_TEMPERATURE_SHIFT_DETECTED`, and `NO_DETECTABLE_SHIFT` becomes `NO_DETECTABLE_TANK_SHIFT`, because one tank per temperature means the contrast cannot isolate temperature. A recovered HSP induction (5.3.1) shows that the tank difference includes a heat response; it cannot exclude other tank effects, and the report states both.

**A1.4 design.tsv (Phyllis, item 2).** The column `isolate` is renamed `biosample_isolate_ARTEFACT_not_fish_id`. Its values are the submitter's BioSample field; 12 larvae gave all 48 libraries (section 2), and libraries cannot be mapped to fish.

**A1.5 Extension E1: mobilisation framework, 28 °C vs 31 °C tank (secondary).** The mobilisation study (`paper/results/te-mobilisation-transcription-20261002/`, PREREGISTRATION sections 3–4 and A1–A3) is applied unchanged except for inputs: its scripts, tier rules (T1 full length, Kimura ≤ 2 % CpG-adjusted, intact ORF ≥ 95 % of the consensus longest ORF, LTR pair rule; T2; SINE-T; T1u; C1 old and fragmented), per-class ORF thresholds, detection rule (≥ 5 unique fragments in ≥ 3 libraries of a tissue), self-driven / passenger / undetermined classes at 3× (2× and 5× sensitivity), and 99 + 900 shuffles (seed 20261003). Inputs are the A1.1 library, `.out` and `.align`, the RefSeq GTF, and the 48 STAR BAMs of section 4. Differences forced by the data: the libraries are 150 bp (OA 101 bp), and libraries cannot be mapped to fish, so the tissue contrast `~ fish + tissue` is not fitted and tissue patterns are descriptive. The warming test is DESeq2 `~ temperature` per tissue on T1/T1u/SINE-T/T2 family sums (unique and rounded 1/NH counts; size factors from gene counts as in section 5.1), BH within tissue, with the exact labeling null of section 5.1 (same function for observed and permuted labelings, same post-exclusion rule). Outcomes per tissue: `TANK_TEMPERATURE_SHIFT_INTACT` (at least one family at FDR < 0.05 and permutation p < 0.0125, reported with its family outcome, for example SELF_TRANSCRIBED_INTACT) or `NO_DETECTABLE_INTACT_SHIFT`. Family outcomes (SELF_TRANSCRIBED_INTACT, TRANSCRIBED_PASSENGER_ONLY, INTACT_NOT_DETECTED, NOT_RESOLVABLE) are reported per temperature group descriptively. As in the mobilisation study, a self-driven class means coverage above the flanks in unstranded short reads, not a demonstrated promoter, and nothing here shows transposition. Compute: Stage 1 about 2 CPU-h; Stage 2 counting about 5–10 CPU-h (the mobilisation run took 8.5 CPU-h for 162 shorter libraries); Stages 3–4 with permutations about 1–3 CPU-h. **E1 total about 8–15 CPU-h.**

**A1.6 Extension E2: warming versus CO2 on matched families (descriptive; supersedes section 5.4).** Families are matched between this library and the OA *A. polyacanthus* library (sha256 152a2a7e…) by the 80-80-80 rule on consensus sequences (≥ 80 % identity over ≥ 80 % of both lengths, ≥ 80 bp), reciprocal best hit, rmblastn from the tetools container. Only one-to-one matches are compared. For brain and liver, the two tissues present in both studies, warming log2 fold changes (this study) are set against the OA CO2 log2 fold changes (developmental and intergenerational, merged OA results): sign concordance and Spearman correlation, with the share of families matched reported first. Muscle and digestive tract have no OA counterpart and are not compared. It is descriptive only; species, stressor, age, rearing design and library differ. Compute under 0.5 CPU-h.

**A1.7 Extension E3: superfamily-level contrast (secondary).** Family counts are summed to RepeatMasker superfamily (for example DNA/TcMar-Tc1, LINE/L2) and tested per tissue with the section 5.1 model, BH within tissue, same labeling null. This trades resolution for power at n = 6 vs 6 and is reported next to the family result, never in place of it. Compute under 1 CPU-h.

**A1.8 Compute after this amendment.** Base items without RepeatMasker: about 62–105 CPU-h (STAR index 2–3, alignment 40–55, TEcount 10–20, DESeq2 and permutations 10–25). Extensions E1–E3: about 9–16 CPU-h. **Total about 71–121 CPU-h**, inside the approved 75–120 at the central estimate. The pilot stop rule is unchanged: if the measured projection, extensions included, passes 120 CPU-h, I stop and ask Michael before going further.

### Amendment 1 addendum (3 October 2026, ~10:00 UTC): Unknown families, the .align file, file identity

Phyllis's science review of Amendment 1 (R1–R3) and Andy's numbers notes. Registered before any RepeatMasker-derived file, index or count is used.

**A1.9 Unknown families: one rule for 5.1, E1 and E3.** The primary BH set, and the set that decides `TANK_TEMPERATURE_SHIFT_DETECTED`, is classified families only: DNA, LINE, SINE, LTR, RC, PLE and Retroposon. Families of class `Unknown` are tested in the same per-tissue model but in a separate BH set labelled "unclassified repeats", reported next to the primary and never counted toward any outcome. Any Unknown family with FDR < 0.05 is reported with the share of its copies that overlap RefSeq exons, so that host-gene fragments and segmental duplications are visible. The same split applies to E3 (superfamily) and to E1, where T1u stays its own stratum and does not enter the E1 outcome. The permutation null of 5.1 counts the primary set only.

**A1.10 The .align file.** The Mobilome run wrote `genome.fa.align` (sha256 9d0105e5…, 2.64 GB) and `genome.fa.out.gff` (sha256 33fe5dfd…) next to the `.out` (sha256 9f51c6e7…), so it was run with `-a` and `-gff`. Its records carry `Kimura (with divCpGMod)`, which 5.3.3 and the E1 young-copy rule use. The exact RepeatMasker command line is not recorded in the Mobilome repository. From the outputs: `-a` and `-gff` were set; simple repeats are present, so `-nolow` was not set; `-s` and `-no_is` cannot be confirmed. This departs from the section 3 options and is accepted: simple and low-complexity repeats are removed when the TE GTF is built, and the `-s` setting changes sensitivity at the margin, not family identity. No RepeatMasker rerun is made.

**A1.11 File identity.** Libraries are identified by sha256, never by file name. The species library is copied to `~/jim_te_warm/ref/aocel-families.fa` (sha256 16f2773127d2b942017a93a6e3cb19e8fdcd27f2e19c9ac05c1984cfb654e038, 1,849 sequences). The OA *A. polyacanthus* library is read only as sha256 152a2a7e…. Every script takes each library, `.out` and `.align` by absolute path, checks its sha256 on start, and exits on a mismatch. These hashes go to `RUN_RECORD.md`; the library hash and family count are not in the Mobilome provenance file and are recorded here for the first time.

**A1.12 Section 5.3.3 rationale.** With the species library, young families are no longer "the families the cross-species library serves worst". They are reported separately because near-identical young copies collect fewer unique reads, so their counts rest more on multi-mapping assignment.

**A1.13 E2 reporting.** Spearman rho is reported with its 95 % CI and n; any p-value is labelled descriptive.

### Amendment 2 (3 October 2026, ~09:55 UTC): compute ceiling raised to 200 CPU-h

The human approved a ceiling of 200 CPU-h on pSS for the corrected base study and the Amendment 1 extensions (Michael, message e3fb7e, 18:37 JST). The stop-and-ask threshold in section 6 and A1.8 moves from 120 to 200 CPU-h. Projected use by stage: download 1–3 (I/O bound, measured); STAR index 2–3; alignment 40–55; TEcount 10–20; DESeq2 and permutations 10–25; E1 mobilisation 8–15; E2 under 0.5; E3 under 1. RepeatMasker is 0, because the existing output is reused (A1.1). **Total about 72–123 CPU-h against 200.** Every process runs under `/usr/bin/time -v`, and cumulative CPU-h against 200 is written to `RUN_RECORD.md` at the end of each stage. If the measured projection passes 200, I stop and ask Michael.

### Amendment 3 (3 October 2026, ~09:55 UTC): TE GTF keeps all repeat classes, as OA actually did

Section 3 says simple repeats, low complexity, rRNA, tRNA, snRNA and satellites are excluded from the TE GTF "as in OA". That is wrong about OA. The OA TE GTF (`apoly.de-novo-te.gtf`, sha256 71bb6589…) keeps every RepeatMasker class (checked on pSS: Simple_repeat 266,254 loci, Low_complexity 35,320, Satellite 2,830, rRNA 1,438, tRNA 178, beside the TE classes), and the strict-TE set is chosen at analysis time; the Mobilome converter `pipelines/repeat-annotation/repeatmasker-out-to-gtf.py` has no class filter. To keep the counting step identical to the frozen OA pipeline, the TE GTF here is built by that converter, unchanged, from the A1.1 `.out`, and keeps all classes. Which features enter which test is already fixed at analysis time: the primary BH set and every outcome use the classified TE classes of A1.9; Unknown families form the separate labelled set; non-TE classes enter no TE test or index and are kept only as count features. Registered before the GTF is built or any read is counted.

**Amendment 3 addendum (3 October 2026, ~10:00 UTC; Phyllis and Andy, before any DESeq2 fit).**

- *DESeq2 input rows, following OA.* The OA full-cohort fit (`pipelines/full-cohort/deseq-te-full-cohort.R`) drops only class `rRNA` rows before fitting, because ribosomal sequence captured by the annotation carried about 6.6 % of TE read mass in the OA pilot; all other rows (genes, TE classes, Unknown, Simple_repeat, Low_complexity, Satellite, tRNA, snRNA) stay in the joint table, enter size factors and the dispersion trend, and only strict-TE rows enter the TE results. The same rule is used here, for the observed fit, all permuted fits, E1 family sums and E3: `rRNA` rows are dropped before fitting; every other row stays in the table; only the A1.9 classified set enters the primary BH set and the outcomes, and Unknown forms its separate labelled set. Non-TE repeat rows enter no BH set, test or index. The rRNA read share per library is reported.
- *A1.10 rationale, corrected.* A1.10 accepted the missing `-nolow` because simple and low-complexity repeats were said to be removed at the GTF build. Under Amendment 3 they are not removed there; they are kept as count features and excluded at analysis time. The acceptance still holds: the OA annotation also contains simple and low-complexity repeats, so the two counting steps match, and no test uses those rows.

### Amendment 4 (4 October 2026, ~08:00 UTC): compute ceiling 250 CPU-h, a parallelism rule and a recompute rule

**Why.** Counting finished on 4 October at 04:02:29 UTC at 202.26 CPU-h (registered per-process measure, `RUN_RECORD.md`), which reached the 200 CPU-h ceiling of Amendment 2 before any analysis step. The cause was oversubscription. Six libraries ran at once inside a 21-core pin, each with multi-threaded STAR, and STAR cost 9,124 CPU-s per library against 2,731 in the pilot. The stop rule bound, and nothing has run since. The human answered at 16:50 JST (07:50:45 UTC, 4 October), verbatim: "increase the ceiling to reach the capacity". Michael's reading, stated to the human and relayed in message c991ae: option (a), a ceiling of 250 CPU-h, enough for the registered work, with host guards unchanged. This amendment uses that reading. Michael's message calls it "Amendment 3"; it is numbered 4 here because Amendment 3 (TE GTF classes) already exists.

**A4.0 Deviation recorded.** The overrun to 202.26 CPU-h is a dated deviation in `RUN_RECORD.md` (2026-10-04, counting stage closed). It records the cause, oversubscription of multi-threaded STAR (9,124 against 2,731 CPU-s per library), and the missed recompute at the 17:31 UTC check. The registered CPU measure is the per-process `/usr/bin/time -v` sum. The per-library wrapper sum (184.87) is recorded for comparison only.

**A4.1 Ceiling.** The stop-and-ask threshold in section 6, A1.8 and Amendment 2 moves from 200 to 250 CPU-h. Every CPU-h spent so far counts toward it. Projection of the remaining registered work, with estimates unchanged from Amendment 2:
- 5.1–5.3 DESeq2 with 3,696 full fits: 10–25 CPU-h.
- E1 mobilisation: 8–15 CPU-h.
- E3: under 1 CPU-h.
- E2: under 0.5 CPU-h.

That gives about 222–244 CPU-h in total. Nothing in the analysis changes: outcome names, the exact 924-labeling null, the feature filter, size factors, BH sets (A1.9), the gene positive control against PMID 37923256, E1–E3 and the order of stages all stay as frozen.

**A4.2 Parallelism rule.**
- DESeq2 fits run single-threaded: `BiocParallel::SerialParam()`, with `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` set.
- At most 10 fits run at once, all inside `taskset -c 0-20`.
- For any multi-threaded step (E1 featureCounts), threads per job × concurrent jobs must not exceed the number of pinned cores.

**A4.3 Recompute rule.**
- At every status check and at the end of each stage, the total is projected again from measured per-process times: spent so far, plus measured cost per unit times the units remaining. The projection is written to `RUN_RECORD.md` with its UTC time.
- The per-fit DESeq2 time is measured on the 4 observed fits before the permutations start, as section 6 already requires. The projection is redone from it and sent to Michael. No permuted fit runs until Phyllis and Andy have passed this amendment.
- If any projection passes 250 CPU-h, work stops at the end of the running unit (a fit batch or a library), and I ask Michael before going further.

**A4.4 Order under the ceiling.** If a projection made under A4.3 shows that the ceiling cannot cover everything, the registered order decides what runs first:
1. Section 5.1, with its permutation null and the gene positive control (5.3.1).
2. Then 5.2–5.4, E3 and E2.
3. Then E1.

Nothing is dropped silently. Whatever does not run is reported as not run, with the reason.

**A4.5 Host guards.** These are unchanged:
- one third of pSS cores stay free (`taskset -c 0-20`);
- MemAvailable stays above 10 % until 2026-10-05 05:31 UTC and above 20 % after;
- `/` keeps at least 15 % and 300 GB free;
- nothing is deleted (BAMs stay) until Michael gives the word after the merge.

**Amendment 4 addendum (4 October 2026, by 08:05 UTC; Phyllis R1–R5, Andy items 1–3, parallel plan at the human's request). Append-only. Registered before any permuted fit.**

- *Erratum, time stamps (Phyllis R1, Andy 3).* Two stated times were later than their commits:
  - The Amendment 4 heading reads "~08:00 UTC"; the commit, e76f808, is from 07:51:39 UTC. Read it as "by 07:51 UTC".
  - The `RUN_RECORD.md` heading "Counting stage closed (recorded 2026-10-04 07:50 UTC)" was committed at 07:40:11 UTC. Read it as "by 07:40 UTC". It is corrected in `RUN_RECORD.md` by an appended note.
- *Erratum, wording of A4.0 and "Why" (Phyllis R2).* Read it as: "Counting finished at 202.26 CPU-h, which passed the 200 CPU-h ceiling of Amendment 2 by 2.26 CPU-h before any analysis step; the stop rule was applied only after the overrun, and nothing has run since." The missed check was at 17:31 UTC on 3 October.
- *Erratum, the missed-check figure (Andy 2).* The `RUN_RECORD.md` sentence "the run could have been stopped at about 120 CPU-h spent" does not follow from the record.
  - At 17:31 UTC on 3 October, 14 libraries were finished.
  - Spend was then about 21.9 CPU-h for the first 8, plus 6 × about 4.5 CPU-h. That is the mean of the last 40 bulk libraries: (200.50 − 4.08 − 15.63) / 40 = 4.52.
  - Six libraries were in flight.
  - So the run could have been stopped at roughly 60–70 CPU-h spent. That is the figure. `RUN_RECORD.md` gets the same correction.
- *Erratum, A4.1 projection (Phyllis R3, Andy 1).* 202.26 + (10 to 25) + (8 to 15) + (under 1) + (under 0.5) = about 220–244 CPU-h, not 222–244. After the observed timing fits, the measured projection in `RUN_RECORD.md` is about 227–234.
- *A4.4 reworded (Phyllis R4).* Replace "the registered order decides what runs first" with "the following priority, fixed here before any analysis output exists, decides what runs first". Add: whether a stage stops or is skipped is decided by the CPU projection alone, never by any result.
- *A4.3 extended (Phyllis R5): observed fits.*
  - The 4 observed fits use exactly the A4.2 settings and the same scoring function as the permuted fits. They ran at 07:53:47–07:54:02 UTC with `SerialParam`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` set to 1, and the `fit_one()` function of `scripts/deseq_fit.R` (index 0). The permuted fits use the same function.
  - Their outputs are frozen by sha256 in `RUN_RECORD.md` and are not reused or edited.
  - **Disclosure.** Collecting the timing printed the four per-tissue summary rows to my session. These are the counts of FDR < 0.05 features per set. So the observed statistics have been seen by me, before the null exists. No full results table was opened.
  - That exposure cannot steer a choice under A4.4. The measured projection (about 227–234) is below 250, so every registered stage runs and nothing is dropped. If a later projection forced a drop, the A4.4 priority applies mechanically, on the CPU projection alone.
- *A4.6 Parallel plan (human, 07:54 UTC: "propoze parallelization and use more machines to make the work efficient if possible").*
  - **Host.** All DESeq2 fits, observed and permuted, run on pSS in the one frozen environment, so the observed and null fits come from one platform and one library build. No other host is used: the fits total about 8 CPU-h, and the E1 counting reads the 144 GB of BAMs that sit on pSS.
  - **Slots.** Inside `taskset -c 0-20` (21 cores), at most 20 compute threads run at once. A single-threaded R or Python process counts as 1. A multi-threaded job counts as its thread count; for example, featureCounts `-T 6` counts as 6. The remaining core and the 11 cores outside the pin stay free.
  - **Concurrency.**
    - E1 Stage 1 (tier calling from the library, `.out`, `.align` and genome) and E1 Stage 2 (featureCounts over the 48 BAMs) do not depend on any DESeq2 output.
    - E3 (superfamily sums of the same count tables) does not depend on the 5.1 permutations.
    - So 5.1 permutations, E3 fits (observed and permuted) and E1 Stages 1–2 run at the same time, sharing the 20 slots.
    - While one E1 featureCounts job (`-T 6`) runs, at most 14 single-threaded fits run. Otherwise up to 20 do.
    - E1 Stages 3–4 follow Stage 2. E2 and the 5.2 and 5.3 summaries follow the 5.1 null.
  - **Unchanged.** Every registered estimate, the 250 CPU-h ceiling, the one-third-free guard, the memory guard and the A4.4 priority (decided by CPU projection only).

### Amendment 5 (6 October 2026; first version committed 13:56:44 UTC at c9f90b0, revised after review): POST-HOC re-evaluation with an rRNA correction and a relaxed gene control

**Status.** This amendment is POST-HOC. The study is complete, and its registered readings were reported and merged (TE_OA PR #181, errata #182).

- **Registered readings stand.** The registered readings of sections 5.1 to 5.4 and Amendments 1 to 4 are not changed, re-labelled or replaced.
- **New readings sit beside them.** Every new reading below carries the prefix `POSTHOC_` and is reported next to its registered counterpart.
- **Not a test of the original hypothesis.** It cannot rescue or overturn any registered result. A `POSTHOC_TANK_TEMPERATURE_SHIFT_DETECTED` reading is exploratory. It was chosen after the registered results and the acute finding were seen, and it is never cited as a confirmatory test or a replication. In any manuscript it appears only with the POSTHOC label and the tank clause.

**Request.** The human, on 6 October 2026 at about 22:50 JST, relayed by the orchestrator (message e69a82), verbatim: "for warming check muscle again for positive control and relax the criteria, reevaluate for all tissues based on correction of mt 16S RNA of gene counts".

**What motivated it (disclosed).** Two earlier results, both already seen by me, prompted this re-evaluation:
1. **Muscle control.** In the registered analysis the muscle gene positive control (5.3.1) did not exceed the 95th percentile of its null: 101 genes at FDR < 0.05 against a null 95th percentile of 402. Brain, digestive tract and liver exceeded theirs.
2. **Acute zebrafish CO2 check** (TE_OA PR #188). There, mitochondrial 16S rRNA took 8 to 69 % of brain and 44 to 83 % of gill gene counts, and in brain it was the main axis of gene-level variation. Under that check's frozen design, the failed brain gene control left brain without a TE reading.

**What I have not seen.** I have not looked at the mitochondrial rRNA share of any warming library. I have not looked at any percentile of the registered gene-count nulls other than the 95th already reported. I have not looked at any count, fit or permutation under the corrections below. Everything in this amendment was fixed before those numbers were opened, at c9f90b0.

**Seen during review (disclosed).** Andy's numbers review of c9f90b0 (message af5105, 13:58 UTC) reported the rank-832 values of the stored registered gene nulls to me: brain 57, liver 73, muscle 92, digestive tract 43, against observed counts of 917, 1,436, 101 and 1,017. I had not opened them. The rule of A5.5 is the one committed at c9f90b0 before that message, and this revision does not change it. No M1 or M2 number existed or had been opened.

**A5.1 rRNA features.** In the RefSeq GTF used for counting (GCF_022539595.1), the rRNA features are every gene with `gene_biotype "rRNA"`. That is 95 genes:
- the two mitochondrial rRNA genes on NC_009065.1: `KEG47_r01` (16S, positions 1,090 to 2,786) and `KEG47_r02` (12S, positions 70 to 1,017);
- 93 nuclear 5S, 5.8S, 18S and 28S rRNA genes.

This list comes from the annotation only. RepeatMasker class `rRNA` rows are already dropped (Amendment 3 addendum).

**A5.2 Measured share (descriptive).** For each of the 48 libraries, from the TEcount multi table:
- the mitochondrial rRNA share, defined as the counts of `KEG47_r01` plus `KEG47_r02` divided by all gene counts (all non-repeat rows);
- the 16S share alone;
- the all-rRNA share.

The shares are reported per library, with the per-tissue means of the 28 °C and 31 °C groups. No test is run on them.

**A5.3 Correction M1, rRNA genes removed (the POST-HOC primary).** This is the registered procedure of 5.1 and 5.3.1, unchanged except for the input rows. The 95 rRNA gene rows of A5.1 are dropped from the joint table before the feature filter, the `poscounts` size factors and every fit, alongside the RepeatMasker `rRNA` rows already dropped.

Everything else is as registered:
- the filter (≥ 10 counts in ≥ 6 libraries);
- design `~ temperature`;
- independent filtering off;
- BH over the A1.9 classified TE families;
- the full 924-labeling null computed per tissue with the same function for observed and permuted labelings;
- the threshold p < 0.0125;
- the single-threaded fits and the parallel rule of A4.

The same fits give the gene count at FDR < 0.05 for the control in A5.5.

**A5.4 Correction M2, covariate (POST-HOC sensitivity).** M1 changes only the input rows, and `poscounts` size factors (a median-of-ratios estimate) are not moved much by a few dominant rows. M1 therefore cannot absorb a composition axis like the one seen in the acute check. M2 adds one covariate to the M1 table and design: `~ mtrRNA + temperature`.
- **Covariate.** `mtrRNA` is the per-library log ratio log((KEG47_r01 + KEG47_r02 + 0.5) / (all non-rRNA gene counts + 0.5)), from the multi table and centred within tissue.
- **Same procedure.** Rows, size factors and BH set are as in M1, and the same 924-labeling null applies. The covariate stays with its library under every labeling.
- **Caveat.** The covariate is measured on the same libraries as the outcome, so it could carry part of a real temperature effect. The two temperature groups are also two tanks. If the mt-rRNA share differs between them, the covariate is partly a tank variable, and M2 absorbs tank as well as composition. For that reason, M2 is a sensitivity and yields no reading of its own. Its observed count and permutation p are reported beside M1, together with the 28 °C and 31 °C means of the covariate, which are read with this in mind.
- **Expected M1 outcome.** Because `poscounts` is robust to a few dominant rows, an M1 result close to the registered one is expected. It would not be evidence that composition was unimportant.

**A5.5 Relaxed gene positive control, per tissue.** The post-hoc control differs from the registered control (5.3.1) in two respects:
1. **Per tissue.** The registered control is study-level: TE nulls are uninformative only if no tissue exceeds its 95th percentile. The post-hoc control is applied to each tissue separately.
2. **90th percentile.** The threshold moves from the 95th percentile to the 90th.

**Consequence.** Under null labels, each tissue's control passes by chance about 10 % of the time instead of 5 %. A passing control is what lets a null be read (`POSTHOC_NO_DETECTABLE_SHIFT` rather than `POSTHOC_UNINFORMATIVE`), so the relaxation makes null readings easier to call informative. A null read through this control is therefore weaker evidence than one read through the registered control, and every such reading says so.

The rule:
- **Rule.** A tissue passes if its observed gene count at FDR < 0.05 is strictly greater than the **90th percentile** of the same count over its 924 labelings, where the registered threshold was the 95th.
- **Percentile.** The percentile is the value at rank ⌈0.9 × 924⌉ = 832 of the 924 ordered counts, the observed labeling included.
- **Where it applies.**
  - (a) To the registered fits, using the stored nulls (no new fit). This part is descriptive only and gives no reading: the registered readings stand, and A5.6 reads only under M1. It does not re-label the registered muscle result.
  - (b) To the M1 fits.
- **Why 90th.** It doubles the tail in which the observed count may sit, from 5 % to 10 %, while keeping the gate a comparison with the tissue's own null. The threshold was not chosen by looking at the muscle null below the 95th percentile, which I had not opened at c9f90b0.
- **Descriptive only.** The heat-shock comparison with the paper (5.3.1) is repeated under M1, descriptively.

**A5.6 Reading rule per tissue (POST-HOC, under M1).**
- **`POSTHOC_TANK_TEMPERATURE_SHIFT_DETECTED`:** at least one classified TE family at FDR < 0.05 and permutation p < 0.0125. It is reported whether or not the tissue's relaxed control passes, with the control result beside it.
- **`POSTHOC_NO_DETECTABLE_SHIFT`:** no family passes, or p ≥ 0.0125, and the tissue's relaxed control (A5.5 b) passes.
- **`POSTHOC_UNINFORMATIVE`:** no family passes, or p ≥ 0.0125, and the tissue's relaxed control fails. A null cannot be read in a tissue whose own gene control fails.

The tank clause of the registered report applies to every post-hoc reading: the 28 °C and 31 °C groups are two tanks, so a shift is a difference between two tanks. The unique-mode sensitivity (5.3.2) is refitted under M1 on the observed labeling only.

**A5.7 Other tests.**
- **E3 superfamily test (A1.7).** Its size factors come from the same gene rows, so it is refitted under M1, with its 924 null, as a secondary. It runs only if the compute rule in A5.8 allows.
- **Global TE index (5.2.2).** Recomputed with the A5.1 rRNA rows removed from the denominator. Descriptive.
- **Not re-evaluated, stated in the report:** E1 (separate unique-fragment family sums), E2 and 5.4 (descriptive comparisons of fold changes), and the 5.2.4 tissue model.

**A5.8 Compute and order.** The study has used 220.04 of 250 CPU-h. No new alignment or counting is done; everything runs from the existing TEcount tables on S.
- **Projection.** The registered per-fit time was 3.6 to 7.9 s. One 924-labeling set over four tissues is about 3,700 fits, or about 4 to 8 CPU-h. The stored permuted fits averaged 5.9 s (family level) and 5.4 s (superfamily level). One M1 set is therefore about 6 CPU-h, and M1 and M2 together are about 12. E3 under M1 adds about 5 to 6. The total is about 17 to 22 CPU-h, or 237 to 242 of 250.
- **Order:** A5.2 shares; M1 observed fits, with the per-fit time measured; M1 nulls with the control; M2; E3; the index.
- **Stop rule.** The projection is redone from the measured M1 fit time before any null starts. A step that would take the study past 250 CPU-h does not start, and I stop and tell the orchestrator.

**A5.9 Outputs.** The work sits on branch `analysis/te-warming-mtcorr-20261006`:
- a dated section in `RUN_RECORD.md`;
- scripts and tables under this folder's `mtcorr/` subfolder;
- a hive report `te-warming-mtcorr-<date>` (md and pdf only), reviewed by Phyllis (science) and Andy (numbers) at one sha, then a PR.

The registered report and tables are not edited.
