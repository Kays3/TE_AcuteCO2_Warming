# Run record: TE warming, A. ocellaris (PRJNA970523)

Host pSS (32 cores), kaisar account, no sudo, every process `taskset -c 0-20` and `/usr/bin/time -v`. Working directory `~/jim_te_warm`. Ceiling 200 CPU-h (Amendment 2).

## Gates and releases (UTC)

- 2026-10-03 09:37 — base pre-registration passed at da4d4a8: Andy (numbers, 09:37:24), Phyllis (science, 09:37:09).
- 2026-10-03 09:38:43 — Michael (82dc75) released (a) download of the 48 poly(A) runs with checksums and (b) `jim_te_warm` environment setup. Not released until both reviewers pass Amendment 1: STAR indexing or alignment, counting, and any look at expression values (no QC beyond read counts and checksums).

- 2026-10-03 09:45:35 — Phyllis PASS (science) and Andy PASS (numbers) on Amendment 1, its addendum and Amendment 2 (29ccc4c / 1a0d79a). Under Michael's 82dc75 ruling, indexing, alignment and counting are released from this time.
- 2026-10-03 ~09:50 — Amendment 3 (TE GTF keeps all classes, as OA did) sent for review; the TE GTF build waits for it.

- 2026-10-03 10:32:12 — Michael (e0065b): full release of the registered chain on pSS under the 200 CPU-h ceiling.
- 2026-10-03 10:33 — DESeq2-row clarification (follow OA: drop rRNA only) sent to Phyllis and Andy; must be settled before the first DESeq2 fit.

## Environment (09:40–09:41 UTC)

- Frozen OA env, read-only by absolute path: `~/KD/workspace/envs/te-ocean-acidification` (STAR 2.7.11b, fastp 1.3.6, samtools 1.23.1, TEcount 2.2.4, DESeq2 1.46.0 on R 4.4.3, TEsorter). Export: `~/jim_te_warm/env_te-ocean-acidification_frozen.yml` (532 lines).
- `jim_te_warm` (micromamba 2.9.0, user-space binary sha256 366cd9cd…, root `~/jim_te_warm/mamba`): sra-tools 3.4.1, seqtk 1.5-r133, Subread 2.0.8. Export: `~/jim_te_warm/env_jim_te_warm.yml`. Create took 14.5 s user + 3.1 s system.

## TE reference files (sha256, verified on pSS)

- species library `~/jim_te_warm/ref/aocel-families.fa`: 16f2773127d2b942017a93a6e3cb19e8fdcd27f2e19c9ac05c1984cfb654e038 (1,849 sequences)
- `genome.fa.out`: 9f51c6e7a99370dd1b982a667db88223f22552db6495dacf5ddba5edfa367f67
- `genome.fa.align`: 9d0105e50ff09b18f0f25501fdaa1659df493e67bcd7d199b0444f3236434440
- `genome.fa.out.gff`: 33fe5dfd1b545e9868e3b15a77ce0226def7385ad9a1994ea35590ceb8ed3cef
- source: `~/KD/workspace/TE_OA_XEON_TRANSFER/mobilome-comparative-2026-09-01/results/GCF_022539595.1/`

## Gene reference (NCBI FTP, md5 verified against md5checksums.txt)

- `GCF_022539595.1_ASM2253959v1_genomic.fna`: sha256 e139d17c9078af5fce5f72272534af9193783af9bb114f3831624aef11df6071 (874,336,546 bytes). Identical to the Mobilome panel manifest hash, so the RepeatMasker coordinates of A1.1 refer to this exact FASTA.
- `GCF_022539595.1_ASM2253959v1_genomic.gtf`: sha256 c63a989abf5eff2a3ee27d99d88e86958bab96a63db10b0153c25d33529f450f, annotation RS_2023_03 (dated 29 March 2023).

## TE GTF (Amendment 3, both reviewers passed 09:49 UTC)

Built 2026-10-03 ~09:52 UTC on pSS with the Mobilome converter `pipelines/repeat-annotation/repeatmasker-out-to-gtf.py` (sha256 bd25f9b2…, TE_Mobilome_Paper main 297d876), unchanged, from the A1.1 `.out` (sha256 verified before the run). Output `~/jim_te_warm/ref/aocel.de-novo-te.gtf`, sha256 d8b3625ebe12af7cb10a115c3bbbd573bdd1cb6dabdd1daa0ae050b8272badb5, 1,432,460 loci (equal to the `.out` locus total). 3.7 s CPU.

Loci by class, this GTF vs the OA GTF (`apoly.de-novo-te.gtf`, 71bb6589…, counted on pSS):

| class | A. ocellaris | OA (A. polyacanthus) |
|---|---|---|
| DNA | 398,887 | 456,652 |
| LINE | 126,433 | 148,630 |
| LTR | 85,285 | 105,985 |
| SINE | 9,648 | 15,507 (+337 `SINE_`) |
| RC / PLE / Retroposon | 4,772 / 4,969 / 653 | 10,040 / 5,388 / 469 |
| Unknown | 484,082 | 562,148 |
| Simple_repeat / Low_complexity / Satellite | 280,114 / 34,724 / 2,208 | 266,254 / 35,320 / 2,830 |
| rRNA / tRNA / snRNA | 395 / 265 / 25 | 1,438 / 178 / 0 |

## STAR index

Launched 09:49:43 UTC, `scripts/run_star_index.sh`: STAR 2.7.11b, `--sjdbOverhang 149 --genomeSAindexNbases 13`, 16 threads. The first launch used the default 14; STAR warned that 14 is too large for an 863 Mb genome and may segfault at mapping, and recommended 13 (min(14, log2(genome)/2 − 1) = 13.7). The pre-registration does not fix this parameter; 13 is STAR's own rule, and the first launch was stopped after 30 s.

## Strandedness and per-library pipeline

- Strandedness test launched 10:33:26 UTC, `scripts/strandedness.sh`, on SRR24489122 (the first brain library in `design.tsv`, 28 °C): full `fasterq-dump`, `seqtk sample -s20261003` 2,000,000 pairs, STAR `--quantMode GeneCounts`.
- **Strandedness result (10:36 UTC): `--stranded reverse` for all 48 libraries.** Over genes in `ReadsPerGene.out.tab`: column 3 = 35,367, column 4 = 1,724,573, ratio column 3 / (column 3 + column 4) = 0.0201, below the registered 0.2 line (unstranded column 1,710,877). Subset uniquely mapped 95.18 %. CPU: seqtk 2 × ~6 s, STAR 172.6 s user.
- Pilot launched 10:38 UTC (`scripts/pilot.sh`): SRR24489122 then SRR24489123 (both brain, 28 °C), sequential, for timing.
- Pilot library 1 (SRR24489122), CPU seconds (user + system) and wall: fasterq-dump 138 s (0:21), fastp 444 s (1:18), STAR 2,731 s (3:21, 11.2 GB RSS), samtools index 20 s, TEcount multi 1,834 s (31:53, single-threaded, 2.1 GB RSS); TEcount uniq pending. Uniquely mapped 95.87 %. BAM 2.96 GB. Projection about 1.9 CPU-h per library, about 93 CPU-h for 48.
- 11:30 UTC: bulk run queued to start automatically when the pilot finishes (`scripts/after_pilot.sh` → `run_all.sh`, 6 libraries at a time, because single-threaded TEcount sets the wall time). pSS only; thinkstation1 offered by the human at 11:25 UTC and not used (see below).
- Per-library pipeline: `scripts/process_one.sh <run> <stranded>` (fasterq-dump, fastp 1.3.6, STAR 2.7.11b, samtools index, TEcount 2.2.4 multi and uniq; checks both GTF sha256 on start; FASTQ deleted after alignment).

## Download

Launched 2026-10-03 09:41:44 UTC with 2 runs in parallel, relaunched 09:47:34 UTC with 4 (`scripts/launch_download.sh`, prefetch resumes partial files) after adding Andy's checks: the spot count must equal runinfo `spots` or the run fails, and the file type (`sra` or `sralite`) is recorded. `.sra` only, `vdb-validate`, sha256, md5 and spot count per run into `sra_checksums.tsv`. FASTQ extraction waits for the processing release.

- 2026-10-03 15:30 UTC: 46 of 48 runs had a `.done` record. SRR24489125 and SRR24489126 failed at 09:47:54 UTC (prefetch exit 3 after about 3 s, "lock exists while copying file", then "transfer interrupted"). That was the window when two download dispatchers overlapped after the 09:47 relaunch. Their partial `.tmp` files (157 MB and 141 MB) were deleted, the failed logs were moved to `logs/download/old/`, and both runs were re-fetched at 15:30 UTC with the unchanged `download_sra.sh` (P=2, `taskset -c 0-20`, list `logs/refetch_2.tsv`). The same spot-count and validation checks apply. `run_all.sh` starts a library only once its `.done` record exists, so the bulk run picks both up without a restart.

## Cumulative CPU-h

Status at 15:33 UTC, 3 October. User plus system time is summed from the `/usr/bin/time -v` files on pSS.

| Stage | CPU-h | Cumulative | Ceiling |
|---|---|---|---|
| Download, validation, checksums (95 time files) | 1.23 | 1.23 | 200 |
| TE GTF build and STAR index | 0.91 | 2.14 | 200 |
| Strandedness test | 0.05 | 2.19 | 200 |
| Pilot, 2 libraries | 4.08 | 6.27 | 200 |
| Bulk, first 6 libraries | 15.63 | 21.90 | 200 |

Bulk libraries run at 2.35–2.78 CPU-h each, above the pilot's roughly 2.0, because six run at once. Wall time is 1:36–1:55 per library. At 2.6 CPU-h for each of the 40 remaining libraries, counting closes at about 126 CPU-h in total. The estimated finish is about 03:30 UTC on 4 October.

### Counting stage closed (recorded 2026-10-04 07:50 UTC)

`run_all.sh` wrote `ALL_DONE 2026-10-04T04:02:29Z`. All 48 libraries have both count tables (`.uniq.cntTable`, `.multi.cntTable`), and `logs/run_all.log` has no `FAIL` line. The last two libraries counted were SRR24489125 and SRR24489126, after their re-fetch.

The registered measure is the per-process sum. Section 6 and Amendment 2 put every process under `/usr/bin/time -v` and count user plus system time against 200 CPU-h. The per-library wrapper (`logs/run_<run>.time`, around `process_one.sh`) gives a lower total, which is recorded for comparison only.

| Stage | CPU-h (per process) |
|---|---|
| Download, validation, checksums (99 time files, incl. re-fetch) | 1.25 |
| TE GTF build and STAR index | 0.46 |
| Strandedness test | 0.05 |
| 48 libraries: fasterq-dump 2.34, fastp 12.64, STAR 121.66, samtools index 0.30, TEcount multi 32.82, TEcount uniq 30.74 | 200.50 |
| **Total, registered measure** | **202.26** |
| Wrapper sum, 46 bulk + 2 pilot (comparison only) | 184.87 |

The 17.4 CPU-h gap between the two measures has no established cause. Both are reported. The registered measure is the one compared with the ceiling.

Correction to the 15:33 UTC status table above. Its "TE GTF build and STAR index" row (0.91) counted `star_index.time` twice, through two overlapping file patterns. The correct value is 0.46.

**Deviation (dated 2026-10-04): the compute ceiling was reached during counting.**
- Counting alone used 202.26 CPU-h of the 200 CPU-h ceiling. The 15:33 UTC projection had been about 126 CPU-h.
- **Cause:** oversubscription. Six libraries ran at once inside the 21-core pin (`taskset -c 0-20`), each with a multi-threaded STAR.
  - STAR averaged 9,124 CPU-s per library, against 2,731 CPU-s in pilot library 1. Pilot wall time was 3:21, against 4:33 for SRR24489126.
  - Threads competing for the same cores inflated user time without speeding the run up.
  - The other steps stayed near pilot cost: TEcount multi 2,461 s, uniq 2,306 s, fastp 948 s per library.
- **Missed check:** the projection was not recomputed at the 17:31 UTC status check, when 14 libraries had finished. A recompute then would have shown the drift, and the run could have been stopped at about 120 CPU-h spent.
- **Status:** the stop rule (section 6, Amendment 2) binds. No analysis step has run. Michael was told at 07:38 UTC with three options, and the decision is with the human (ASK ME board, OHV-35). pSS has run no job of mine since 04:02:29 UTC.
- **Storage:** nothing is deleted. The BAMs take 144 GB in `out/`. The `.sra` files and FASTQ were deleted per library after alignment by `process_one.sh`, as registered. `sra_checksums.tsv` keeps the sha256, md5 and spot count of all 48.

## Analysis

- 2026-10-04 07:50:45 UTC: the human raised the ceiling ("increase the ceiling to reach the capacity"). Michael read it as 250 CPU-h, and that reading is Amendment 4 (e76f808). Review requests went to Phyllis and Andy at 07:51 UTC. No permuted fit runs before both pass.
- `scripts/deseq_fit.R` implements 5.1, 5.2.3 and 5.3.1 as one procedure for the observed and permuted labelings:
  - the joint gene + TE table, with `rRNA` rows dropped;
  - the filter (≥ 10 counts in ≥ 6 libraries) and `poscounts` size factors, computed once per tissue and held fixed;
  - `DESeq()` with default settings (`fitType` parametric), serial;
  - `results()` with `independentFiltering = FALSE`;
  - BH applied separately within the classified-TE set, the Unknown set and the gene set.
  - Labelings follow the `combn(12, 6)` order over the tissue's runs sorted by run id. The chosen set is the 31 °C group, and index 0 is the observed labeling.
- **Implementation choices fixed before any permuted fit:**
  - **No LFC shrinkage.** The registered text names none, and apeglm is absent on pSS. CIs are MLE ± 1.96 SE.
  - **5.2.3 statistic.** The text says "median |LFC| of TE families against genes of matched base mean, 10 deciles". It is implemented as: median |LFC| of classified TE families ÷ the median, over those families, of the gene median |LFC| in the family's gene-baseMean decile.
- **Launcher note.** The frozen environment's `Rscript` is a shell wrapper that drops `--args`. The first launch, at 07:53:12 UTC, therefore failed the argument check after 0.2 s in all 4 processes, with no computation. The script is called as `R --no-echo --no-restore --file=deseq_fit.R --args …` instead.
- **Observed fits**, multi mode, 07:53:47–07:54:02 UTC: 4 processes, one per tissue, `taskset -c 0-20`, `/usr/bin/time -v`.
  - Process CPU: brain 11.3 s, liver 10.6 s, muscle 14.6 s, digestive tract 12.1 s. That includes table building. Maximum RSS is about 0.8 GB.
  - DESeq2 fit time alone: 3.6–7.9 s.
- **Projection after the timing fits (A4.3):**
  - Spent: 202.27 CPU-h.
  - Permutations, 4 × 924 fits at ≤ 8 s: about 8 CPU-h.
  - Unique-mode observed fits, the 48-library tissue model and E2: under 0.5.
  - E3 with its labeling null, a similar table: about 8.
  - E1: 8–15.
  - **Total about 227–234 CPU-h against 250.**

**Corrections (appended by 08:06 UTC, 4 October; Amendment 4 addendum).**
- The heading "Counting stage closed (recorded 2026-10-04 07:50 UTC)" should read "by 07:40 UTC", the time of commit 23085b4.
- The deviation sentence "the run could have been stopped at about 120 CPU-h spent" is wrong. At 17:31 UTC on 3 October, 14 libraries were finished and spend was roughly 60–70 CPU-h. That is 21.9 CPU-h for the first 8, plus 6 × 4.52 CPU-h, the mean of the last 40 bulk libraries, plus six libraries in flight.
- In the 15:33 UTC table, the cumulative column is 0.45 CPU-h high from the "TE GTF build and STAR index" row onward, after the 0.91 → 0.46 correction. No later figure depends on it.
- **Observed fits, frozen by sha256 (Phyllis R5).** Their summary rows were printed in the agent's session when timing was collected (disclosed in the Amendment 4 addendum).

| File | sha256 |
|---|---|
| analysis/observed_brain_multi.tsv | 30c31c6e85e3560d7e1951317a935b94ddc3c4a991f08ecd521f9d80f6d9bd4c |
| analysis/observed_digestive_tract_multi.tsv | b4157483dcc7df496b5f21f889545e176e8e296ace1601406b7000a21f776b36 |
| analysis/observed_liver_multi.tsv | b70907d30522866f7a33fd189c4aa3f6e4b7cc90c1b5e3be1d53b7936ce14136 |
| analysis/observed_muscle_multi.tsv | c9e9064818e9c3e6232536ee313ea931f5583700759da953f1aff66834fc790f |
| analysis/perm_{brain,digestive_tract,liver,muscle}_multi_0000_0000.tsv | a662f3c0…, 0e931423…, 491c7f11…, fde6b1fd… |
| analysis/prep_{brain,digestive_tract,liver,muscle}_multi.rds | dffe0551…, edc45561…, f12a26f1…, 344570b3… |

### Analysis stage (recorded by 09:30 UTC, 4 October)

- **5.1 and E3 nulls.** `scripts/run_perms.sh`, launched 08:01:11 UTC, at most 14 single-threaded jobs. 44 family chunks, 4 E3 observed fits and 44 E3 chunks, each of 84 labelings. All chunk outputs were complete before the 09:05 UTC summary. `summarize_nulls.py` asserts 924 distinct labelings per tissue and level, and checks that the observed labeling's refit (index 0) equals its count inside the null.
- **Launcher flaw (disclosed).** In `run_perms.sh` and `run_e1_stage4.sh`, the final `wait` sat outside the pipeline subshell that owns the jobs. The "DONE" lines in the status logs (PERM_DONE 08:49:16, E1S4_DONE 09:12:09) were therefore written when the last job was launched, not when it finished. No result depends on those lines: every summary checks completeness itself (924 labelings per tissue and measure). The last E1 chunk wrote its output at 09:14 UTC, before the E1 summary. Both launchers were fixed at 09:13 UTC.
- **E1 environment addition.** `jim_te_e1` (micromamba, pSS) holds the mobilisation study's versions: python 3.12.14, pandas 2.3.3, numpy 2.5.3, pysam 0.24.1, biopython 1.88, subread 2.0.8 (featureCounts) and edlib 1.3.9.post1 (pip). pSS had no pandas or edlib. This is an environment addition, not a method change; featureCounts is Subread 2.0.8, as A1.2 registers.
- **E1 inputs.** The committed mobilisation scripts are used unchanged, checked by sha256 on start: `stage1_candidates.py` 3acc6acd, `stage2_regions.py` ba99035a. The A3 version is byte-identical to the dgx run. featureCounts options are those of `s2one.sh` / `s2extra_one.sh` (`-T 6 -p --countReadPairs -s 0 -O --minOverlap 1`; 1/NH pass `-M --fraction`).
  - `e1_stage3_calls.py` differs from `stage3_calls.py` in one line, the tissue tuple. Sample ids are `<tissue>_<run>`, with digestive tract written `digestivetract` because Stage 3 splits on "_".
  - `e1_stage4_warming.R` replaces the CO2 designs with `~ temperature` per tissue. It uses the same labeling enumeration as 5.1, size factors from gene counts (poscounts) and BH within the classified and Unknown sets.
  - **E1 outcome rule, fixed before the E1 null was read:** a tissue is `TANK_TEMPERATURE_SHIFT_INTACT` on a measure (unique or 1/NH) only if both conditions hold on that measure. The two measures are reported separately and never combined.
- **Stage 1 (A. ocellaris).** 1,010,533 copies; tiered loci T1 328, T1u 6, SINE-T 145, T2 6,181, C1 76,360.
- **E2.** The OA library is read as sha256 152a2a7e. A second file named `apoly-families.fa` on pSS has sha ca2da39f and was not used.
  - rmblastn 2.17.1 (tetools 2.00) needs explicit scoring, so it runs with `-reward 1 -penalty -2 -gapopen 5 -gapextend 2 -evalue 1e-10`. The first E2 attempt stopped at the option check.
  - The containers run under the Docker daemon, so their CPU is not in the `/usr/bin/time` files. The whole E2 step took 6 s of wall time and is counted as under 0.01 CPU-h.
  - **"Merged OA results"** is implemented as the mean, over the two parental phenotypes, of the full-set developmental (and, separately, intergenerational) log2 fold change in `response-patterns.tsv.gz` (sha256 2699acba of the uncompressed content).
- **5.2.2.** `index_build.py` is the CO2-validation builder with only paths changed.
- **5.3.3.** Median CpG-adjusted Kimura per family from the `.align`: 1,844 families, 1,323,398 records.
- **A1.9.** Exon-overlap share for the 35 Unknown families at FDR < 0.05.

| Analysis step | CPU-h (per process) |
|---|---|
| Observed timing fits | 0.014 |
| 5.1 and E3 nulls (92 processes) | 11.340 |
| E1 environment, Stage 1, region build | 0.016 |
| E1 featureCounts (144 passes) | 4.616 |
| Batch 2 (E1 Stage 3, E2, index, Kimura, unique-mode observed fits, 5.2.4 tissue model) | 0.033 |
| E1 Stage 4 with null (48 processes) | 1.757 |
| **Analysis total** | **17.78** |
| **Study total (202.26 + 17.78)** | **220.04 of 250** |

Every process ran inside `taskset -c 0-20` with at most 20 compute threads at once. MemAvailable stayed above 90 % at every check. Nothing was deleted: the BAMs (144 GB) stay until Michael's word.

**Correction (appended by 09:40 UTC, 4 October; Phyllis R3, PR #181).** The analysis-stage note says the E1 outcome rule was "fixed before the E1 null was read". Commit order does not support that. The rule is in 951f6b0, committed 09:24:35 UTC, after the E1 Stage 4 null finished (by 09:14 UTC). Read it as: the per-measure rule was written down after the E1 null was computed. The outcome is identical under any combination of the two measures (per measure, both required, or either): liver passes on both (3 unique, 7 1/NH, both p = 0.0065) and no other tissue passes on either.

**Exclusion check (section 4; Phyllis R1).** No library was excluded. Uniquely mapped rates from the 48 STAR `Log.final.out` files are 57.22 to 95.95 % (median 93.06 %), all above the 50 % line; TEcount produced both tables for all 48. Per-library input pairs, unique and multimapped rates are in `tables/star_qc_48.tsv`. Digestive tract maps lowest, and lower in the 31 °C tank (median 77.4 %, minimum 57.2 % in SRR24489106) than in the 28 °C tank (median 85.9 %). This is stated in the report as a possible compositional confound of the digestive-tract comparison. It was not examined further.

**Source tables and a numerical note (appended by 09:50 UTC, 4 October; Andy, PR #181).**
- Every number quoted in the report now has a committed source in `tables/`:
  - `heatshock_dnaj_genes.tsv`: heat-shock and DnaJ genes per tissue, observed fit;
  - `tested_families_by_class.tsv`: classified families tested per tissue. DNA is 48–51 % of them;
  - `te_results_observed_multi.tsv.gz`: all classified and Unknown TE rows of the four observed fits, with CIs;
  - `tissue_model_te_pca.tsv`: PC1 30.6 %, PC2 19.1 %;
  - `e2_rbh_808080.tsv`: the 167 one-to-one 80-80-80 pairs.
- `index_5_2_2.tsv` now prints p in 4 significant digits; the liver Welch p is 1.373e-05.
- **Mirror refits.** A labeling and its mirror should give identical |log2 fold change|. The upper-tail counts of the 5.2.3 ratio are odd in some rows, for example 441 of 924 for brain at family level, so mirror pairs do not always tie exactly. The likely cause is floating-point difference between the two refits, since the reference level flips and dispersion fitting is iterative; it was not traced further. The FDR counts that decide the outcomes tie between mirrors in every tissue checked (the observed labeling's count equals its index-0 refit). The ratio is descriptive, and every row is far from any threshold.
- **Precision on the mirror note above.** "Tie between mirrors in every tissue checked" originally rested only on the observed-against-index-0 check. It has now been checked directly, by pairing each labeling with its complement in `null_distributions.tsv.gz`:
  - The classified-TE FDR count ties in all 462 mirror pairs, for every tissue at both family and superfamily level.
  - The gene FDR count ties in all pairs except one liver pair at family level.
  - So the outcome counts are mirror-exact, and floating-point differences show up only in the continuous ratio and in one gene count.
- Correction (Andy, PR #181): the DNA share of tested families is 48.1-51.6 %, i.e. 48-52 %, not 48-51 %.
- Erratum (appended 4 October, after the merge at 1ce33ac; text only): the report's Summary and Interpretation point 3 said the liver E1 shift sits in families with no intact copies. That omits the one exception stated in Results: ltr-1_family-2 (1/NH measure) has 2 intact copies, INTACT_NOT_DETECTED. Both now read 'every passing family except one has no intact copy'. Interpretation point 2 changes from 'does not exceed' to 'is not detectably larger than', matching the descriptive 5.2.3 wording in Results. No result, table or figure changed.

## Amendment 5 run (post-hoc mt-rRNA correction; appended 7 October, 03:07 JST (18:07Z))

Design: `PREREGISTRATION.md` Amendment 5 at 5b09082 (passed by Andy and Phyllis). No new alignment; every fit reads the existing TEcount tables on workstation S. Scripts and tables are in `mtcorr/`.

- **Observed fits.** 8 M1 observed fits (4 tissues, multi and uniq) finished 23:04 JST, 6 October (14:04Z).
- **Nulls.** Each tissue was fitted on all 924 labelings for three jobs:
  - M1 family: done 23:30 JST (14:30Z);
  - M2 family: done 23:50 JST (14:50Z);
  - M1 superfamily (E3): done 00:14 JST, 7 October (15:14Z).
- **Completeness.** Every tissue-by-job set has labelings 0–924. No guard stop occurred, and no error appears in the 148 job logs.
- **Summary.** `scripts/summarize_mtcorr.py` was run on workstation C, because the S environment has no scipy. Before the run, the glob was tightened so the registered family nulls no longer pick up the `_sf_` superfamily files.
- **Index check (A5.7).** The strict-TE numerator and the all-count denominator match `tables/te_index_per_sample.tsv` exactly (maximum difference 0) before rRNA removal.

| Item | CPU-h (user + sys) |
|---|---|
| M1 family, observed and null | 5.911 |
| M2 family, observed and null | 4.534 |
| M1 superfamily (E3), observed and null | 5.526 |
| M1 uniq observed fits | 0.012 |
| **Amendment 5 total** | **15.98** |
| **Study total (220.04 + 15.98)** | **236.02 of 250** |

The projection was 17–22 CPU-h. Every process ran inside `taskset -c 0-20` with the memory guard on. Nothing was deleted.

Outputs:
- `mtcorr/tables/a5_results.tsv` (readings);
- `a52_shares*.tsv`;
- `a57_index*.tsv`;
- `m1_families_<tissue>.tsv`;
- `null_mtcorr.tsv.gz`.
- Heat-shock comparison under M1 (A5.5, descriptive; appended 7 October, 03:10 JST (18:10Z)). `scripts/heatshock_m1.py` reads the 251 registered heat-shock and DnaJ gene rows in the M1 observed fits. All rows are present. The count at FDR 0.05 is unchanged in every tissue: brain 8, digestive tract 11, liver 11, muscle 4. The largest change in log2 fold change is 0.00016. Outputs: `mtcorr/tables/heatshock_m1_{genes,summary}.tsv`.
