# TE expression under acute CO2 and ocean warming

This repository brings together two completed transposable element (TE) expression analyses. Both were run under frozen pre-registrations and passed science and numbers review. Nothing here is a new analysis. Every file under `analyses/` and `reports/` is a byte-for-byte copy of its source, and `SOURCES.md` lists where each came from. `MANIFEST.sha256` holds the checksums.

| Analysis | Data | Folder | Report |
|---|---|---|---|
| Acute CO2 (5 days, about 1,000 µatm) in adult zebrafish brain and gill | PRJNA1102381, 26 libraries (Sourisse et al. 2025) | `analyses/zf-acute-co2/` | `reports/md/zf-acute-20261006.md` |
| Warming (28 vs 31 °C) in post-larval *Amphiprion ocellaris*, four tissues | PRJNA970523, 48 libraries (Moore et al. 2024) | `analyses/warming-aocellaris/` | `reports/md/te-warming-aocellaris-20261003.md` |
| Warming, post-hoc re-check after removing the rRNA genes (2 mitochondrial, 93 nuclear) from gene counts | same 48 libraries | `analyses/warming-aocellaris/mtcorr/` | `reports/md/te-warming-mtcorr-20261007.md` |

There is no acute warming experiment. The warming study has no acute arm, and the only acute exposure here is CO2.

## Readings and their limits

These readings are taken from the reports. Each limit belongs with its reading and should be quoted with it.

**Acute CO2, zebrafish.**
- Gill: `NO_DETECTABLE_TE_SHIFT`. 0 of 870 TE subfamilies are at FDR 0.05 (permutation p 1.0). With 6 against 6 libraries, and mitochondrial 16S rRNA taking 44 to 83 % of gene counts, the test had power only for broad or strong shifts.
- Brain: `UNINFORMATIVE`. The registered brain gene positive control failed, so the brain result cannot be read as evidence about TEs, either as a null or otherwise.

**Warming, *A. ocellaris*.**
- Each temperature was one tank of 6 larvae from a single breeding pair, so every result is a difference between two tanks. No analysis here can separate temperature from tank.
- Brain (8 families), digestive tract (22) and liver (24) read `TANK_TEMPERATURE_SHIFT_DETECTED`. The shifts are mostly downward in the 31 °C tank and modest, with a median |log2 fold change| of 0.22 to 0.72. Brain has fewer than 10 passing families and is flagged as underpowered.
- Muscle failed the registered gene positive control, so its 0 families is not a negative result.
- Post-hoc (Amendment 5): removing the rRNA genes (2 mitochondrial, 93 nuclear) from gene counts leaves the TE calls unchanged. Those genes hold only 0.03 to 0.30 % of gene counts, so this does not test whether library composition more broadly drives the shifts; the covariate sensitivity cannot answer that in the digestive tract, where the covariate separates the two tanks. Muscle passes a relaxed (90th percentile) control by 7 genes, so its null is weak evidence, and the aggregate TE index in muscle is still lower in the 31 °C tank. These post-hoc readings do not upgrade or change the registered ones.

Nothing here bears on transposition or mobilisation.

## Layout

- `analyses/<name>/`: pre-registration, registration or amendments, run record, scripts, compact result tables. Copied from TE_Ocean_Acidification `paper/results/`. The written reports are only under `reports/`.
- `reports/md`, `reports/pdf`, `reports/png`: the written reports with their figures. The directory layout follows the source, so the image links in the Markdown files resolve.
- `reports/py`, `reports/tsv`: the figure script and tables for the mito-rRNA report.

## Data and reruns

Raw reads are public in the European Nucleotide Archive (ENA) and the NCBI Sequence Read Archive (SRA) under the two project accessions above. Counts and alignments stay on the analysis hosts and are not in this repository. Each `RUN_RECORD.md` gives the host, tool versions, input checksums and compute used.

The scripts record the paths of the original runs (home-relative working directories and, in three scripts, an absolute R path on the server). They are kept unedited so that the files match their source. A rerun on another machine has to change those paths.
