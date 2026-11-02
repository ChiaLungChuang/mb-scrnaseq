# Proposal: which 8 tumors to re-quantify in phase 1

**Status: accepted by Chia-Lung on 2 Nov 2026: the 8 as proposed** (plan
M2 step 2). Recorded in `docs/results/phase1/samples.tsv`.

Made on 1 Nov 2026 from public metadata only:

- the ENA run table for PRJNA649773 (file names, sizes, MD5s);
- GEO GSE155446's sample annotations (SOFT file: titles, "tumor
  subgroup", protocol text);
- the first 4 MB of each human BAM (header, and the first 50,000-79,000
  reads' barcode and UMI tags);
- the **first line only** of `GSE155446_human_raw_counts.csv.gz` (cell
  IDs, i.e. sample labels and barcodes).

No count value and nothing from the cell metadata file was read
(CLAUDE.md rule 1). The full table of all 30 human runs is
[`gse155446_runs.tsv`](gse155446_runs.tsv).

## What the 30 human runs are

1. **30 libraries, 28 tumors.** Two GEO samples are second libraries:
   - GEO "934-repeat" (SRR12354571) is labelled **943** everywhere else:
     BAM `1_943_possorted_genome_bam.bam`, read group `1_943`, published
     label `943` (1,757 cells). 67% of its first reads carry barcodes of
     published label 943, and 0.2% carry those of 934. Tumor 934 has its
     own library (SRR12354570, label `934`).
   - GEO "966-recurrence" (SRR12354593) is published as **966-2**.
   The cohort count of 28 (1 WNT, 9 SHH, 7 Group 3, 11 Group 4) matches
   once these two are folded into tumors 934 and 966. Whether 943 is a
   repeat library of tumor 934 or a mislabelled separate tumor cannot be
   told from metadata; GEO and the cohort count say repeat.
2. **Every run maps to exactly one published label,** checked by barcode
   overlap, not by names: 20-69% of a BAM's first reads carry barcodes of
   its label's published cells, and at most 8% those of any other label
   (`frac_reads_in_label` and `next_label` in the TSV). The rest are reads
   from barcodes the authors did not keep as cells.
3. **Read-group names do not always carry the label.** Four BAMs are named
   `Foreman_1` to `Foreman_4` (tumors 801, 877, 996, 1070). The rule "map
   by read group or file name" in CLAUDE.md and the runbook does not work
   for them; the TSV is the mapping to use.
4. **Same reference and Cell Ranger version everywhere.** All 30 headers
   name `refdata-cellranger-GRCh38-3.0.0` and `STAR_2.5.1b`; GEO's
   processing text says Cell Ranger 3.0.1. This confirms
   `cellranger_compat: cr3` and the GRCh38-3.0.0 reference for every run.
5. **Chemistry:** v3 (12-nt UMI) for all but **1355**, which is v2 (10-nt
   UMI). Barcodes are 16 nt everywhere.
6. **cDNA reads are 151 nt.** The pipeline builds its STAR index with
   `star_sjdb_overhang: 100`. That matches Cell Ranger only if Cell
   Ranger's own index used STAR's default of 100; check
   `star/genomeParameters.txt` in the downloaded reference (runbook step
   4) before locking. If it says otherwise, that is a decision for you, as
   it changes alignments slightly.
7. **Four tumors share barcodes beyond chance.** Published labels 1238,
   1397, 1416 and 1433 share 52-148 cell barcodes pairwise, where chance
   on the 6.8M v3 whitelist predicts fewer than one. They were probably
   processed together. Phase 1 is unaffected (each label is compared only
   with its own BAM); phase 2 might want to note it.

## Proposed 8

Principles, in order:

1. One library per tumor (no 934 with 943, no 966 with 966-2).
2. Every subgroup, roughly in proportion to the cohort.
3. Both chemistries: 1355 is the only v2 run.
4. Every naming batch. The six naming patterns (`name_pattern` in the
   TSV) look like separate processing batches, and batch-specific
   differences are what criteria C1-C4 exist to catch.
5. Within those constraints, smaller BAMs: download, bamtofastq and
   STARsolo time grow with reads. (Concordance should not depend on
   depth, but this is a bias to state, not hide.)

| Run order | Label | Subgroup | Chem. | SRR | BAM GB | Published cells | Naming batch | Why |
|---|---|---|---|---|---|---|---|---|
| 1 | 1355 | Group 3 | v2 | SRR12354589 | 6.4 | 2,159 | `<n>_<label>` | Pilot: smallest BAM, and the only v2 run |
| 2 | 1195 | Group 4 | v3 | SRR12354584 | 6.6 | 1,030 | `<label>_<n>` | Smallest of the largest batch |
| 3 | 1325 | SHH | v3 | SRR12354588 | 8.2 | 1,510 | `<label>_<n>` | Small, and the most cells of the small SHH BAMs |
| 4 | 831 | SHH | v3 | SRR12354566 | 12.0 | 620 | `<n>_<label>` | Second SHH, another batch |
| 5 | 925 | Group 3 | v3 | SRR12354569 | 12.8 | 1,170 | `<n>_<label>_010819` | Smallest of its batch |
| 6 | 1066 | Group 4 | v3 | SRR12354576 | 15.5 | 1,333 | `Foreman_<label>_03142019_<n>` | Group 4 from this batch; 898 is SHH and 1155 is heavily filtered (see alternatives) |
| 7 | 1128 | WNT | v3 | SRR12354579 | 26.2 | 832 | `Foreman_<label>` | The only WNT tumor |
| 8 | 996 | Group 4 | v3 | SRR12354574 | 31.1 | 704 | `Foreman_<n>` | Smallest of the batch whose read group has no label |

WNT 1, SHH 2, Group 3 2, Group 4 3; both chemistries; all six naming
batches; about 119 GB of BAMs in total.

**Alternatives,** if you prefer:

- **Cheaper:** replace 996 (31.1 GB) with 934 (Group 4, 10.6 GB, batch
  `<n>_<label>`). Saves about 20 GB but leaves the `Foreman_<n>` batch
  (the largest BAMs) untested.
- **A label-trap test:** add or swap in 943 ("934-repeat", 23.1 GB), to
  show that the pipeline's samplesheet uses the verified label rather than
  the GEO title.
- **A cell-calling stress test:** 1155 (Group 4, 15.3 GB). Only 458
  published cells, and 20% of its reads carry their barcodes, so the
  authors filtered it hard; the `published_cells_called_by_starsolo`
  metric would show how far STARsolo's own cell calling differs.

## Decisions for Chia-Lung

1. Accept the 8, or swap any of them (for example the alternatives above).
2. Nothing else here is a decision yet. Item 6 (STAR overhang) becomes one
   only if the reference's `genomeParameters.txt` does not say 100.

## Follow-ups once decided

- Commit the chosen 8 as `docs/results/phase1/samples.tsv` (plan M2 step
  5), taking labels from `gse155446_runs.tsv`.
- Update the label-mapping rule in CLAUDE.md and runbook step 7: map by
  this table. Read groups `Foreman_1` to `Foreman_4` carry no label, and
  the SRA/GEO title is wrong for 943.
- Phase 2 (PR #4, D7): 934 and 943, and 966 and 966-2, must be grouped as
  one patient each in the classifier's folds, or one library of each pair
  dropped.

## How the table was made

Reproducible with `curl` and `samtools` from the sources above:
`scientific_name`, `submitted_ftp`, `submitted_bytes` and `submitted_md5`
from the ENA filereport API; `tumor subgroup` and titles from
`GSE155446_family.soft.gz`; `samtools view -H` and the CR/UR/CB tags of
the first reads from a byte-range download (bytes 0-4,000,000) of each
BAM; the cell IDs from the counts file's first line, read from a
3 MB byte-range download and discarded after use.
