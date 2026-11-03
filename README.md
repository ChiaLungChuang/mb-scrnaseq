# mb-scrnaseq

Re-quantification of childhood medulloblastoma single-cell RNA-seq
(GEO [GSE155446](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155446):
28 tumors, 10x Genomics 3') from the deposited Cell Ranger BAMs, with a
pre-registered test of whether the result reproduces the published count
matrix. Built with Nextflow and containers, with an execution profile for
the St. Jude HPC cluster (LSF + Singularity).

**Status:** phase 1 (reads to counts, QC, concordance) is built and passes
its tests in CI (GitHub Actions), with the containers, on two synthetic
datasets with known ground truth: the R unit and accuracy tests, the script
tests, lint, and the whole pipeline under Docker
(`nf-test test --profile +docker`, with the R image built from
`containers/r-analysis` in the same run). The R image is published on
GHCR (`mb-scrnaseq-r:0.1.0`, pullable without a login) and the pipeline
pins it by digest. Nothing has been run on the St. Jude cluster (the
LSF/Singularity profile is untested) or on the real data. Phase 2 is
planned, not started: the single-cell analysis
of all 28 tumors in R (Seurat), plus a Python module that compares scVI
with Harmony integration and tests whether tumor subgroup can be predicted
from malignant-cell expression, validated by patient.

## Why re-quantify instead of using the published matrix

A published matrix is one run of one pipeline. Re-deriving it from the reads
checks that the data, the sample labels and the processing are what the
paper says they are, and gives a pipeline that can process new data the
same way. Disagreement is informative: a wrong chemistry, whitelist,
reference version or sample mapping shows up as a failed criterion rather
than as plausible but wrong biology later.

## Pipeline

```
Cell Ranger BAM ─► bamtofastq ─► read check ─► STARsolo ─┬─► per-sample QC
                                  (chemistry             │   (MAD flags + thresholds)
                                   + read layout)        └─► concordance vs published matrix
                                                              (pre-registered criteria C1-C4)
```

| Step | Tool | Notes |
|---|---|---|
| `BAMTOFASTQ` | 10x bamtofastq 1.4.1 | Recovers the original reads from each BAM |
| `CHECK_READS` | `bin/check_reads.sh` | Confirms v2/v3 chemistry from read lengths; handles both bamtofastq layouts (cDNA as R2 or R3) |
| `STAR_GENOMEGENERATE` | STAR 2.7.11b | Built once; reusable via `--star_index` |
| `STARSOLO` | STAR 2.7.11b | Options from STAR's docs for matching Cell Ranger 3.x or ≥4 (`--cellranger_compat`, required); EmptyDrops cell calling |
| `SAMPLE_QC` | R, scuttle | Flags cells by median ± 3 MAD per sample; records thresholds; removes nothing |
| `PREP_PUBLISHED` | R | Streams the dense GEO CSV into a sparse matrix in chunks; detects orientation |
| `CONCORDANCE` | R | Per sample: barcode recovery, UMI totals, pseudobulk and per-cell correlation |

All 30 GSE155446 human BAM headers (read from their first 4 MB) name Cell Ranger's
`refdata-cellranger-GRCh38-3.0.0` reference and STAR 2.5.1b, i.e. Cell
Ranger 3.x, so this dataset is run with `--cellranger_compat cr3` and the
GRCh38-3.0.0 reference. `scripts/fetch_inputs.sh bams` reports both for
every BAM it downloads.

## Pre-registration

The thresholds that decide whether the re-quantification reproduces the
published data are in
[`preregistration/acceptance_criteria.yaml`](preregistration/acceptance_criteria.yaml).
On real data the pipeline refuses to run the concordance step unless:

- the file has a top-level `status: locked`, is tracked by git and unmodified;
- the commit that last changed it is on a remote branch or tag (pushed, not
  only committed locally);
- the pipeline code has no uncommitted changes.

Every run records the criteria file's SHA-256, its last commit, a SHA-256 of
the code (`lib/R`, `bin/`) and the git commit in
`pipeline_info/run_provenance.json`, and the verdict file states whether the
run was pre-registered. Test runs bypass the gate with
`--allow_unlocked_criteria true` and are marked `preregistered: false`.

This makes the order of events checkable, not tamper-proof: git commit dates
are set by whoever commits. A tagged GitHub release of the locking commit,
made before the run, gives an independent timestamp.

## Running

Synthetic test (about a minute):

```bash
nextflow run . -profile test,docker        # or test,singularity
```

CI runs this synthetic test with Docker on every push (through nf-test; see
Tests). The `singularity` and `stjude` profiles have not been run yet.

Real data on the St. Jude cluster: see [docs/RUNBOOK.md](docs/RUNBOOK.md).

Outputs land in `--outdir`: `starsolo/<sample>/` (raw and filtered count
matrices, `Summary.csv`, STAR log), `qc/`, `concordance/`, `reads_check/`,
`reference/` (the STAR index, when built) and `pipeline_info/` (provenance,
software versions, Nextflow report, timeline and trace).

## Tests

Three layers, each checked by planting bugs and confirming the tests fail
(see `docs/DEVELOPMENT_PLAN.md`, "Testing rules").

- **R library** (`Rscript tests/run_tests.R`): unit tests for chosen cases,
  plus accuracy tests (`tests/testthat/test-accuracy.R`) that compare every
  function with an independent implementation on random inputs: the
  concordance metrics with a dense base-R version of their definitions, QC
  metrics and thresholds with `scuttle` and with median ± n MAD by hand, the
  CSV and Matrix Market readers with whole-file reads.
- **Scripts** (`bash tests/shell/run.sh`): `check_reads.sh` on hand-built
  FASTQ folders (both layouts, v2/v3, lanes and chunks, every refusal);
  `fetch_inputs.sh` with stand-in `curl` and `singularity` (MD5, size,
  truncation, missing 10x headers, dry run, SRA fallback, label parsing,
  whitelists, reference layout); `check_versions.sh` (version, image tag
  and digest agreement); `hpcf_preflight.sh` (the cluster checks,
  with stand-in LSF, modules, registries and download hosts, including
  that it changes nothing unless asked to submit its one test job); the LSF directives the St. Jude profile writes
  (Nextflow with a stand-in `bsub`: per-core memory, queue routing); a
  check that the phase 1 sample record
  (`docs/results/phase1/samples.tsv`) agrees with the barcode-verified run
  table; and a check that the committed test data is
  exactly what its generators produce.
- **Pipeline** (`nf-test test --profile +docker`), on two synthetic
  datasets with known answers:
  - `tests/main.nf.test` (`make_test_data.py`): every cell's UMI total
    matches the truth, QC flags exactly the planted low-quality cells, both
    read layouts work, swapped labels fail barcode recovery, and wrong
    chemistry, unlocked criteria, a missing `--cellranger_compat` and
    duplicate samples are refused.
  - `tests/accuracy.nf.test` (`make_accuracy_data.py`): every cell x gene
    count equals Cell Ranger's rules exactly, for v3 and v2 chemistry, on
    reads built to trip them: PCR duplicates, barcode and UMI sequencing
    errors, UMIs shared between genes, spliced, minus-strand, antisense,
    intronic, intergenic and multi-mapping reads, and two genes sharing a
    symbol.
- **Results report** (`bin/phase1_report.R`, for `docs/results/phase1/README.md`):
  unit tests with worked-out compute in `tests/testthat/test-report.R`, and
  `tests/check_phase1_report.R`, which CI runs on the synthetic runs'
  outputs (passing and failing verdicts).
- `.github/workflows/ci.yml` runs all of the above on every push to
  GitHub, with the pipeline tests on the R image built from that commit,
  plus `nextflow lint`, ShellCheck, `scripts/check_versions.sh` (release
  tag, `manifest.version` and the R image tag agree; the digest is well
  formed), and a job that pulls the pinned R image without logging in,
  checks that its tag still points to the pinned digest, and runs the R
  tests inside it.

## Known limitations

- STAR 2.7.11b is not the STAR 2.5.1b inside Cell Ranger 3.x, so alignments
  can differ slightly. The acceptance criteria are set to tolerate that, not
  to demand identical counts.
- Genes are matched to the published matrix by symbol, because that is all
  it has. Symbols that were made unique (`RGS5`, `RGS5.1`) are matched the
  same way; this relies on STARsolo and Cell Ranger listing genes in the
  same (GTF) order.
- The R image is pinned by digest, and provenance records the full
  reference. The STAR and bamtofastq images (biocontainers) are still
  referenced by tag only: a re-pushed biocontainers tag would change what a
  re-run uses.
- The R steps' memory and time were measured on synthetic inputs at the
  real dimensions (`tests/scale/README.md`), not on the real data; compare
  with `peak_rss` in `pipeline_info/trace.tsv` after the first real run.

## Data

Riemondy KA, *et al.* Neoplastic and immune single-cell transcriptomics
define subgroup-specific intra-tumoral heterogeneity of childhood
medulloblastoma. *Neuro-Oncology* 24(2):273-286, 2022. Raw data:
SRA PRJNA649773.

## Acknowledgements

`conf/stjude.config` is adapted from the
[nf-core institutional config for St. Jude](https://github.com/nf-core/configs/blob/master/conf/stjude.config) (MIT).
