# mb-scrnaseq

[![CI](https://github.com/ChiaLungChuang/mb-scrnaseq/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ChiaLungChuang/mb-scrnaseq/actions/workflows/ci.yml?query=branch%3Amain)

A Nextflow pipeline that re-quantifies childhood medulloblastoma
single-cell RNA-seq (GEO
[GSE155446](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155446):
28 tumors, 10x Genomics 3') from the authors' deposited Cell Ranger BAMs,
and tests, against criteria written down before any result exists, whether
the result reproduces the published count matrix.

**Status: built and tested on synthetic data only.** CI runs the test
suites on every push to `main` and every pull request. The acceptance
criteria are locked. Nothing has run on St. Jude compute or on the real
data yet. [Details](#status).

## What it shows

- **A pre-registered comparison.** Five thresholds (C1-C4: barcode
  recovery, per-cell UMI totals, pseudobulk and per-cell correlation) are in
  [`preregistration/acceptance_criteria.yaml`](preregistration/acceptance_criteria.yaml).
  On real data the pipeline refuses to compare until that file is locked,
  committed and pushed, and every run records the file's hash and commit,
  so the order of events can be checked ([below](#pre-registration)).
- **Tests with known answers.** The whole pipeline runs on two synthetic
  datasets whose correct output is known exactly. In one, every cell x gene
  count must equal Cell Ranger's counting rules, on reads built to trip
  them (duplicates, barcode and UMI errors, multi-mapping, antisense and
  intronic reads); the other checks UMI totals, QC flags and every refusal.
  Tests are checked by planting the bug each one targets and confirming it
  fails ([the list](docs/DEVELOPMENT_PLAN.md#testing-rules)).
- **Containers and CI.** Every step runs in a container at a fixed version.
  The R image is built from [`containers/r-analysis`](containers/r-analysis),
  published on
  [GHCR](https://github.com/ChiaLungChuang/mb-scrnaseq/pkgs/container/mb-scrnaseq-r)
  and pinned by digest. [CI](.github/workflows/ci.yml) runs the R, script,
  lint and whole-pipeline tests in Docker, and the R tests inside the
  published image.
- **Two routes on St. Jude compute**: a cluster profile (LSF +
  Singularity) with a test of the LSF directives it writes, and a
  standalone-server profile (local executor) that caps every task at the
  machine share it is given and refuses to start without enough memory or
  disk. One preflight script serves both. Prepared, not yet run.

## Status

As of 10 Nov 2026:

| Part | State | Evidence |
|---|---|---|
| Phase 1 pipeline: BAM -> FASTQ -> STARsolo -> QC -> concordance | built; passes on synthetic data | [CI runs on `main`](https://github.com/ChiaLungChuang/mb-scrnaseq/actions/workflows/ci.yml?query=branch%3Amain) |
| R image `ghcr.io/chialungchuang/mb-scrnaseq-r:0.1.0` | published, pullable without a login, pinned by digest | [release run](https://github.com/ChiaLungChuang/mb-scrnaseq/actions/runs/41257603921); the `pinned-r-image` CI job |
| St. Jude compute: the Research HPC cluster (LSF + Singularity) or a department server (no scheduler) | both prepared, not run; waiting on the HPCF account and answers, and on the department's reply about its machines | [runbook](docs/RUNBOOK.md#two-routes) |
| Real data | 8 of the 28 tumors chosen; criteria locked (10 Nov 2026); inputs not yet downloaded | [`samples.tsv`](docs/results/phase1/samples.tsv), [proposal](docs/proposals/phase1_samples.md), [lock](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/22) |
| Phase 2: all 28 tumors in R (Seurat), plus scVI vs Harmony and a subgroup classifier in Python | designed, not started | [draft PR #4](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/4) |

The [development plan](docs/DEVELOPMENT_PLAN.md) keeps the full status
table and what comes next.

## Where to look

- [`docs/DEVELOPMENT_PLAN.md`](docs/DEVELOPMENT_PLAN.md): milestones,
  status, testing rules and the planted-bug list.
- [`docs/RUNBOOK.md`](docs/RUNBOOK.md): the real run on St. Jude compute
  (the cluster, or a department server), step by step, with memory and
  disk estimates.
- [`preregistration/`](preregistration/): the acceptance criteria.
- [`docs/proposals/`](docs/proposals/): why these 8 tumors, and all 30
  sequencing runs mapped to the published sample labels by barcode overlap.
- [`main.nf`](main.nf), [`modules/local/`](modules/local/),
  [`lib/R/`](lib/R/), [`bin/`](bin/): the pipeline and its R code.
- [`tests/`](tests/): the synthetic datasets and their generators, the
  test suites, and a full-size rehearsal of the R steps
  ([`tests/scale/`](tests/scale/README.md)).

## Why re-quantify instead of using the published matrix

A published matrix is one run of one pipeline. Re-deriving it from the
reads checks that the data, the sample labels and the processing are what
the paper says, before any biology is built on them: a wrong chemistry,
whitelist, reference or sample mapping shows up as a failed criterion
rather than as plausible but wrong results later.

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

All 30 human BAM headers of GSE155446 (their first 4 MB) name Cell
Ranger 3.x (STAR 2.5.1b) and `refdata-cellranger-GRCh38-3.0.0`, so the
dataset is run with `--cellranger_compat cr3` and that reference;
`scripts/fetch_inputs.sh bams` reports both for every BAM it downloads.

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

Synthetic test (about a minute; CI runs it through nf-test):

```bash
nextflow run . -profile test,docker        # or test,singularity
```

The `singularity`, `apptainer` and `stjude` profiles have not been run yet;
`server` runs in CI with Docker. Real data on St. Jude compute (the
cluster, or a department server): [docs/RUNBOOK.md](docs/RUNBOOK.md).

Outputs land in `--outdir`: `starsolo/<sample>/` (raw and filtered count
matrices, `Summary.csv`, STAR log), `qc/`, `concordance/`, `reads_check/`,
`reference/` (the STAR index, when built) and `pipeline_info/` (provenance,
software versions, Nextflow report, timeline and trace).

## Tests

Three layers. Each test is checked by planting the bug it targets and
confirming that it fails ([the list](docs/DEVELOPMENT_PLAN.md#testing-rules)).

- **R library** (`Rscript tests/run_tests.R`): unit tests, plus accuracy
  tests that compare each function with an independent implementation on
  random inputs (the concordance metrics with a dense base-R version, QC
  thresholds with `scuttle` and by hand, the CSV and Matrix Market readers
  with whole-file reads). A skipped test fails the run.
- **Scripts** (`bash tests/shell/run.sh`): every script with stand-ins for
  the network, containers and LSF. `check_reads.sh` on hand-built FASTQ
  folders, `fetch_inputs.sh` (MD5, size, truncation, 10x headers, SRA
  fallback), `check_versions.sh`, `hpcf_preflight.sh` (including that it
  changes nothing unless asked to submit its one test job, and its
  `--local` mode for a server), the LSF directives the St. Jude profile
  makes Nextflow write, the caps the standalone-server profile puts on
  each task and on the machine, the phase 1 sample
  record against the run table, and the committed test data against its
  generators.
- **Pipeline** (`nf-test test --profile +docker`), on two synthetic
  datasets with known answers:
  - `tests/main.nf.test`: every cell's UMI total matches the truth, QC flags
    exactly the planted low-quality cells, both read layouts work, swapped
    labels fail barcode recovery, and wrong chemistry, unlocked criteria, a
    missing `--cellranger_compat` and duplicate samples are refused.
  - `tests/accuracy.nf.test`: every cell x gene count equals Cell Ranger's
    rules, for v3 and v2 chemistry, on reads built to trip them: PCR
    duplicates, barcode and UMI errors, UMIs shared between genes, spliced,
    minus-strand, antisense, intronic, intergenic and multi-mapping reads,
    and two genes sharing a symbol.
  - The results report (`bin/phase1_report.R`) is checked on both runs'
    outputs, and every R step's logged thread count equals its CPUs.
  - `tests/server.nf.test`: the standalone-server profile's launch checks
    (too little memory, more than the machine has, no disk for one tumor).
    CI also runs both suites under `-profile server`.

[CI](.github/workflows/ci.yml) runs all of this on every push to `main` and
every pull request, with the pipeline tests on the R image built from that
commit, plus `nextflow lint`, ShellCheck, a version check (release tag,
`manifest.version` and the R image tag agree) and a job that pulls the
pinned R image without logging in and runs the R tests inside it.

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

## License and acknowledgements

MIT ([LICENSE](LICENSE)). `conf/stjude.config` is adapted from the
[nf-core institutional config for St. Jude](https://github.com/nf-core/configs/blob/master/conf/stjude.config) (MIT).
