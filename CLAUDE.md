# CLAUDE.md: working on mb-scrnaseq

Read this first in every session, then `docs/DEVELOPMENT_PLAN.md` for what
to do next. The plan's status table says where the project is.

## What this is

A Nextflow pipeline that re-quantifies childhood medulloblastoma 10x
scRNA-seq (GEO GSE155446) from the deposited Cell Ranger BAMs and tests,
against pre-registered criteria, whether the result reproduces the
published count matrix (phase 1). Phase 2 is not started: the main
single-cell analysis of all 28 tumors in R (Seurat + Bioconductor), plus a
Python module (scvi-tools, scikit-learn) for two questions, an integration
comparison and a tumor-subgroup classifier. It is a portfolio project for Chia-Lung Chuang (GitHub `ChiaLungChuang`), so
honesty about what has and has not been run matters as much as the code.

## Commands

Run all four before every push; CI runs the same.

```bash
Rscript tests/run_tests.R                          # R unit + accuracy tests (lib/R)
bash tests/shell/run.sh                            # script tests + test-data consistency (needs samtools, python3, nextflow)
nextflow lint . && shellcheck -x bin/*.sh scripts/*.sh tests/shell/*.sh tests/scale/*.sh
nf-test test --profile +docker                     # pipeline on both synthetic datasets
```

From M5, the Python module adds `pytest python/tests` and `ruff check python`
(not there yet).

Without Docker, `--profile +local_tools` runs with tools on `PATH`
(STAR 2.7.11b, bamtofastq 1.4.1, samtools, R 4.4 with Matrix, data.table,
jsonlite, yaml, scuttle). Versions: Nextflow >= 25.10 (CI uses 26.10.4),
nf-test 0.9.5, R 4.4.3 / Bioconductor 3.20 (`containers/r-analysis`).

Regenerate synthetic data after changing a generator, then commit it:
`python3 tests/data/make_test_data.py` and
`python3 tests/data/make_accuracy_data.py`. `tests/shell/test_test_data.sh`
fails if committed data and generators disagree.

## Rules that are not negotiable

1. **Blinding.** Never print, summarise, plot or compute anything from the
   *values* of the published data (`GSE155446_human_raw_counts.csv.gz`,
   `GSE155446_human_cell_metadata.csv.gz`) before the criteria that judge
   them are locked. Names, dimensions and labels are fine. Phase 1 criteria
   lock before the first concordance run; the cell metadata's values stay
   blind until the phase 2 criteria are locked.
2. **Locked thresholds do not move.** `preregistration/*.yaml` with
   `status: locked` is never edited to make a result pass. If a criterion is
   wrong in principle, propose a change to Chia-Lung, and report the old
   and new verdicts side by side.
3. **Every test must be able to fail.** For each new or changed test,
   plant the bug it targets, show the test fails, revert. Say in the PR
   which bug you planted and what failed. A test that passes against
   broken code is worse than none.
4. **Docs claim only what was run.** README status, the runbook and any CV
   text must not say "run on the cluster", "containerised and tested",
   "reproduces" etc. until it has happened, with the evidence committed.
5. **Propose, don't pick.** At the decision points in the plan, and for
   any choice that changes scientific results (parameters, filters,
   references, criteria), write the options with trade-offs and stop for
   Chia-Lung's decision. Engineering choices with no effect on results
   you may make; note them in the PR.
6. **No data in git** beyond the tiny synthetic sets in `tests/data`. No
   BAMs, FASTQs, real count matrices or real per-cell outputs. Curated
   small result files (verdict, summary table, provenance, trace) go in
   `docs/results/` once real runs exist.
7. **The cluster is Chia-Lung's.** Do not ssh to or run anything on the St.
   Jude HPC cluster unless he says HPCF/IT policy allows it. Default: you
   prepare commands and fixes, he runs them and pastes logs back.

## How to work

- One branch and one PR per task; CI green before merge; small commits
  with messages that say why. Use `gh` for PRs and for CI logs
  (`gh run view --log-failed`).
- Tests first for behaviour changes, including bug fixes: reproduce the
  bug in a test, then fix it.
- When a PR is done, update the status table in `docs/DEVELOPMENT_PLAN.md`
  and say what is next and what needs Chia-Lung.
- End-of-session report: what changed, test counts before/after, planted
  bugs and results, open decisions. Keep it short.

## Things learned the hard way (keep them true)

- Command-line booleans arrive as strings: `--flag false` is the string
  "false", which Groovy treats as true. Parse with `boolParam()` in main.nf.
- Strict syntax (`nextflow lint`): no `import`, use fully qualified class
  names (`new org.yaml.snakeyaml.Yaml()`); `onComplete:` lives in the entry
  workflow; helper functions are top-level `def`s.
- Nextflow stages input folders as symlinks: use `find -L`.
- `-resume` does not see changes in `bin/` or inside a staged folder such
  as `lib/R`: each step that runs that code takes a SHA-256 fingerprint of
  only the code it runs as a `val` input (see `codeFingerprint()`). Keep
  it scoped: a fingerprint that is too wide re-runs STARsolo for every R
  edit.
- STARsolo's raw matrix has a column for every whitelist barcode (6.8M for
  v3). "Present in raw" means nothing; "at least one UMI" is the test.
- bamtofastq names the cDNA read R2 or R3 depending on the BAM header.
- `samtools view -b` embeds the command line (with local paths) in @PG;
  generators use `--no-PG`.
- data.table < 1.15 reads `fread(text = x)`, when `x` is one string without
  a newline, as a file name, or as a shell command if it has a space. Parse
  lines with `fread_lines()` (lib/R/io.R), which always adds the newline.
- `methods::new("dgCMatrix", ...)` fails in a fresh session until something
  has loaded the Matrix namespace; the unit tests always have, `bin/`
  scripts may not. Look the class up with `asNamespace("Matrix")`.
- `r_container` is pinned as `name:tag@sha256:...`. `docker tag` cannot take
  a digest, so CI tests the image built from the commit through the
  `local_r_image` profile, and checks the pinned one in `pinned-r-image`.
  Releases: bump without a digest, tag, then pin (docs/RUNBOOK.md, 0).

## Settled decisions (do not reopen without asking)

- Nextflow DSL2. Phase 1 stays in R: do not port it to Python (decided
  31 Oct 2026; its language is infrastructure and barely shows).
- Phase 2 splits by language, each used where it is strongest. The main
  analysis is R (Seurat + Bioconductor): it matches the published analyses
  of this dataset and the R-leaning roles the project targets. A Python
  module on the same data answers two questions only: scVI versus Harmony
  integration, and whether tumor subgroup can be predicted from
  malignant-cell expression (scikit-learn, validated by patient). Do not
  duplicate the R analysis in Python.
- Phase 1 re-quantifies 8 of the 28 tumors from the BAMs; phase 2 analyses
  all 28 from the published matrix, after phase 1.
- STARsolo 2.7.11b with the Cell Ranger 3.x options (`cellranger_compat
  cr3`) and the GRCh38-3.0.0 reference: the BAM headers checked name
  Cell Ranger 3.x (STAR 2.5.1b, refdata-cellranger-GRCh38-3.0.0).
- BAMs come from ENA (size + MD5 verified); SRA prefetch is the fallback.
- Criterion C1 counts a published barcode as recovered only with >= 1 UMI.
- Samples are mapped to published labels by
  `docs/proposals/gse155446_runs.tsv` (all 30 human runs, mapped by barcode
  overlap), never by the SRA/GEO title or a read-group name alone:
  SRR12354571 (GEO "934-repeat") is published label 943, most likely a
  second library of tumor 934; "966-recurrence" is 966-2; read groups
  `Foreman_1` to `Foreman_4` carry no label. 30 libraries, 28 tumors.
- The 8 phase 1 tumors (accepted 2 Nov 2026) are in
  `docs/results/phase1/samples.tsv`, checked against the run table by
  `tests/shell/test_phase1_samples.sh`. Label 1355 (SRR12354589, 6.4 GB) is
  the only v2 run and the pilot.
