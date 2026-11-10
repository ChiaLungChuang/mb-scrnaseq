# Development plan

The order of work from here to a finished project, written so Claude Code
can follow it a milestone at a time. Rules for *how* to work (blinding,
locked thresholds, testing, honesty) are in [`CLAUDE.md`](../CLAUDE.md);
this file says *what* to do and when to stop.

Markers:

- **[you]**: only Chia-Lung can do this (accounts, the cluster, scientific
  or publishing decisions). Claude Code prepares and then waits.
- **[decide]**: a decision point. Claude Code writes the options with
  trade-offs and a recommendation, then stops until Chia-Lung chooses.

## Status

Claude Code updates this table at the end of every PR.

| Milestone | State | Evidence |
|---|---|---|
| Pre-M0: phase 1 code and tests | done | tests pass locally (tools on `PATH`) and in CI |
| M0: repository on GitHub, CI green, image published | done | Repository public (step 7; GitHub reports `visibility: public` on 3 Nov 2026), README written for visitors with the CI badge. CI green on `main` (latest checked: [41547964441](https://github.com/ChiaLungChuang/mb-scrnaseq/actions/runs/41547964441)): R, script, lint and pipeline tests in Docker, synthetic data only. `v0.1.0` published `ghcr.io/chialungchuang/mb-scrnaseq-r:0.1.0` ([run 41257603921](https://github.com/ChiaLungChuang/mb-scrnaseq/actions/runs/41257603921)); `r_container` pins it by digest (`sha256:9ce6480a…`); the `pinned-r-image` job pulls it with no login and runs the R tests inside it. |
| M1: cluster smoke test | prepared, not run | `scripts/hpcf_preflight.sh` ([#11](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/11); compute-node session, ENA and download hosts, lsf.conf memory settings and an opt-in job submission test in [#17](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/17)); `conf/stjude.config` compared with nf-core's ([#9](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/9)) and aligned with HPCF's policy, with its LSF directives tested against a stand-in `bsub` ([#16](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/16)); runbook follows HPCF's policy (nothing beyond `cd` and `bsub` on login nodes; two paths for the head process). Nothing has run on the cluster. Next **[you]**: HPCF's answers to the questions in runbook section 1 (item 6), then section 3 (preflight with `--submit-test` in an interactive session, then `nf-test test --profile +stjude`) |
| M2: real inputs and phase 1 pre-registration | in progress | 8 tumors chosen (2 Nov 2026, [#10](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/10)) and recorded in `docs/results/phase1/samples.tsv` (steps 2 and the record of 5). Not done: inputs and BAMs downloaded **[you]**, `bams_summary.tsv` read, criteria locked **[you]** |
| M3: phase 1 run on 8 tumors | prepared, not run | results report `bin/phase1_report.R` ([#13](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/13)), tested on synthetic runs only; raw trace with CPU and memory fields. The R steps rehearsed on synthetic inputs at real dimensions (`tests/scale/README.md`): peaks of at most 4.8 GB (`PREP_PUBLISHED`, after removing two copies) and 5.8 GB (`CONCORDANCE`), requests resized; peak disk for the 8 tumors estimated at 350-445 GB (runbook section 9). Not run on real data |
| M4: phase 2 design and pre-registration | in progress | design drafted for decision in [draft PR #4](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/4); nothing decided, no criteria drafted |
| M5: phase 2 implementation (R analysis + Python module) and run | not started | known-answer data generator drafted ([draft PR #14](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/14)); waits on the M4 decisions |
| M6: write-up | not started | |

What exists before M0: the phase 1 pipeline (BAM → bamtofastq → read
check → STARsolo → QC → concordance), the St. Jude LSF/Singularity
profile, the R container recipe, CI configuration, the runbook, and the
three test layers below. Before M0 none of it had run in CI or in a
container. Since 31 Oct 2026 it has, on synthetic data; none of it has
run on the cluster or on real data.

## Testing rules

Every change ships with tests in the layer it touches, and every test is
shown to fail when the bug it targets is planted (CLAUDE.md rule 3).

| Layer | Where | Kind of test |
|---|---|---|
| R library | `tests/testthat/` | unit tests for chosen cases; accuracy tests against an independent implementation on random inputs (`test-accuracy.R`) |
| Scripts | `tests/shell/` | case tables on hand-built inputs; network and containers replaced by stand-ins |
| Pipeline | `tests/*.nf.test` | whole runs on synthetic data whose correct output is known exactly |
| Python module (from M5) | `python/tests/` | pytest unit and property tests (hypothesis) against independent implementations; known-answer data |

Known-answer data is generated, never hand-edited: `tests/data/make_test_data.py`
(basic) and `tests/data/make_accuracy_data.py` (every read type that can be
miscounted, v2 and v3). New analysis steps in phase 2 need the same: a
generator with planted structure (cell types, copy-number changes) and a
test that the step recovers it.

Mutation checks done so far (each made the tests fail):

| Planted bug | Caught by |
|---|---|
| empty whitelist columns counted as recovered (the old C1 bug) | `test-accuracy.R`, `test-concordance.R`, `main.nf.test` swapped labels |
| per-cell Pearson off by one; UMI ratio inverted; duplicate symbols summed when they should be made unique | `test-accuracy.R` |
| rows lost at CSV chunk boundaries; QC threshold on the wrong scale | `test-accuracy.R` |
| `_R1_` replaced in the folder name; cDNA always assumed R2 | `test_check_reads.sh` |
| no `samtools quickcheck`; index not excluded from the prefetch rename | `test_fetch_inputs.sh` |
| STARsolo with Exact UMI dedup, no multi-gene filter, no barcode correction, or unstranded counting | `accuracy.nf.test` |
| scuttle not loadable, so testthat skips the QC tests instead of failing them | `tests/run_tests.R` (any skip fails the run) |
| `fread(text =)` given a bare line, which data.table < 1.15 reads as a file name or a command (header and chunks; chunks alone) | `test-io.R` (a strict `fread` stand-in, so it fails on any data.table version) |
| version check reads the digest as the image tag; digest format unchecked; release tag ignored; an already pinned version published again; a profile override read as the pipeline's image | `test_check_versions.sh` |
| a package fails to install but the install step passes (simulated: install functions that only warn) | `tests/install_r_packages.R` (loads every package) |
| R tasks run in the published image instead of the commit's build, or not in the local image at all (simulated task scripts) | CI step "The R steps ran in the image built from this commit" |
| preflight: a missing queue only warns; versions compared as text; Java 11 accepted; a probe file written; queue list hard-coded; 1 GB floor; HTTP 401 counted as reachable; module load not checked; the submit test run without `--submit-test`, EXIT counted as PASS, no polling, no `bkill` on timeout, output files written; ENA's answer not checked; GET instead of HEAD on the download hosts; no connection counted as reachable; LSF_ENVDIR ignored; lsf.conf unfiltered; the session check always passing | `test_hpcf_preflight.sh` (stand-ins under `env -i`, including a changes-nothing check) |
| St. Jude profile: `perTaskReserve` false or removed (rusage not per core); `perJobMemLimit` true (-M per job); large_mem from 512 GB instead of above 1 TB; short queue only up to 10 minutes | `test_stjude_lsf.sh` (Nextflow with a stand-in `bsub`, with and without an lsf.conf) |
| pre-registration gate: a run with `--allow_unlocked_criteria true` and locked, pushed criteria marked pre-registered (the old behaviour); the test-run note added whatever the flag; pre-registered whatever the problems | `main.nf.test` (criteria in a scratch git repository pushed to a scratch remote; with the gate enforced, the run passes on a clean checkout) |
| standalone server: no per-task cap (`resourceLimits`); no CPU or memory cap on what runs at once (local executor); the Apptainer profile leaving Singularity on | `test_server_config.sh` (a tiny pipeline under the server config; task pairs that fit the share's memory but not its cores, and the reverse) |
| server launch check: minimum memory, machine memory, cores or disk not checked; a missing share not caught; the check never called | `server.nf.test` |
| preflight `--local`: memory threshold not read from `nextflow.config`; the cluster checks still run; Docker without its daemon counted as a runtime; `--submit-test` allowed with it; the cache variable always Singularity's | `test_hpcf_preflight.sh` |
| R thread limits: a step without them (data.table took 2 threads on a 1-CPU task); one variable missing; a limit hard-coded instead of `${task.cpus}`; an R script that does not log; `log_threads()` reporting a fixed or wrong value | `main.nf.test` (each R task's logged threads equal its CPUs), `test_r_threads.sh`, `test-io.R` |
| published-matrix reader: TRUE/FALSE read as 1/0 (the old behaviour); the NA check removed, so a missing value becomes 0; the sparse class created before Matrix is loaded (fresh session); `rbind_csc` with no row offset, positions off by one, counts not advanced, rows offset by one, or no column check | `test-io.R`; the offsets and positions also by `test-accuracy.R` |
| phase 1 sample record: wrong label, chemistry or MD5; both libraries of patient 934; pilot not first; a sample dropped (and, in the checker: patient folding, MD5 comparison or pilot rule removed) | `test_phase1_samples.sh` |
| results report: %cpu not divided by 100; cached tasks counted; criteria refs unfiltered; hash consistency check removed; PASS/FAIL swapped; powers of 1000; "any" for "all"; human-formatted trace read as NA | `test-report.R`; PASS/FAIL swap and dropped tasks also by `check_phase1_report.R` on the synthetic runs |

## M0. Repository on GitHub, CI green, image published

Goal: the code lives on GitHub and CI proves it works in containers.

1. **[you] [decide]** Repository visibility. Recommendation: create it
   private, make it public once CI is green, so the first public commit
   has a passing badge. Create it on github.com, or let Claude Code run
   `gh repo create ChiaLungChuang/mb-scrnaseq --private --source . --remote origin --push`.
2. Push `main`. CI has never run, so expect failures. Likely ones:
   R package installation time or Bioconductor version on the runner; the
   R image build; Docker file permissions with `-u $(id -u):$(id -g)`;
   biocontainers images missing a tool a script uses; ShellCheck or
   samtools version differences; actions that warn about deprecated
   Node.js versions (update them to current major versions).
3. Fix each failure in its own PR, with the cause in the description.
   Never make CI green by weakening a test (looser tolerance, skipped
   case, removed assertion). If a test is wrong, say why in the PR.
4. **[you]** When `main` is green: `git tag v0.1.0 && git push --tags`.
   The `publish-image` job pushes `ghcr.io/chialungchuang/mb-scrnaseq-r:0.1.0`.
5. **[you]** On GitHub: Packages → `mb-scrnaseq-r` → Package settings →
   visibility Public (Singularity on the cluster pulls without logging in).
6. Claude Code checks that an anonymous pull works, then in a follow-up PR
   pins `r_container` by digest (`...:0.1.0@sha256:...`) and adjusts the
   CI version check to compare the tag part. That closes the "tags, not
   digests" limitation in the README.
7. **[you] [decide]** Make the repository public. Done (public by 3 Nov
   2026).

Done when: every CI job is green on `main`, the image is public and
pinned, and the README status says containers and CI pass (on synthetic
data), and still says nothing has run on real data or the cluster.

## M1. Cluster smoke test

Goal: the same tests pass on St. Jude HPCF (LSF + Singularity).

1. **[you]** HPCF account (login host `hpc.stjude.org`), a project or
   scratch path, and HPCF's answers to the questions in runbook section 1
   (compute-node internet, job submission from jobs, `tmux` on login
   nodes, a queue for the head job, how `-M` is enforced).
2. **[you] [decide]** Whether Claude Code may run commands on the cluster.
   Check HPCF/IT policy first. Default: no; you run, paste logs back.
3. **[you]** In an interactive session (login nodes are only for `cd`
   and `bsub`): runbook section 2 (micromamba, Nextflow 26.10.4, nf-test
   0.9.5, Singularity cache, clone), then section 3:
   `scripts/hpcf_preflight.sh --submit-test`, then
   `nf-test test --profile +stjude`. This runs both synthetic test files
   through LSF and Singularity with the same assertions as CI, including
   the exact-count accuracy test.
4. For each failure, paste `.nextflow.log`, the failing task's
   `.command.err` and `bhist -l <job>`. Claude Code fixes
   `conf/stjude.config` (queues, memory reservation, binds, `TMPDIR`,
   pull timeouts) with the reason in the PR.
5. Commit `docs/results/m1/trace.tsv` from a passing run as evidence.

Done when: `nf-test test --profile +stjude` passes on the cluster and the
README says the cluster profile is tested on synthetic data.

## M2. Real inputs and phase 1 pre-registration

Goal: inputs downloaded and verified, samples chosen, criteria locked and
public before any real concordance result exists.

1. **[you]** Runbook section 4: `geo`, `whitelists`, `reference
   GRCh38-3.0.0`.
2. **[decide]** Which 8 tumors. Claude Code proposes a set that spans the
   subgroups (cohort: 1 WNT, 9 SHH, 7 Group 3, 11 Group 4), using the SRA
   run table and the paper, not the GEO cell metadata file (blinded).
   Include SRR12354589 (tumor 1355, v2, 6.4 GB) as the pilot.
3. **[you]** `DRY_RUN=1 scripts/fetch_inputs.sh bams ...` to check sizes,
   then download.
4. Claude Code reads `bams_summary.tsv`. If any BAM names a reference
   other than GRCh38-3.0.0 or a STAR version other than 2.5.1b,
   **[decide]** before going on.
5. Build the samplesheet from `bams_summary.tsv` and
   `published_labels.tsv` (labels from names only). Commit the record as
   `docs/results/phase1/samples.tsv` (sample, SRR, label, chemistry,
   subgroup); the runnable samplesheet with cluster paths stays on the
   cluster.
6. **[you] [decide]** Review `preregistration/acceptance_criteria.yaml`.
   Change thresholds now or never. Set `status: locked`, commit, push, tag
   `prereg-phase1`, and make a GitHub release from the tag. In the commit
   message, disclose that the code review read the published matrix's
   dimensions, labels, gene names and barcodes, but no count values.

Done when: the locked criteria are on GitHub with a release, the sample
record is committed, and every BAM passed its checks.

## M3. Phase 1 run

Goal: the verdict for 8 tumors, reported whatever it is.

1. **[you]** Pilot: the 1355 sample alone, full pipeline (runbook 9).
2. Claude Code reads `trace.tsv` (`peak_rss`, `realtime`) and right-sizes
   the resource labels. This is engineering, not science: no effect on
   results, so no decision needed; note it in the PR.
3. If the pilot fails a criterion: follow runbook 10 to find the cause.
   Thresholds do not move. A code change is allowed only for a
   demonstrable bug, with a failing test written first. Every such change
   is logged in `docs/results/phase1/log.md` with the verdict before and
   after it. That log is the record that results were not tuned.
4. **[you]** All 8 samples, with `-resume`.
5. Commit to `docs/results/phase1/`: `concordance_summary.verdict.json`,
   `concordance_summary.tsv`, `run_provenance.json`,
   `software_versions.yml`, `trace.tsv`, the QC threshold JSONs, and a
   `README.md`: commit, criteria hash, dates, the verdict per sample and
   criterion, the explanation of every failure, and compute used (CPU
   hours, peak memory) from the trace.
6. Update the README status. **[you]** Fill in the CV bullet from these
   numbers. Tag `v0.2.0`.

Done when: results and report are committed and the README states the
verdict.

## M4. Phase 2 design and pre-registration

Goal: an agreed analysis plan for all 28 tumors (published matrix) and
locked criteria, before any value in the cell metadata is read.

Already decided (31 Oct 2026): the main analysis is R (Seurat +
Bioconductor); a Python module answers two questions on the same data
(below); phase 1 is not ported. Everything else is for the design doc.

1. Claude Code writes `docs/phase2_design.md`. For each step: options,
   trade-offs, a recommendation, how it will be tested with known-answer
   data, and a compute estimate.

   R, main analysis:
   - QC: the published matrix is already filtered by the authors;
     re-filter with the phase 1 MAD rules or not.
   - Normalisation: LogNormalize or SCTransform (v2).
   - Annotation: marker-based or reference-based (e.g. SingleR) and which
     reference.
   - Copy-number inference to separate malignant cells: inferCNV,
     copykat, or numbat (which needs allele counts from BAMs, so only the
     8 re-quantified tumors); which cells serve as the normal reference.
   - Tumor programs: per-tumor NMF and meta-programs, or not.

   Python module:
   - **Integration comparison.** Harmony (run in the R analysis) against
     scVI (scvi-tools), with no integration as the baseline. Use the cells
     where integrating across patients makes sense: non-malignant cells
     (immune, stromal). Malignant cells are patient-specific, so
     integrating them can erase real differences; the design says how
     they are handled. Metrics for batch mixing (e.g. iLISI, graph
     connectivity) and for keeping cell types apart (e.g. ARI or NMI
     against the published labels, cLISI), for example from
     scib-metrics. A decision rule fixed in advance, e.g. prefer scVI only
     if mixing improves and label ARI falls by no more than an amount set
     now. scVI training is random: fix seeds and report the spread over
     several.
   - **Subgroup classifier.** Predict each tumor's subgroup from its
     malignant cells' expression with scikit-learn. Validate by patient
     (leave-one-patient-out or grouped folds); never split one patient's
     cells between training and test. There is one WNT tumor, so WNT cannot
     be both learned and tested: the classifier is SHH / Group 3 /
     Group 4 on 27 tumors, and WNT is reported separately. Baselines: the
     majority class and a known-marker rule. Chance level from permuting
     subgroup labels across patients. Also run the leaky version (cells
     split at random) once and report it next to the real one, to show
     how much leakage inflates accuracy.
   - **Language boundary:** how data moves from R to Python (Matrix Market
     plus TSV, or `.h5ad`), with a round-trip test that the matrix and
     metadata arrive identical.
   - **Python layout:** a `python/` package (`pyproject.toml`), `bin/*.py`
     entry points, pytest with hypothesis, ruff, a second image
     `mb-scrnaseq-py` published and pinned like the R one, a CI job.

   Wiring: new processes and a subworkflow in this pipeline, or a
   separate entry workflow.
2. **[you] [decide]** Choose.
3. Claude Code drafts `preregistration/phase2_criteria.yaml` from the
   published cell metadata's column **names** only: agreement of
   annotation and malignant calls with the published labels, the
   integration decision rule, and the classifier's claim rule (e.g. say
   "subgroup is predictable from malignant-cell expression" only if
   patient-level balanced accuracy beats both the permutation null at
   p < 0.05 and the marker baseline). Each threshold has its rationale.
4. **[you] [decide]** Review, lock, push, release `prereg-phase2`. Only
   after this may any value in the cell metadata be read.

Done when: the design is approved and the phase 2 criteria are locked and
released.

## M5. Phase 2 implementation (R analysis + Python module) and run

1. Known-answer data first. Generators with planted cell types, planted
   patient batch effects and planted copy-number gains and losses (for
   the R steps), planted subgroup signal (for the classifier), and a
   negative control with no subgroup signal. Tests that each step
   recovers what was planted, written before the step. The negative
   control is the leakage test: the classifier must stay at chance on it,
   and switching to random cell-level splits must make that test fail.
2. R analysis: lib/R functions, `bin/*.R` entry points, Nextflow
   processes. Extend the R image (Seurat, Harmony, the chosen annotation
   and CNV packages), versions pinned; release and pin its digest.
3. Python module: the `python/` package, `bin/*.py` entry points,
   Nextflow processes, the `mb-scrnaseq-py` image (scvi-tools,
   scikit-learn, anndata/scanpy, scib-metrics), pinned. CI job for pytest
   and ruff. The language-boundary round-trip test runs in CI.
4. **[you]** Run on the cluster. scVI on CPU is probably enough for about
   40,000 cells; the design doc costs the gpu queue as an option. Commit
   curated outputs and a report to `docs/results/phase2/` as in M3, with
   figures as PDF and TIFF.

Done when: CI is green for both languages, and results and the phase 2
verdict are committed.

## M6. Write-up

1. Final README: what the project shows, with links to both result
   reports and the pre-registration releases.
2. A short report (methods, what reproduced, what did not and why).
3. **[you]** Final CV bullets from the committed numbers: the R analysis
   and the Python module each get their own line, with what was measured.

## Parked (do not start without asking)

- Porting phase 1 to Python: decided against (31 Oct 2026).
- The next Python-first project, aimed at Python-required roles, belongs
  in its own repository, not here (candidate: a flow-cytometry classifier
  on public FlowCAP-II data with scikit-learn/XGBoost and patient-level
  validation).
- Re-quantifying all 28 tumors from the BAMs: revisit after M3, by cost.
- Running Cell Ranger itself: not planned (STARsolo is the method).
- Docker on the cluster: not available; Singularity only.

## Known risks

- `CONCORDANCE` memory on a full-size sample is unmeasured (M3 step 2).
- STAR 2.7.11b vs the 2.5.1b inside Cell Ranger 3.x: alignments can
  differ slightly; the criteria are meant to tolerate that.
- Gene matching is by symbol and relies on STARsolo and Cell Ranger
  listing genes in the same GTF order (true in a test with STAR 2.7.11b;
  real data will show it through `gene_match_fraction`).
- The classifier has 27 tumors for 3 classes, and the Group 3 / Group 4
  boundary is biologically blurred: expect errors there and a wide
  confidence interval. The point is a correctly validated estimate, not a
  high number.
- Anonymous GHCR pulls from the cluster could be rate-limited; the runbook
  sets a Singularity cache so each image is pulled once.
