# Runbook: first real run on the St. Jude cluster

This is the order to do things in, from no account to a concordance verdict.
Steps marked **(you)** need your accounts, judgement or St. Jude-specific
details that are not public; everything else is a command to run.

## How the cluster fits together (two minutes)

- You `ssh` into a **login node**. It is shared, so you do not compute on it;
  you only launch things from it.
- Nextflow runs on the login node as a lightweight **head process**. For every
  step it writes a small job script and submits it to **LSF** (`bsub`), the
  scheduler that places jobs on **compute nodes** with the CPUs, memory and
  time each step asked for (see `conf/base.config`). The `stjude` profile picks
  the LSF queue for you: `short`, `standard`, or `large_mem`.
- Each job runs inside a **Singularity** container, pulled once from the
  image names in the modules. That is how the exact STAR, bamtofastq and R
  versions follow the pipeline onto the cluster. Docker is not available there.
- The head process has to stay alive for the whole run (hours), so start it
  inside `tmux` and detach. If it dies, re-run the same command with
  `-resume` and finished steps are skipped.
- Intermediate files live in `work/`. Put it on project or scratch storage,
  not your home directory: the FASTQs regenerated from 8 BAMs run to roughly
  the size of the BAMs themselves.

## 0. Before the cluster

1. **(you)** Push this repository to GitHub as `ChiaLungChuang/mb-scrnaseq`.
   CI runs the R unit tests, lint, and the full pipeline on the synthetic
   data. Fix anything red before going further.
2. Tag the first release: `git tag v0.1.0 && git push --tags`. CI then builds
   the R image and publishes it as `ghcr.io/chialungchuang/mb-scrnaseq-r:0.1.0`.
3. **(you)** On GitHub, open Packages → `mb-scrnaseq-r` → Package settings and
   set visibility to **Public**, so Singularity can pull it without a login.

`r_container` in `nextflow.config` pins the R image by digest
(`...:0.1.0@sha256:...`), so re-pushing a tag cannot change what a run uses.
For a later release (e.g. `v0.2.0`): in one PR, set `manifest.version` and
the `r_container` tag to the new version **without** a digest; when it is
merged, tag and push (`git tag v0.2.0 && git push --tags`), and CI publishes
the image; then pin the new digest (`docker buildx imagetools inspect
ghcr.io/chialungchuang/mb-scrnaseq-r:0.2.0` prints it) in a follow-up PR.
`scripts/check_versions.sh` (run by CI) refuses a tag that disagrees with
`manifest.version`, a malformed digest, and publishing a version that is
already pinned.

## 1. Account and storage

4. **(you)** Request an HPCF account through IT, if you do not have one.
5. **(you)** Find out from HPCF documentation or your lab: the login hostname,
   and which project or scratch path you can use for large data. Note both;
   below they are `<login-host>` and `<project-dir>`.

## 2. One-time setup on the login node

```bash
ssh <you>@<login-host>

# micromamba in your home directory (small; no admin rights needed)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# open a new shell, then:
micromamba create -y -n nf -c conda-forge -c bioconda nextflow=26.10.4 nf-test=0.9.5
micromamba activate nf
nextflow -version        # should print 26.10.4

# Singularity image cache; there is no shared one on HPCF
mkdir -p <project-dir>/singularity_cache
echo 'export NXF_SINGULARITY_CACHEDIR=<project-dir>/singularity_cache' >> ~/.bashrc
source ~/.bashrc

cd <project-dir>
git clone https://github.com/ChiaLungChuang/mb-scrnaseq.git
```

The cluster's own `nextflow` module is older (23.10) than this pipeline
needs (25.10 or newer), which is why it is installed with micromamba.

## 3. Smoke test on the cluster

```bash
cd <project-dir>/mb-scrnaseq
tmux new -s mbtest
micromamba activate nf
nf-test test --profile +stjude
```

This runs both synthetic test suites through LSF and Singularity with the
same checks as CI, including the exact-count accuracy test. It takes
minutes of compute but may sit in the queue first. It proves the
scheduler, the containers and the paths work, and give the same numbers,
before any real data is involved. (`nextflow run . -profile test,stjude`
runs the pipeline alone, without the checks.) Detach with
`Ctrl-b d`; reattach with `tmux attach -t mbtest`. Watch jobs with `bjobs`.

## 4. Fetch the inputs

```bash
module load singularity/4.1.1
D=<project-dir>/gse155446
scripts/fetch_inputs.sh geo        $D/geo
scripts/fetch_inputs.sh whitelists $D/whitelists
scripts/fetch_inputs.sh reference  $D/reference GRCh38-3.0.0    # see step 6
```

`geo` also writes `$D/geo/published_labels.tsv`: every sample label in the
published matrix with its number of cells, taken from the column names
(the file's first line) only. No count values are read.

## 5. Choose and download the samples **(you)**

Open the SRA Run Selector for project **PRJNA649773** and download the
metadata table. Pick 8 human tumors that span the subgroups (the cohort has
1 WNT, 9 SHH, 7 Group 3 and 11 Group 4). Put their SRR accessions in
`runs.txt`, one per line. Check sizes first (BAMs are 6-25 GB each), then
download:

```bash
DRY_RUN=1 scripts/fetch_inputs.sh bams $D/bams runs.txt   # lists file names, sizes, MD5s
scripts/fetch_inputs.sh bams $D/bams runs.txt
```

BAMs come from ENA, which keeps the files as submitted. Each one is checked
against ENA's size and MD5, with `samtools quickcheck` (truncation) and for
the 10x header lines bamtofastq needs. `$D/bams/bams_summary.tsv` lists, per
BAM, the read-group sample name, the Cell Ranger reference and the STAR
version from its header.

To test the whole chain cheaply first, start with SRR12354589 (6.4 GB).

## 6. Match the reference and Cell Ranger version

Read the `reference` and `star_version` columns of `bams_summary.tsv`. The
two BAMs checked so far say `refdata-cellranger-GRCh38-3.0.0` and
`STAR_2.5.1b`, which is Cell Ranger 3.x. Then use reference `GRCh38-3.0.0`
(step 4) and `cellranger_compat: cr3` (step 9). If any BAM says otherwise,
stop and work out why before going on. A different gene annotation from the
authors' lowers gene matching and is a confound, not a pipeline result.

## 7. Fill the samplesheet, looking only at names

Copy `assets/samplesheet_gse155446.csv` and add one row per BAM.

- `published_label`: the label in `$D/geo/published_labels.tsv` that
  matches the BAM's read-group sample or file name in `bams_summary.tsv`
  (e.g. file `1_943_possorted_genome_bam.bam` → label `943`). Do **not**
  use the SRA run title: SRR12354571 is titled "934-repeat" but its BAM is
  tumor 943, and both 934 and 943 are in the published matrix. A wrong
  label fails C1 (barcode recovery), but it costs a run.
- `chemistry`: v2 or v3. Expect v3 for all but label `1355` (SRR12354589),
  whose barcodes all lie in the v2 whitelist and whose UMIs are 10 nt. If
  unsure, enter your best guess: the pipeline stops with the right answer
  from the read lengths.

Do not open the published counts themselves before step 8.

## 8. Lock the acceptance criteria **(you)**

Read `preregistration/acceptance_criteria.yaml`. Change any threshold you
disagree with now, with the reason in the commit message. Then set
`status: locked`, commit and **push**. Optionally tag the commit and make a
GitHub release from it (e.g. `prereg-phase1`): that gives the lock a
timestamp you did not set yourself.

On real data the pipeline refuses the concordance step unless the criteria
are locked, committed, unmodified and pushed, and the pipeline code has no
uncommitted changes (untracked files such as your params file and
samplesheet are fine). Every run records the criteria file's hash and
commit, and the verdict says whether the run was pre-registered.

## 9. Real run

```bash
cp assets/params_gse155446.yaml my_params.yaml    # fill in the paths
tmux new -s mbreal
micromamba activate nf
nextflow run . -profile stjude -params-file my_params.yaml -w <project-dir>/work -resume
```

The first run also builds the STAR index (about 32 GB of RAM and roughly an
hour). It is copied to `<outdir>/reference/star_index` (about 30 GB); for
later runs, set `star_index` to that path and remove
`genome_fasta`/`genome_gtf`.

## 10. Read the result

- `concordance/concordance_summary.verdict.json`: overall pass/fail, and
  `preregistered` (must be `true` for a reportable result).
- `concordance/concordance_summary.tsv`: each sample × criterion, with the
  observed value and the rule.
- `pipeline_info/run_provenance.json`: commit, code hash, parameters,
  container engine, and the criteria file's hash and lock state.
- `qc/*.qc_thresholds.json`: per-sample QC thresholds actually applied.
- `starsolo/<sample>/`: the count matrices and STARsolo's `Summary.csv`.
- `pipeline_info/report.html`, `timeline.html` and `trace.tsv`: resources
  used per step, which is what you quote when someone asks what the run
  cost. Check `peak_rss` for `CONCORDANCE`; it has not been measured on
  full-size data.

If a criterion fails, do not edit the threshold. Find out why: wrong
reference version (low `gene_match_fraction`), wrong sample-to-label
mapping (low `barcode_recovery`), different cell calling (compare
`published_cells_called_by_starsolo`). Then report the failure and the
explanation together.

## When something breaks

| Symptom | Likely cause |
|---|---|
| `bsub: command not found` | Not on a login node, or LSF not in PATH there |
| Singularity pull times out | First pull of a large image; re-run with `-resume` |
| `ghcr.io ... unauthorized` | The R image package is still private (step 3 in section 0) |
| `CHECK_READS` stops: chemistry needs 26/28 | Samplesheet chemistry is wrong for that sample; fix and `-resume` |
| `label '...' not found` | `published_label` does not match `published_labels.tsv` |
| `Refusing to run the concordance check` | The message lists what is missing: lock, commit, push, or uncommitted code |
| Job killed, exit 130-145 | Out of memory or time; retried automatically with more, twice. Confirm with `bhist -l <job id>` |
| One sample fails, others keep going | Expected: no new jobs start, running ones finish, then the run stops. Fix the cause and `-resume` |
| A step re-runs after `-resume` | You changed code or the criteria; the steps that use it re-run by design. R changes re-run only the R steps; a change to `bin/check_reads.sh` also re-runs STARsolo |
