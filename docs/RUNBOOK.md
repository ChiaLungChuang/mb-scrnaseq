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
nextflow run . -profile test,stjude
```

This runs the synthetic test through LSF and Singularity. It takes minutes
of compute but may sit in the queue first. It proves the scheduler, the
containers and the paths work before any real data is involved. Detach with
`Ctrl-b d`; reattach with `tmux attach -t mbtest`. Watch jobs with `bjobs`.

## 4. Fetch the inputs

```bash
module load singularity/4.1.1
D=<project-dir>/gse155446
scripts/fetch_inputs.sh geo        $D/geo
scripts/fetch_inputs.sh whitelists $D/whitelists
scripts/fetch_inputs.sh reference  $D/reference <version>   # see step 6
```

## 5. Choose the samples **(you)**

Open the SRA Run Selector for project **PRJNA649773** and download the
metadata table. Pick 8 human tumors that span the subgroups (the cohort has
1 WNT, 9 SHH, 7 Group 3 and 11 Group 4). Put their SRR accessions in
`runs.txt`, one per line, then:

```bash
scripts/fetch_inputs.sh bams $D/bams runs.txt
```

Each BAM is checked for the 10x header lines before it is accepted.

## 6. Match the reference to the paper **(you)**

On the GSE155446 GEO page, read the samples' *data processing* field for the
Cell Ranger version and reference used. Use `GRCh38-3.0.0` if it names
Cell Ranger 3.x with the 3.0.0 reference, `GRCh38-2020-A` if it names
2020-A. A different gene annotation than the authors' lowers gene matching
and is a confound, not a pipeline result.

## 7. Fill the samplesheet, looking only at names

Copy `assets/samplesheet_gse155446.csv`, add one row per BAM. For
`published_label`, look at the published matrix's **column names only**:

```bash
zcat $D/geo/GSE155446_human_raw_counts.csv.gz | head -n 1 | tr ',' '\n' | sed -n '2,6p'
```

Do not open the counts themselves before step 8. `chemistry` is v2 or v3;
if unsure, enter your best guess and the pipeline will stop with the right
answer from the read lengths.

## 8. Lock the acceptance criteria **(you)**

Read `preregistration/acceptance_criteria.yaml`. Change any threshold you
disagree with now, with the reason in the commit message. Then set
`status: locked`, commit and push. The pipeline will not run the
concordance step on real data until this is done, and every run records the
file's hash and commit.

## 9. Real run

```bash
cp assets/params_gse155446.yaml my_params.yaml    # fill in the paths
tmux new -s mbreal
micromamba activate nf
nextflow run . -profile stjude -params-file my_params.yaml -w <project-dir>/work -resume
```

The first run also builds the STAR index (about 32 GB of RAM and roughly an
hour). It is published to `<outdir>/reference/star_index`; for later runs,
set `star_index` to that path and remove `genome_fasta`/`genome_gtf`.

## 10. Read the result

- `concordance/concordance_summary.verdict.json`: overall pass/fail.
- `concordance/concordance_summary.tsv`: each sample × criterion, with the
  observed value and the rule.
- `pipeline_info/run_provenance.json`: commit, parameters, container engine,
  and the criteria file's hash and lock state.
- `qc/*.qc_thresholds.json`: per-sample QC thresholds actually applied.
- `pipeline_info/report.html` and `timeline.html`: resources used per step,
  which is what you quote when someone asks what the run cost.

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
| `label '...' not found` | `published_label` does not match the published cell IDs |
| Job killed, exit 137/140 | Out of memory or time; it retries with more automatically, twice |
