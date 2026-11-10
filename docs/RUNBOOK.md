# Runbook: first real run on St. Jude compute

This is the order to do things in, from no account to a concordance verdict.
Steps marked **(you)** need your accounts, judgement or St. Jude-specific
details that are not public; everything else is a command to run.

## Two routes

Phase 1 runs on St. Jude compute (decided 10 Nov 2026), by one of two
routes. Nothing has run on either yet.

| Route | Where | Sections | Waiting on |
|---|---|---|---|
| 1. The Research HPC cluster | LSF + Singularity, `hpc.stjude.org` (`-profile stjude`) | 0-10 | **(you)** the HPCF account request, and HPCF's answers to the questions in section 1 (meeting upcoming) |
| 2. A department server | Developmental Neurobiology machines offered by the NeuroImage Analysis Lab, if they are standalone servers without a scheduler (`-profile server`) | "Route 2" at the end; it reuses 4-8 and 10 | **(you)** the lab's reply on what the machines are. A machine shared through a scheduler needs that scheduler's profile instead (not written yet) |

Both routes share the inputs and samplesheet (sections 4-7), the locked
criteria and the pre-registration gate (section 8), the disk estimate
(section 9) and the results report (section 10). Google Cloud and Colab
are not being pursued for now.

## How the cluster fits together (two minutes)

HPCF's rules (HPCF wiki, "Introduction to the HPCF cluster") shape every
step below.

- You `ssh` into **`hpc.stjude.org`** with your AD (St. Jude) credentials.
  That puts you on a **login node**, which is only for changing directories
  and submitting jobs: even `cp` or `cat` belongs on a compute node. So the
  first thing after logging in is an **interactive session**, a shell on a
  compute node:
  ```bash
  hpcf_interactive                                  # or:
  bsub -P interactive -q interactive -Is bash
  ```
  The interactive queue has a **16-hour run limit**, and `hpcf_interactive`
  reports a **48-hour accumulated CPU-time limit**. Every step in this
  runbook that does more than `cd` or `bsub` runs in such a session or in a
  batch job.
- Nextflow runs on a compute node as a lightweight **head process**: in an
  interactive session, or as a batch job of its own (section 9 gives both).
  For every step it writes a small job script and submits it to **LSF**
  (`bsub`), the scheduler that places jobs on compute nodes with the CPUs,
  memory and time each step asked for (see `conf/base.config`). So the
  compute node the head process runs on must be allowed to submit jobs;
  `scripts/hpcf_preflight.sh --submit-test` checks that (section 3).
- The `stjude` profile picks the queue for you: `short` (30-minute hard
  limit), `standard` (the default), or `large_mem` (above 1 TB, which no
  step of this pipeline requests). HPCF also has `heavy_io` and
  `interactive`. Per user, at most 50,000 jobs pending and 5,000 running;
  the profile keeps at most 100 in the queue.
- Memory is requested **per core**: LSF multiplies `rusage[mem]` by the
  job's slots. The profile divides each step's memory by its CPUs, so an
  8-CPU, 48 GB step asks for 6144 MB per slot (`conf/stjude.config`,
  "Memory").
- Each job runs inside a **Singularity** container. The head process pulls
  each image once into `NXF_SINGULARITY_CACHEDIR`, from the node it runs on,
  so that node needs `singularity` on its `PATH` (`module load
  singularity/4.1.1`) and access to ghcr.io and quay.io. That is how the
  exact STAR, bamtofastq and R versions follow the pipeline onto the
  cluster. Docker is not available there.
- The head process has to stay alive for the whole run (hours). If it
  dies, re-run the same command with `-resume` and finished steps are
  skipped.
- Intermediate files live in `work/`. Put it on project or scratch storage,
  not your home directory: the FASTQs regenerated from 8 BAMs run to 1 to
  1.7 times the size of the BAMs. The whole 8-tumor run needs about
  350-445 GB at its peak (section 9, "Disk space").
- Downloads (setup, inputs, images) need outbound internet. Whether compute
  nodes have it is **pending HPCF's answer** (the preflight shows what one
  compute node can reach). If they do not, those steps run wherever HPCF
  says downloads should run, never on a login node.

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

## 1. Account, storage and questions for HPCF

4. **(you)** Request an HPCF account through IT, if you do not have one.
   The login host is `hpc.stjude.org`, with your AD credentials.
5. **(you)** Find out from HPCF documentation or your lab which project or
   scratch path you can use for large data, with about 0.5 TB free (section
   9, "Disk space"). Below it is `<project-dir>`.
6. **(you)** Ask HPCF (or find in their documentation) the answers this
   runbook marks **pending HPCF's answer**:
   1. Do compute nodes have outbound HTTPS (GitHub, conda-forge and
      bioconda, ghcr.io, quay.io, ENA, GEO, 10x Genomics, labshare.cshl.edu)?
      If not, where should downloads and Singularity image pulls run (a
      data-transfer node?), given that login nodes are not for them?
   2. May a job (an interactive session or a batch job) submit other jobs?
      Nextflow's head process does this for every step.
   3. Is `tmux` on a login node acceptable for keeping an interactive
      session alive across disconnects? Do repeat logins to
      `hpc.stjude.org` reach the same login node?
   4. For a long, low-CPU batch job running Nextflow's head process (1 core,
      8 GB, up to a few days): which queue, and what run limit?
   5. How is `-M` enforced: per process (LSF's default), per job
      (`LSB_JOB_MEMLIMIT=Y`) or by cgroups (`LSB_RESOURCE_ENFORCE`)? STAR is
      one multi-threaded process that can use the whole 48 GB of an 8-slot
      job; a per-process limit of 6144 MB would kill it
      (`conf/stjude.config`, "Memory"). The preflight prints the lsf.conf
      settings that answer this.
   6. Should the BAM download (about 119 GB, mostly I/O) run in `heavy_io`?

## 2. One-time setup, in an interactive session

```bash
ssh <AD user>@hpc.stjude.org     # login node: cd and bsub only
hpcf_interactive                 # or: bsub -P interactive -q interactive -Is bash
# now on a compute node

# micromamba in your home directory (small; no admin rights needed)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# open a new shell (exec bash), then:
micromamba create -y -n nf -c conda-forge -c bioconda nextflow=26.10.4 nf-test=0.9.5
micromamba activate nf
nextflow -version        # should print 26.10.4

# Singularity image cache (there is no shared one on HPCF), a log folder for
# batch jobs, and a fixed Java heap for Nextflow's head process (its batch
# job in section 9 reserves 8 GB)
mkdir -p <project-dir>/singularity_cache <project-dir>/logs
echo 'export NXF_SINGULARITY_CACHEDIR=<project-dir>/singularity_cache' >> ~/.bashrc
echo "export NXF_OPTS='-Xms1g -Xmx4g'" >> ~/.bashrc
source ~/.bashrc

cd <project-dir>
git clone https://github.com/ChiaLungChuang/mb-scrnaseq.git
```

The cluster's own `nextflow` module is older (23.10) than this pipeline
needs (25.10 or newer), which is why it is installed with micromamba.

This step downloads (micromamba, conda packages, the repository). If the
compute node has no outbound internet, run it where HPCF says downloads
should run (**pending HPCF's answer**, section 1, question 1). This assumes
home and project storage are shared with the compute nodes, as on most
clusters.

## 3. Smoke test on the cluster

First the preflight, in an interactive session, so that what it reports is
what a compute node sees. It checks:

- that it runs inside an LSF job;
- LSF and the queues `conf/stjude.config` sends jobs to, and lsf.conf's
  memory settings (printed for question 5, not judged);
- that the Singularity module loads;
- the git, Java, Nextflow and nf-test versions;
- the Singularity cache, and free space and quota where `work/` will live;
- that every container image the pipeline uses can be reached anonymously;
- that the download hosts of `scripts/fetch_inputs.sh` answer (one ENA
  query for the pilot run, then HEAD requests: nothing is downloaded).

With `--submit-test` it also submits one 1-minute job (`true`, queue
`short`, no output files) from this node and waits for it to finish, up to
15 minutes: that is question 2, tested. Without the option it writes no
file, pulls no image and submits no job. Paste the whole report back; fix
every FAIL before the smoke test (WARN lines are worth reading but do not
block).

```bash
ssh <AD user>@hpc.stjude.org
hpcf_interactive
cd <project-dir>/mb-scrnaseq
micromamba activate nf
scripts/hpcf_preflight.sh --submit-test <project-dir>   # the directory that will hold work/
```

If the registry and download checks FAIL here, compute nodes have no
outbound internet: the image pulls and the downloads of sections 4 and 5
then wait for **HPCF's answer** to question 1. If the submit test FAILs,
the head process cannot run on a compute node; that is question 2.

Then the smoke test. It is short, so an interactive session is enough.
To survive a dropped connection, start `tmux` on the login node first and
the session inside it (**pending HPCF's answer**, question 3; without
`tmux`, keep the connection open, or submit it as a batch job like the
real run in section 9, with `nf-test test --profile +stjude` as the
command):

```bash
ssh <AD user>@hpc.stjude.org
tmux new -s mbtest          # on the login node; only holds the terminal
hpcf_interactive            # the rest runs on a compute node
cd <project-dir>/mb-scrnaseq
micromamba activate nf
module load singularity/4.1.1    # the head process pulls the images
nf-test test --profile +stjude
```

This runs both synthetic test suites through LSF and Singularity with the
same checks as CI, including the exact-count accuracy test. It takes
minutes of compute but may sit in the queue first; that should fit well
inside the session's 16-hour limit. It proves the scheduler, the
containers and the paths work, and give the same numbers, before any real
data is involved. (`nextflow run . -profile test,stjude` runs the pipeline
alone, without the checks.) Detach with `Ctrl-b d`; reattach with `tmux attach -t mbtest`
from the same login node (`hostname` before detaching tells you which).
Watch jobs with `bjobs`.

## 4. Fetch the inputs

In an interactive session, if the preflight showed the download hosts
reachable from it; otherwise where HPCF says downloads should run
(**pending HPCF's answer**, question 1). Not on a login node.

```bash
cd <project-dir>/mb-scrnaseq
module load singularity/4.1.1
D=<project-dir>/gse155446
scripts/fetch_inputs.sh geo        $D/geo
scripts/fetch_inputs.sh whitelists $D/whitelists
scripts/fetch_inputs.sh reference  $D/reference GRCh38-3.0.0    # see step 6
```

`geo` also writes `$D/geo/published_labels.tsv`: every sample label in the
published matrix with its number of cells, taken from the column names
(the file's first line) only. No count values are read.

## 5. Download the 8 samples **(you)**

The 8 tumors were chosen on 2 Nov 2026 (proposal:
`docs/proposals/phase1_samples.md`) and are recorded, in run order, in
`docs/results/phase1/samples.tsv`: label, SRR, chemistry, subgroup, BAM
file name, MD5 and size, taken from the barcode-verified run table
`docs/proposals/gse155446_runs.tsv` (all 30 human runs). Together about
119 GB. Check sizes first, then download; the pilot (1355, SRR12354589,
6.4 GB) can go alone first:

```bash
awk -F'\t' 'NR > 1 { print $2 }' docs/results/phase1/samples.tsv > runs.txt
DRY_RUN=1 scripts/fetch_inputs.sh bams $D/bams runs.txt   # lists file names, sizes, MD5s
scripts/fetch_inputs.sh bams $D/bams runs.txt
```

Where to run the download (**pending HPCF's answer**, questions 1 and 6):

- **In the interactive session** (under `tmux`, as in section 3), if
  compute nodes reach ENA. 119 GB takes about 3.3 hours at 10 MB/s, inside
  the 16-hour limit; if the session ends, re-run the same command:
  finished files are skipped and a partial one resumes.
- **As a batch job**, which needs no open terminal. From the interactive
  session, after `module load singularity/4.1.1` (LSF copies this shell's
  environment into the job):
  ```bash
  bsub -q standard -n 1 -M 4096MB -R "rusage[mem=4096MB]" -W 24:00 -J mb_fetch \
    -o <project-dir>/logs/mb_fetch.%J.log scripts/fetch_inputs.sh bams $D/bams runs.txt
  ```
  `heavy_io` may be the right queue instead of `standard` (question 6).
- **Wherever HPCF says downloads go**, if compute nodes have no internet.

BAMs come from ENA, which keeps the files as submitted. Each one is checked
against ENA's size and MD5, with `samtools quickcheck` (truncation) and for
the 10x header lines bamtofastq needs. `$D/bams/bams_summary.tsv` lists, per
BAM, the read-group sample name, the Cell Ranger reference and the STAR
version from its header.

The file names and MD5s in `bams_summary.tsv` must match those in
`samples.tsv`.

## 6. Match the reference and Cell Ranger version

Read the `reference` and `star_version` columns of `bams_summary.tsv`. All
30 human BAM headers, read from their first 4 MB for the proposal, say
`refdata-cellranger-GRCh38-3.0.0` and `STAR_2.5.1b`, which is Cell Ranger
3.x. Then use reference `GRCh38-3.0.0` (step 4) and `cellranger_compat:
cr3` (step 9). If any downloaded BAM says otherwise, stop and work out why
before going on. A different gene annotation from the authors' lowers gene
matching and is a confound, not a pipeline result.

Also read `star/genomeParameters.txt` in the downloaded reference. The
pipeline builds its index with `star_sjdb_overhang: 100` (the cDNA reads
are 151 nt). If Cell Ranger's index says a different `sjdbOverhang`, stop:
which to use is a decision for Chia-Lung, because it changes alignments
slightly.

## 7. Fill the samplesheet, looking only at names

Copy `assets/samplesheet_gse155446.csv` and add one row per BAM, in the run
order of `docs/results/phase1/samples.tsv`.

- `sample`: the `sample` column of `samples.tsv` (the GSM accession).
- `published_label`: the `published_label` column of `samples.tsv`, which
  comes from `docs/proposals/gse155446_runs.tsv`, where every run was mapped
  to its label by barcode overlap. Check that the label is in
  `$D/geo/published_labels.tsv`. Do **not** derive it from the SRA/GEO
  title or from the read-group name alone: SRR12354571 (GEO "934-repeat")
  is published label 943, most likely a second library of tumor 934 (both
  934 and 943 are in the published matrix), and the read groups
  `Foreman_1` to `Foreman_4` carry no label at all (tumors 801, 877, 996,
  1070). A wrong label fails C1 (barcode recovery), but it costs a run.
- `chemistry`: the `chemistry` column of `samples.tsv`: v3 for all but
  label `1355` (SRR12354589), whose UMIs are 10 nt. The pipeline checks it
  against the read lengths and stops if it is wrong.

Do not open the published counts themselves before step 8.

## 8. The locked acceptance criteria

Done on 10 Nov 2026: `preregistration/acceptance_criteria.yaml` is
`status: locked`, thresholds unchanged from the proposal ([#22](https://github.com/ChiaLungChuang/mb-scrnaseq/pull/22)).
On a fresh clone the pre-registration gate passes (a stub run with
`--allow_unlocked_criteria false`). **(you)** Make the GitHub release
`prereg-phase1` from that commit: it gives the lock a timestamp you did
not set yourself. Pull the clone you run from to get the locked file
(GitHub access from that machine: **pending HPCF's answer**, question 1,
on the cluster).

On real data the pipeline refuses the concordance step unless the criteria
are locked, committed, unmodified and pushed, and the pipeline code has no
uncommitted changes (untracked files such as your params file and
samplesheet are fine). Every run records the criteria file's hash and
commit, and the verdict says whether the run was pre-registered.

## 9. Real run

The head process runs on a compute node, in one of two ways. Both start
from an interactive session:

```bash
ssh <AD user>@hpc.stjude.org
tmux new -s mbreal               # path A only (question 3)
hpcf_interactive
cd <project-dir>/mb-scrnaseq
cp assets/params_gse155446.yaml my_params.yaml    # once; fill in the paths
micromamba activate nf
module load singularity/4.1.1    # the head process pulls the images
```

**Path A: in the interactive session, under `tmux`.** Simple, and you see
progress live, but the session's limits (16 hours run, 48 hours CPU time)
end the run if it is not finished; `-resume` then continues from the last
finished steps, and steps that were running start again. Good for the
pilot sample alone.

```bash
nextflow run . -profile stjude -params-file my_params.yaml -w <project-dir>/work -resume
```

**Path B: as a batch job.** No terminal to keep open and no 16-hour limit;
recommended for all 8 samples. Submit it from the same interactive
session, which LSF copies into the job (`micromamba activate nf`, the
Singularity module, `NXF_SINGULARITY_CACHEDIR`, `NXF_OPTS`), then end the
session:

```bash
bsub -q standard -n 1 -M 8192MB -R "rusage[mem=8192MB]" -W 72:00 -J mb_head \
  -o <project-dir>/logs/mb_head.%J.log \
  nextflow run . -profile stjude -params-file my_params.yaml -w <project-dir>/work -resume
exit
```

The queue and the 72-hour limit are **pending HPCF's answer** (question
4). The job needs one core and 8 GB (Nextflow's Java heap is capped at
4 GB by `NXF_OPTS`, section 2). Follow it from an interactive session with
`bjobs -w` and `tail -f <project-dir>/logs/mb_head.<job id>.log`. Both
paths need a compute node that may submit jobs (question 2, the
preflight's `--submit-test`).

The first run also builds the STAR index (about 32 GB of RAM and roughly an
hour). It is copied to `<outdir>/reference/star_index` (about 30 GB); for
later runs, set `star_index` to that path and remove
`genome_fasta`/`genome_gtf`.

### Disk space

Nextflow keeps every task's outputs in `work/` until you delete them, so
at the end of the run everything below is on disk at once. "Measured" means
measured in this repository on synthetic data; the rest are estimates.

| What | Where | Size | Basis |
|---|---|---|---|
| 8 BAMs | `$D/bams` | 119 GB | ENA's sizes, `docs/results/phase1/samples.tsv` |
| FASTQs from bamtofastq | `work/` | 119-202 GB | 1.0-1.7 times the BAMs (below) |
| STAR index | `work/` | ~30 GB | estimate for a human STAR index |
| its published copy | `<outdir>/reference/` | ~30 GB | `publishDir` copies it |
| Cell Ranger reference, tarball and extracted | `$D/reference` | ~11 + ~30 GB | estimate; mostly Cell Ranger's own STAR index, which this pipeline does not use |
| STARsolo matrices, in `work/` and `<outdir>` | both | 7-18 GB | measured: 0.4-1.2 GB per sample at real dimensions (`tests/scale/README.md`), twice |
| Singularity images | the cache | ~1 GB | the R image is 472 MB (its registry manifest); the others are small |
| GEO files, published matrix, QC and concordance outputs | | ~1-2 GB | measured at real dimensions: 0.1-0.3 GB each for the CSV.gz and the RDS |
| **Peak, end of the run** | | **~350-445 GB** | |

The FASTQ size is the uncertain part. bamtofastq on the two synthetic test
BAMs wrote FASTQ.gz files 1.68 times the BAM size; with realistic base
qualities and Cell Ranger's alignment tags added to the same reads, 1.01
times. Real BAMs differ from both: they are sorted by position, which
makes the BAM compress better (a higher ratio), and they carry more tags
(a lower one). The pilot gives the real ratio: compare its BAM with
`du -sh` of its `BAMTOFASTQ` task directory in `work/`, and redo the FASTQ
line before the other 7.

**Is 0.5 TB enough?** For one clean run, yes. At the top of the range only
about 55 GB is left, and a failed task leaves its partial output in
`work/` (up to ~53 GB for a retried bamtofastq of the largest BAM, 996).
Before the 8-sample run:

- After section 6, delete the reference tarball and the extracted
  reference's `star/` folder (Cell Ranger's own index): about 36 GB, by
  estimate. The pipeline reads only `fasta/` and `genes/`.
- When a task was retried and the retry succeeded, delete the failed
  attempt's directory (the path is in the error message); `-resume` uses
  only the successful one.

What does not lower the peak:

- `cleanup = true` deletes `work/` only when the whole run has succeeded,
  and `-resume` no longer works after it (Nextflow 26.10's documentation).
  The pilot-then-all-samples plan relies on `-resume`.
- Per-task scratch (`process.scratch`) runs each task on the node's local
  disk, but the declared outputs (the FASTQs, the STARsolo matrices) are
  still copied back to `work/` when the task ends. Only temporary files
  move, and STAR deletes its own. It would also need local disks large
  enough for the largest bamtofastq output (~53 GB; **pending HPCF's
  answer**).
- (`publishDir` with `mode: 'link'` would save the index's second copy,
  ~30 GB, when `<outdir>` and `work/` are on the same file system, but
  fails the run when they are not. Not changed.)

Once the verdict is final and its small files are committed (plan M3), the
FASTQs in `work/` can be deleted; the published index and matrices stay.

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
  cost. Compare `peak_rss` and `realtime` of the R steps with the
  rehearsal on synthetic data at real dimensions (`tests/scale/README.md`:
  at most 4.8 GB for `PREP_PUBLISHED` and 5.8 GB for `CONCORDANCE`).

Then write the results report from the output directory, in an
interactive session, with the pinned R image (it has every package the
script needs):

```bash
module load singularity/4.1.1
IMG=docker://$(grep -oP "^\s*r_container\s*=\s*'\K[^']+" nextflow.config)
singularity exec -B <project-dir> "$IMG" Rscript bin/phase1_report.R --lib lib/R \
  --outdir <project-dir>/results --out docs/results/phase1/README.md
```

`docs/results/phase1/README.md` then holds the verdict per sample and
criterion, whether the run was pre-registered (a run that was not is
marked as not reportable), the commit and the criteria and code hashes,
QC flags, and compute per process (CPU hours used and reserved, run time,
peak memory) from the trace. It refuses an output directory whose files
disagree (for example, a verdict from one run and provenance from
another). Copy the small result files it summarises next to it (plan M3
step 5).

If a criterion fails, do not edit the threshold. Find out why: wrong
reference version (low `gene_match_fraction`), wrong sample-to-label
mapping (low `barcode_recovery`), different cell calling (compare
`published_cells_called_by_starsolo`). Then report the failure and the
explanation together.

## Route 2: a standalone server (no scheduler)

For one Linux server or workstation with no job scheduler, such as a
departmental machine. Nextflow then runs every step on that machine with
its local executor (`-profile server`, `conf/server.config`), in
Singularity, Apptainer or Docker. Steps 4 to 8 and 10 above are the same;
what changes is setup, the preflight, the smoke test and how the run is
started.

What the machine needs:

- **Memory: more than 48 GB.** The STAR index build and STARsolo ask for
  48 GB (`conf/base.config`). The pipeline refuses to start when
  `--server_memory` is below 48 GB or above what the machine has, so a
  machine with exactly 48 GB (about 47 GB usable) is not enough; 64 GB is.
- **Cores: 8 or more** for the STAR steps to run at their requested
  speed. With fewer, every task is capped to `--server_cpus` and runs
  slower.
- **Disk**: the same as the cluster run (section 9, "Disk space": about
  350-445 GB at the peak for 8 tumors). The pipeline refuses to start when
  the work directory's file system cannot hold one tumor: the FASTQs of
  the largest BAM (1.7 times it), a STAR index built there and its copy,
  and `--server_disk_margin` (20 GB).
- **A container runtime**: Singularity or Apptainer (no root needed to
  run), or Docker with your user in the `docker` group.

`--server_cpus` and `--server_memory` are the share of the machine the run
may use. On a shared machine, agree it with the people who run it, and
leave memory for them and the system. Every task is capped at that share
(`process.resourceLimits`), and the local executor starts a task only when
its cores and memory fit in the share.

### Setup

As in section 2, without `hpcf_interactive`: micromamba with Nextflow
26.10.4 and nf-test 0.9.5, the clone, a log folder, and an image cache:

```bash
mkdir -p <project-dir>/singularity_cache <project-dir>/logs
echo 'export NXF_SINGULARITY_CACHEDIR=<project-dir>/singularity_cache' >> ~/.bashrc   # Singularity
echo 'export NXF_APPTAINER_CACHEDIR=<project-dir>/singularity_cache' >> ~/.bashrc     # or Apptainer
echo "export NXF_OPTS='-Xms1g -Xmx4g'" >> ~/.bashrc
```

Docker needs no cache variable.

### Preflight and smoke test

```bash
cd <project-dir>/mb-scrnaseq
micromamba activate nf
scripts/hpcf_preflight.sh --local <project-dir>     # the directory that will hold work/
```

`--local` skips the LSF and module checks. It reports the cores, FAILs
below 48 GB of memory and without Singularity, Apptainer or a usable
Docker, and warns if a scheduler is on the `PATH`, in which case its
profile is the better route. The tools, cache, storage, registries and
download hosts are checked as on the cluster. Paste the report back.

Then the synthetic test suites on this machine (2 cores and 6 GB, set by
the test profile; replace `singularity` with `apptainer` or `docker`):

```bash
tmux new -s mbtest
nf-test test --profile +server,singularity
```

### Inputs and the lock

Sections 4 to 8 unchanged: `scripts/fetch_inputs.sh` for GEO, the
whitelists, the reference and the BAMs; the reference check; the
samplesheet; the locked criteria (the pre-registration gate works the
same way here). Downloads need outbound HTTPS from this machine; the
preflight's network checks show whether it has it.

### Real run

```bash
cd <project-dir>/mb-scrnaseq
cp assets/params_gse155446.yaml my_params.yaml    # once; fill in the paths
tmux new -s mbreal
micromamba activate nf
nextflow run . -profile server,singularity -params-file my_params.yaml \
  --server_cpus 16 --server_memory '96.GB' -w <project-dir>/work -resume
```

(16 cores and 96 GB are an example: use your share.) At launch the
pipeline writes `pipeline_info/server_check.json` with the share, the
machine's cores and memory, the free space and one tumor's needs, or
refuses to start and says which check failed. Detach with `Ctrl-b d`; if
the run stops, re-run the same command and finished steps are skipped.

### The report

Section 10 unchanged, with the container runtime you use:

```bash
IMG=$(grep -oP "^\s*r_container\s*=\s*'\K[^']+" nextflow.config)
singularity exec -B <project-dir> "docker://$IMG" Rscript bin/phase1_report.R --lib lib/R \
  --outdir <project-dir>/results --out docs/results/phase1/README.md
# or: apptainer exec ... (same arguments), or
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -v <project-dir>:<project-dir> -w "$PWD" \
  "$IMG" Rscript bin/phase1_report.R --lib lib/R \
  --outdir <project-dir>/results --out docs/results/phase1/README.md
```

## When something breaks

| Symptom | Likely cause |
|---|---|
| `Refusing to start on this server: ...` | Route 2: the share (`--server_cpus`, `--server_memory`) is missing, below 48 GB or more than the machine has, or the work directory's disk cannot hold one tumor; the message names the check |
| `bsub: command not found` | LSF is not set up in this shell; on a compute node, start from `hpcf_interactive` or a batch job, which carry LSF's environment |
| The preflight's submit test fails, or Nextflow's jobs are refused | This compute node may not submit jobs: ask HPCF (section 1, question 2) |
| The run stops after 16 hours in an interactive session | The interactive queue's run limit: `-resume` in a new session, or use the batch job (section 9, path B) |
| Singularity pull fails with `singularity: command not found` | The head process needs `module load singularity/4.1.1` in its shell |
| Singularity pull times out | First pull of a large image; re-run with `-resume` |
| `ghcr.io ... unauthorized` | The R image package is still private (step 3 in section 0) |
| `CHECK_READS` stops: chemistry needs 26/28 | Samplesheet chemistry is wrong for that sample; fix and `-resume` |
| `label '...' not found` | `published_label` does not match `published_labels.tsv` |
| `Refusing to run the concordance check` | The message lists what is missing: lock, commit, push, or uncommitted code |
| Job killed, exit 130-145 | Out of memory or time; retried automatically with more, twice. Confirm with `bhist -l <job id>` |
| One sample fails, others keep going | Expected: no new jobs start, running ones finish, then the run stops. Fix the cause and `-resume` |
| A step re-runs after `-resume` | You changed code or the criteria; the steps that use it re-run by design. R changes re-run only the R steps; a change to `bin/check_reads.sh` also re-runs STARsolo |
