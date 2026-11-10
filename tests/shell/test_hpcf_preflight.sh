#!/usr/bin/env bash
# Tests for scripts/hpcf_preflight.sh. LSF, Java, Nextflow, nf-test, the
# module system, curl, df and sleep are stand-ins; the script runs with `env -i`
# and a PATH holding only those stand-ins and the core tools it needs, so a
# real java or curl on the test machine cannot leak in.
# shellcheck disable=SC2016  # stand-in bodies are written out unexpanded; they expand when run
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
PF="$REPO/scripts/hpcf_preflight.sh"

SBX=$TMP/bin            # the only PATH the script sees
SING=$TMP/singbin       # added to PATH by the stand-in `module load`
CALLS=$TMP/calls.log    # every stand-in logs its command line here
WORK=$TMP/project/work
LSFENV=$TMP/lsf         # LSF_ENVDIR: holds a stand-in lsf.conf
mkdir -p "$SBX" "$SING" "$TMP/home" "$TMP/cache" "$WORK" "$LSFENV"
for t in bash sed awk grep sort head tail cut date dirname git id hostname cat tr; do
  ln -s "$(command -v "$t")" "$SBX/$t"
done

stub() { # stub <dir> <name> <body>: a logging stand-in command
  printf '#!/bin/bash\necho "%s $*" >> "$CALLS"\n%s\n' "$2" "$3" > "$1/$2"
  chmod +x "$1/$2"
}
BSUB_OK='[ "${STUB_BSUB:-ok}" = refuse ] && { echo "User permission denied. Job not submitted." >&2; exit 255; }
echo "Job <4321> is submitted to queue <short>."'
stub "$SBX" bsub "$BSUB_OK"
# bjobs: the job's state on each successive call, from STUB_JOBSTATE ("PEND RUN DONE").
stub "$SBX" bjobs '
read -r -a seq <<< "$STUB_JOBSTATE"
n=$(cat "$CALLS.bjobs" 2>/dev/null || echo 0); echo $((n + 1)) > "$CALLS.bjobs"
st=${seq[$n]:-${seq[-1]}}
echo "JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME     SUBMIT_TIME"
echo "$1    tester  $st  short      node001     node002     mb_preflight Nov  2 10:00"'
stub "$SBX" bkill 'echo "Job <$1> is being terminated"'
stub "$SBX" sleep ':'
stub "$SBX" lsid 'echo "IBM Spectrum LSF Standard 10.1.0.13, Apr 15 2022"; echo "My cluster name is hpcf"'
stub "$SBX" bparams 'printf "System default queues for automatic queue selection:\n    QUEUE_NAME = standard\n\n    RESOURCE_RESERVE_PER_TASK = Y\n    MBD_SLEEP_TIME = 10\n"'
printf 'LSF_SERVER_HOSTS="hpc1 hpc2"\nLSB_JOB_MEMLIMIT=Y\nLSB_RESOURCE_ENFORCE="memory cpu"\nLSF_UNIT_FOR_LIMITS=MB\n' > "$LSFENV/lsf.conf"
stub "$SBX" bqueues '
echo "QUEUE_NAME      PRIO STATUS          MAX JL/U JL/P JL/H NJOBS  PEND   RUN  SUSP"
for e in $STUB_QUEUES; do
  if [ "${e%%=*}" = "$1" ]; then printf "%-15s 30   %-15s -    -    -    -    0      0      0     0\n" "$1" "${e#*=}"; exit 0; fi
done
echo "$1: No such queue" >&2; exit 255'
stub "$SBX" java 'echo "openjdk version \"$STUB_JAVA\" 2024-10-15" >&2'
stub "$SBX" nextflow '[ "$1" = -version ] && printf "\n      N E X T F L O W\n      version %s build 1\n" "$STUB_NXF"'
stub "$SBX" nf-test 'echo "nf-test $STUB_NFT"'
stub "$SBX" tmux ':'
stub "$SBX" df 'printf "Filesystem 1024-blocks Used Available Capacity Mounted on\nfs 9999999999 1 %s 1%% /project\n" "$STUB_DF_KB"'
stub "$SBX" mmlsquota 'echo "Block Limits | File Limits: 812G used of 2T"'
# curl: STUB_NET=down cuts every host; STUB_CODES="host=code ..." sets one
# host's HTTP status (000: no connection); STUB_ENA=nobam makes ENA's
# filereport list no BAM. Prints the status for -w, ENA's TSV otherwise.
stub "$SBX" curl '
url=${!#}
host=${url#https://}; host=${host%%/*}
code=200
for e in $STUB_CODES; do [ "${e%%=*}" = "$host" ] && code=${e#*=}; done
[ "$STUB_NET" = down ] && code=000
case "$url" in
  *ghcr.io/token*) [ "$code" = 000 ] && exit 7; echo "{\"token\":\"anon\"}"; exit 0 ;;
esac
[ "$STUB_NET" = private ] && [ "$host" = ghcr.io ] && code=401
if [ "$code" = 000 ]; then printf 000; echo "curl: (7) Failed to connect to $host port 443" >&2; exit 7; fi
case "$url" in
  */ena/portal/api/filereport*)
    acc=$(sed -n "s/.*accession=\([A-Z0-9]*\).*/\1/p" <<< "$url")
    printf "run_accession\tsubmitted_ftp\tsubmitted_bytes\n"
    if [ "${STUB_ENA:-ok}" = nobam ]; then printf "%s\t\t\n" "$acc"
    else printf "%s\tftp.sra.ebi.ac.uk/vol1/run/SRR123/%s/4_1355_possorted_genome_bam.bam\t6400000000\n" "$acc" "$acc"; fi
    exit 0 ;;
esac
printf %s "$code"'
stub "$SING" singularity 'echo "singularity-ce version 4.1.1"'
cat > "$TMP/modinit.sh" <<EOF
module() {
  echo "module \$*" >> "\$CALLS"
  if [ "\$1" = load ] && [ "\$2" = "\${STUB_MODULE:-singularity/4.1.1}" ]; then PATH="$SING:\$PATH"; return 0; fi
  echo "Lmod has detected the following error: The following module(s) are unknown: \"\$2\"" >&2
  return 1
}
EOF

# pf [VAR=value ...]: run the preflight with good defaults (inside an LSF
# job), overridden by the arguments. PF_OPTS=--submit-test pf ... adds the option.
pf() {
  local opts=()
  read -r -a opts <<< "${PF_OPTS:-}"
  rm -f "$CALLS.bjobs"
  run env -i PATH="$SBX" HOME="$TMP/home" USER=tester CALLS="$CALLS" \
    MODULES_INIT="$TMP/modinit.sh" NXF_SINGULARITY_CACHEDIR="$TMP/cache" \
    LSB_JOBID=777 LSB_QUEUE=interactive LSF_ENVDIR="$LSFENV" \
    STUB_QUEUES="short=Open:Active standard=Open:Active large_mem=Open:Active" \
    STUB_JAVA=21.0.5 STUB_NXF=26.10.4 STUB_NFT=0.9.5 STUB_DF_KB=2000000000 STUB_NET=ok \
    STUB_JOBSTATE=DONE PREFLIGHT_POLL=1 \
    "$@" bash "$PF" ${opts[@]+"${opts[@]}"} "$WORK"
}
line() { grep -F -- "$1" <<< "$out" | head -n 1 || true; }   # the report line mentioning $1 (empty if none)
expect_line() { # expect_line <PASS|WARN|FAIL> <fragment> <name>
  local l; l=$(line "$2")
  if [[ "$l" == *"$1  "* ]]; then ok "$3"; else bad "$3" "expected $1 on the line with '$2'; got: '${l:-<none>}'"; fi
}

echo "hpcf_preflight.sh"

snapshot() { { find "$TMP/home" "$TMP/cache" "$TMP/project" -printf '%p %s %T@\n' | sort; git -C "$REPO" status --porcelain; } > "$1"; }
snapshot "$TMP/before"; : > "$CALLS"
pf
expect_rc 0 "everything in place: exit 0"
expect_out "summary: " "prints a summary line"
if ! grep -q '  FAIL  ' <<< "$out"; then ok "everything in place: no FAIL"; else bad "everything in place: no FAIL" "$out"; fi
snapshot "$TMP/after"
if cmp -s "$TMP/before" "$TMP/after"; then ok "changes nothing: home, cache, work directory and repository untouched"
else bad "changes nothing: home, cache, work directory and repository untouched" "$(diff "$TMP/before" "$TMP/after")"; fi
if ! grep -qE '^(bsub|bjobs|bkill) |^singularity (pull|build|exec|run)|^nextflow (run|pull|clone|self-update)' "$CALLS"; then
  ok "changes nothing: submits no job, pulls no image, runs no pipeline"
else bad "changes nothing: submits no job, pulls no image, runs no pipeline" "$(cat "$CALLS")"; fi

for q in short standard large_mem; do expect_line PASS "queue $q" "queue $q (read from conf/stjude.config): checked"; done
expect_line PASS "module load singularity/4.1.1" "singularity module from conf/stjude.config: loads"
expect_line PASS "ghcr.io/chialungchuang/mb-scrnaseq-r:" "pinned R image (nextflow.config): checked"
expect_line PASS "quay.io/biocontainers/star:" "STAR image (modules/): checked"
expect_line PASS "quay.io/biocontainers/sra-tools:" "sra-tools image (scripts/fetch_inputs.sh): checked"
expect_out "Block Limits" "quota tool output is shown"
expect_line PASS "inside LSF job 777" "inside an LSF job (interactive session): PASS"
expect_out "LSB_JOB_MEMLIMIT=Y" "lsf.conf (from LSF_ENVDIR): per-job memory limit setting shown"
expect_out 'LSB_RESOURCE_ENFORCE="memory cpu"' "lsf.conf: cgroup enforcement setting shown"
expect_out "RESOURCE_RESERVE_PER_TASK=Y" "bparams -a: reservation setting shown"
if ! grep -qE 'LSF_SERVER_HOSTS|MBD_SLEEP_TIME' <<< "$out"; then ok "lsf.conf and bparams: only the memory settings are shown"
else bad "lsf.conf and bparams: only the memory settings are shown" "$out"; fi
expect_line PASS "ENA portal API: lists the BAM of the pilot run SRR12354589" "ENA API: queried for the pilot run (samples.tsv run_order 1)"
expect_line PASS "ENA file server: reachable (https://ftp.sra.ebi.ac.uk/vol1/run/SRR123/SRR12354589/4_1355_possorted_genome_bam.bam)" "ENA file server: checked through the BAM URL ENA lists"
expect_line PASS "GEO (published matrix folder): reachable (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155446/suppl/)" "GEO folder (fetch_inputs.sh): checked"
expect_line PASS "10x reference GRCh38-3.0.0: reachable (https://cf.10xgenomics.com/" "10x reference (fetch_inputs.sh): checked"
expect_line PASS "STARsolo whitelists: reachable (https://labshare.cshl.edu/shares/dobin/dobin/STARsolo/10X/whitelists/737K-august-2016.txt.gz)" "whitelist (fetch_inputs.sh): checked"
# Downloads would be curl calls without -I to the file hosts.
if ! grep -E '^curl ' "$CALLS" | grep -E 'ftp\.sra\.ebi|ncbi\.nlm|10xgenomics|labshare' | grep -vqE '(^| )-I( |$)'; then
  ok "download hosts: HEAD requests only, nothing downloaded"
else bad "download hosts: HEAD requests only, nothing downloaded" "$(grep '^curl ' "$CALLS")"; fi

pf LSB_JOBID=
expect_line WARN "not inside an LSF job" "not inside an LSF job (login node): WARN"
pf LSF_ENVDIR="$TMP/nolsf"
expect_line WARN "lsf.conf not readable at $TMP/nolsf/lsf.conf" "no lsf.conf: WARN"
mkdir -p "$TMP/lsf2"; echo 'LSF_SERVER_HOSTS=hpc1' > "$TMP/lsf2/lsf.conf"
pf LSF_ENVDIR="$TMP/lsf2"
expect_out "none of LSB_JOB_MEMLIMIT" "lsf.conf without memory settings: says so"

pf STUB_QUEUES="short=Open:Active standard=Open:Active"
expect_rc 1 "a queue the profile uses is missing: exit 1"
expect_line FAIL "queue large_mem" "a queue the profile uses is missing: FAIL names it"
pf STUB_QUEUES="short=Open:Active standard=Closed:Adm large_mem=Open:Active"
expect_line WARN "queue standard" "a closed queue: WARN"

rm "$SBX/bsub"; pf; stub "$SBX" bsub "$BSUB_OK"
expect_line FAIL "bsub: not on PATH" "bsub missing: FAIL"

pf STUB_MODULE=singularity/3.8.0
expect_line FAIL "module load singularity/4.1.1" "the singularity module does not load: FAIL"
pf MODULES_INIT=/nonexistent
expect_line FAIL "module command not found" "no module system: FAIL"

pf STUB_JAVA=11.0.2
expect_line FAIL "java: '11.0.2'" "Java 11: FAIL (needs 17)"
pf STUB_JAVA=1.8.0_292
expect_line FAIL "java: '1.8.0_292'" "Java 1.8: FAIL"
pf STUB_JAVA=17.0.9
expect_line PASS "java: 17.0.9" "Java 17: PASS"

pf STUB_NXF=23.10.0
expect_line FAIL "nextflow: 23.10.0" "Nextflow 23.10 (the cluster module): FAIL"
pf STUB_NXF=9.10.0
expect_line FAIL "nextflow: 9.10.0" "Nextflow 9.10.0: FAIL (versions compared numerically, not as text)"
pf STUB_NXF=25.10.2
expect_line WARN "nextflow: 25.10.2" "Nextflow 25.10.2: WARN (works, not CI's version)"
rm "$SBX/nf-test"; pf; stub "$SBX" nf-test 'echo "nf-test $STUB_NFT"'
expect_line FAIL "nf-test: not on PATH" "nf-test missing: FAIL"

pf NXF_SINGULARITY_CACHEDIR=
expect_line FAIL "NXF_SINGULARITY_CACHEDIR: not set" "cache variable unset: FAIL"
pf NXF_SINGULARITY_CACHEDIR="$TMP/nope"
expect_line FAIL "no such directory" "cache directory missing: FAIL"
if [ "$(id -u)" -ne 0 ]; then
  mkdir -p "$TMP/ro"; chmod 555 "$TMP/ro"
  pf NXF_SINGULARITY_CACHEDIR="$TMP/ro"
  expect_line FAIL "not writable" "cache directory not writable: FAIL"
  chmod 755 "$TMP/ro"
else
  echo "  note  'cache not writable' needs a non-root user (root can write anywhere); CI runs it"
fi

pf STUB_DF_KB=5000000
expect_line FAIL "work directory $WORK: 4 GB free" "under 10 GB free: FAIL"
pf STUB_DF_KB=200000000
expect_line WARN "work directory $WORK: 190 GB free" "190 GB free: WARN (smoke test only)"
pf HOME="$TMP/project"
expect_line WARN "under \$HOME" "work directory under \$HOME: WARN"
rm "$SBX/mmlsquota"; pf; stub "$SBX" mmlsquota 'echo "Block Limits | File Limits: 812G used of 2T"'
expect_line WARN "quota: no mmlsquota" "no quota tool: WARN"

pf STUB_NET=down
expect_line FAIL "no connection to quay.io" "no outbound HTTPS: FAIL"
pf STUB_NET=private
expect_line FAIL "HTTP 401 (not anonymously pullable" "R image not public: FAIL"
pf STUB_CODES="www.ebi.ac.uk=000 ftp.sra.ebi.ac.uk=000"
expect_rc 1 "ENA blocked, registries open: exit 1"
expect_line FAIL "ENA portal API: no connection to www.ebi.ac.uk" "ENA API blocked: FAIL"
expect_line FAIL "ENA file server: no connection to ftp.sra.ebi.ac.uk" "ENA file server blocked: FAIL"
expect_line PASS "quay.io/biocontainers/star:" "ENA blocked: registries still checked on their own"
pf STUB_ENA=nobam
expect_line FAIL "ENA portal API: no BAM listed for the pilot run" "ENA lists no BAM: FAIL"
pf STUB_CODES="cf.10xgenomics.com=404"
expect_line FAIL "10x reference GRCh38-3.0.0: HTTP 404" "10x reference gone (HTTP 404): FAIL"

echo "hpcf_preflight.sh --submit-test"
: > "$CALLS"
PF_OPTS=--submit-test pf
expect_rc 0 "job runs: exit 0"
expect_line PASS "bsub from this node: job 4321 submitted" "submits from this node: PASS"
expect_line PASS "job 4321 ran and exited 0" "the job finishes (DONE): PASS"
if grep -qxF 'bsub -q short -W 1 -n 1 -J mb_preflight -o /dev/null -e /dev/null true' "$CALLS"; then
  ok "one 1-minute job running true on short, no output files"
else bad "one 1-minute job running true on short, no output files" "$(cat "$CALLS")"; fi
PF_OPTS=--submit-test pf STUB_JOBSTATE="PEND PEND RUN DONE"
expect_line PASS "job 4321 ran and exited 0 (seen after 3 s)" "pending, then running, then done: waits for it"
PF_OPTS=--submit-test pf STUB_JOBSTATE=EXIT
expect_rc 1 "job fails: exit 1"
expect_line FAIL "job 4321 exited with an error" "the job fails (EXIT): FAIL"
: > "$CALLS"
PF_OPTS=--submit-test pf STUB_JOBSTATE=PEND PREFLIGHT_SUBMIT_TIMEOUT=3
expect_line WARN "job 4321 still PEND after 3 s, killed" "the job never starts: WARN after the timeout"
if grep -qxF 'bkill 4321' "$CALLS"; then ok "the job never starts: it is killed"; else bad "the job never starts: it is killed" "$(cat "$CALLS")"; fi
PF_OPTS=--submit-test pf STUB_BSUB=refuse
expect_line FAIL "bsub from this node failed" "bsub refused: FAIL"
expect_out "User permission denied" "bsub refused: LSF's message is shown"
PF_OPTS=--bogus pf
expect_rc 2 "unknown option: exit 2"

echo "hpcf_preflight.sh --local (a server without a scheduler)"
stub "$SBX" nproc 'echo "${STUB_NPROC:-32}"'
printf 'MemTotal:       263868372 kB\nMemFree:        1000 kB\n' > "$TMP/meminfo_256"
printf 'MemTotal:        32775532 kB\n' > "$TMP/meminfo_32"
mkdir -p "$TMP/rt_sing" "$TMP/rt_app" "$TMP/rt_docker" "$TMP/rt_docker_down" "$TMP/rt_none"
stub "$TMP/rt_sing" singularity 'echo "singularity-ce version 4.1.1"'
stub "$TMP/rt_app" apptainer 'echo "apptainer version 1.3.4"'
stub "$TMP/rt_docker" docker '[ "$1" = info ] && echo "27.3.1"'
stub "$TMP/rt_docker_down" docker '[ "$1" = info ] && { echo "Cannot connect to the Docker daemon" >&2; exit 1; }'
lpf() { PF_OPTS=--local pf PATH="$SBX:$TMP/$1" LSB_JOBID= PREFLIGHT_MEMINFO="$TMP/meminfo_256" "${@:2}"; }

: > "$CALLS"
lpf rt_sing
expect_rc 0 "local, everything in place: exit 0"
expect_line PASS "memory: 251 GB" "memory read from /proc/meminfo: 251 GB passes"
expect_line PASS "cores: 32" "cores reported"
expect_line PASS "singularity: singularity-ce version 4.1.1" "Singularity found"
if ! grep -qE '^(Session|LSF|Singularity)$' <<< "$out" && ! grep -q 'not inside an LSF job' <<< "$out"; then
  ok "local: no session, LSF or module checks"
else bad "local: no session, LSF or module checks" "$out"; fi
if ! grep -qE '^(module|bsub|bjobs|bqueues|bparams|lsid) ' "$CALLS"; then ok "local: runs no LSF or module command"
else bad "local: runs no LSF or module command" "$(cat "$CALLS")"; fi
expect_line WARN "bsub is on PATH" "a scheduler on PATH: WARN (use its profile instead)"
expect_line PASS "ghcr.io/chialungchuang/mb-scrnaseq-r:" "local: registries still checked"
lpf rt_sing PREFLIGHT_MEMINFO="$TMP/meminfo_32"
expect_rc 1 "local, 31 GB of memory: exit 1"
expect_line FAIL "memory: 31 GB, below the 48 GB" "31 GB: FAIL (48 GB from nextflow.config)"
lpf rt_none
expect_line FAIL "no container runtime" "no container runtime: FAIL"
lpf rt_app NXF_SINGULARITY_CACHEDIR= NXF_APPTAINER_CACHEDIR="$TMP/cache"
expect_line PASS "apptainer: apptainer version 1.3.4" "Apptainer found"
expect_line PASS "NXF_APPTAINER_CACHEDIR=$TMP/cache: writable" "Apptainer only: its cache variable is the one checked"
lpf rt_docker NXF_SINGULARITY_CACHEDIR=
expect_line PASS "docker: server 27.3.1" "Docker with a reachable daemon: PASS"
expect_out "image cache: none needed" "Docker only: no image cache variable needed"
lpf rt_docker_down
expect_line WARN "docker: installed, but the daemon is not reachable" "Docker without its daemon: WARN"
expect_line FAIL "no container runtime" "Docker without its daemon and nothing else: FAIL"
lpf rt_sing STUB_NPROC=4
expect_line WARN "cores: 4" "4 cores: WARN (the STAR steps ask for 8)"
PF_OPTS="--local --submit-test" pf
expect_rc 2 "--local with --submit-test: exit 2"

finish
