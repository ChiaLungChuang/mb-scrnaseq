#!/usr/bin/env bash
# Checks to run before the smoke test (docs/RUNBOOK.md, section 3). Prints
# one report to paste back.
#
# On the cluster, run it on a compute node, inside an interactive session:
# HPCF's login nodes are only for changing directories and submitting jobs,
# and the network checks must describe the nodes that will do the work.
#   hpcf_interactive          (or: bsub -P interactive -q interactive -Is bash)
#   scripts/hpcf_preflight.sh [--submit-test] [work_dir]
# On one server without a scheduler (-profile server; runbook, "Route 2"):
#   scripts/hpcf_preflight.sh --local [work_dir]
#
#   --local        a standalone server: skip the LSF and module checks;
#                  check the machine's cores and memory (FAIL below
#                  server_min_memory in nextflow.config, 48 GB) and that
#                  Singularity, Apptainer or Docker is there
#   --submit-test  also submit one job from this node (`true`, queue short,
#                  1-minute limit, no output files) and wait for it to
#                  finish: Nextflow's head process submits every task from
#                  the node it runs on. Waits up to PREFLIGHT_SUBMIT_TIMEOUT
#                  seconds (default 900), then kills the job. Not with
#                  --local.
#   work_dir       where Nextflow's work/ directory will live (default: the
#                  current directory)
#
# Checks: that it runs inside an LSF job; LSF commands and the queues
# conf/stjude.config routes to; lsf.conf's memory settings (reported for
# HPCF to interpret, not judged); the Singularity module the profile loads;
# git, Java, Nextflow and nf-test versions; NXF_SINGULARITY_CACHEDIR; free
# space and quota for the work directory; anonymous access to every
# container image the pipeline uses (manifest requests to ghcr.io and
# quay.io); the download hosts of scripts/fetch_inputs.sh (one ENA
# filereport query for the pilot run, then HEAD requests for its BAM, the
# GEO folder, the 10x reference and a whitelist; nothing is downloaded).
# The SRA prefetch fallback's hosts are not checked.
#
# Without --submit-test it changes nothing: no file is written, no image is
# pulled, no job is submitted, and the module is loaded only inside a
# subshell. With it, the one job is the only change. Queue names, the
# module, the images, the download URLs and the versions are read from this
# repository's own files, so the checks follow them.
#
# Exit status: 1 if any check FAILs, 2 for a usage error, otherwise 0.
set -uo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
submit_test=0 local_mode=0
while [ $# -gt 0 ]; do
  case "$1" in
    --submit-test) submit_test=1; shift ;;
    --local) local_mode=1; shift ;;
    -*) echo "unknown option: $1 (usage: $0 [--local | --submit-test] [work_dir])" >&2; exit 2 ;;
    *) break ;;
  esac
done
if [ "$local_mode" = 1 ] && [ "$submit_test" = 1 ]; then
  echo "--submit-test is for the cluster; a server without a scheduler has nothing to submit to" >&2; exit 2
fi
work_dir=${1:-$PWD}

npass=0 nwarn=0 nfail=0
pass() { npass=$((npass + 1)); printf '  PASS  %s\n' "$*"; }
warn() { nwarn=$((nwarn + 1)); printf '  WARN  %s\n' "$*"; }
fail() { nfail=$((nfail + 1)); printf '  FAIL  %s\n' "$*"; }
info() { printf '  info  %s\n' "$*"; }
section() { printf '\n%s\n' "$*"; }
indent() { sed 's/^/        /'; }

# version_ge A B: true if version A >= version B (numeric, field by field).
version_ge() { [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n 1)" = "$2" ]; }

# Free space in whole GB on the filesystem holding a path.
free_gb() { df -Pk "$1" 2>/dev/null | awk 'NR == 2 { printf "%d", $4 / 1048576 }'; }

# --- what the repository's config says -------------------------------------
queues=$(sed -n '/queue = {/,/^    }/p' "$REPO/conf/stjude.config" | grep -oP "'\K[A-Za-z0-9_]+(?=')" | sort -u)
sing_module=$(grep -oP 'module load \K\S+' "$REPO/conf/stjude.config" | head -n 1)
nxf_min=$(grep -oP "nextflowVersion\s*=\s*'>=\K[^']+" "$REPO/nextflow.config")
nxf_ci=$(grep -oP "NXF_VER: '\K[^']+" "$REPO/.github/workflows/ci.yml")
nft_ci=$(grep -oP "NFT_VER: '\K[^']+" "$REPO/.github/workflows/ci.yml")
images=$( {
  grep -oP "^\s*r_container\s*=\s*'\K[^']+" "$REPO/nextflow.config"
  grep -hoP "container\s+'\K[^']+" "$REPO"/modules/local/*/main.nf
  grep -oP '^\w+_IMAGE=docker://\K\S+' "$REPO/scripts/fetch_inputs.sh"
} | sort -u)

echo "mb-scrnaseq HPCF preflight"
info "host $(hostname 2>/dev/null || echo '?'), user ${USER:-?}, $(date -u +%Y-%m-%dT%H:%M:%SZ)"
info "repository $REPO at $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo '?')"
info "work directory $work_dir"

# --- The cluster: session, LSF, the Singularity module -------------------------
cluster_checks() {
section "Session"
if [ -n "${LSB_JOBID:-}" ]; then
  pass "inside LSF job $LSB_JOBID (queue ${LSB_QUEUE:-?}), not on a login node"
else
  warn "not inside an LSF job: HPCF's login nodes are only for cd and bsub. Run this inside hpcf_interactive (or bsub -P interactive -q interactive -Is bash); the network checks below describe this node only"
fi

# --- LSF ---------------------------------------------------------------------
section "LSF"
if command -v lsid >/dev/null; then info "$(lsid 2>&1 | head -n 1)"; fi
for c in bsub bjobs bqueues; do
  if p=$(command -v "$c"); then pass "$c: $p"; else fail "$c: not on PATH"; fi
done
if command -v bqueues >/dev/null; then
  for q in $queues; do
    out=$(bqueues "$q" 2>&1)
    status=$(awk -v q="$q" '$1 == q { print $3 }' <<< "$out")
    if [ -z "$status" ]; then
      fail "queue $q (used by conf/stjude.config): not visible: $(head -n 1 <<< "$out")"
    elif [ "$status" = "Open:Active" ]; then
      pass "queue $q: $status"
    else
      warn "queue $q: $status (jobs may wait or be refused)"
    fi
  done
fi

# How LSF applies -M and rusage[mem] (conf/stjude.config, "Memory"). Shown,
# not judged: what the values mean for a multi-threaded STAR is HPCF's call.
mem_keys='LSB_JOB_MEMLIMIT|LSB_MEMLIMIT_ENFORCE|LSB_MEMLIMIT_ENF_CONTROL|LSB_RESOURCE_ENFORCE|LSF_UNIT_FOR_LIMITS|RESOURCE_RESERVE_PER_TASK|RESOURCE_RESERVE_PER_SLOT'
mem_settings() { grep -E "^[[:space:]]*($mem_keys)[[:space:]]*=" | sed -E 's/^[[:space:]]+//; s/[[:space:]]*=[[:space:]]*/=/'; }
lsf_conf=${LSF_ENVDIR:-/etc}/lsf.conf
if [ -r "$lsf_conf" ]; then
  found=$(mem_settings < "$lsf_conf")
  if [ -n "$found" ]; then
    info "memory settings in $lsf_conf:"
    indent <<< "$found"
  else
    info "memory settings in $lsf_conf: none of ${mem_keys//|/, } set (LSF defaults)"
  fi
else
  warn "lsf.conf not readable at $lsf_conf (LSF_ENVDIR=${LSF_ENVDIR:-unset}): ask HPCF how -M is enforced for a multi-threaded job"
fi
if command -v bparams >/dev/null; then
  found=$(bparams -a 2>/dev/null | mem_settings)
  if [ -n "$found" ]; then info "memory settings from bparams -a:"; indent <<< "$found"; fi
fi

# --- Singularity --------------------------------------------------------------
section "Singularity"
if ! type module >/dev/null 2>&1; then
  for f in ${MODULES_INIT:-} /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
           /usr/share/lmod/lmod/init/bash "${MODULESHOME:-/nonexistent}/init/bash"; do
    # shellcheck disable=SC1090  # site-specific init file
    if [ -r "$f" ]; then . "$f" >/dev/null 2>&1; break; fi
  done
fi
if ! type module >/dev/null 2>&1; then
  fail "module command not found (needed for 'module load $sing_module' in every job)"
else
  # In a subshell: the module changes nothing outside this check.
  if out=$( (module load "$sing_module" && singularity --version) 2>&1 ); then
    pass "module load $sing_module: $(tail -n 1 <<< "$out")"
  else
    fail "module load $sing_module failed:"
    indent <<< "$out"
  fi
fi
}

# --- One server without a scheduler (--local) ---------------------------------
machine_checks() {
  section "Machine (no scheduler: -profile server)"
  local cores mem_kb mem_gb min_gb found=0
  min_gb=$(grep -oP "^\s*server_min_memory\s*=\s*'\K[0-9]+" "$REPO/nextflow.config")
  cores=$(nproc 2>/dev/null || true)
  if [ -z "$cores" ]; then warn "cores: unknown (no nproc)"
  elif [ "$cores" -lt 8 ]; then warn "cores: $cores (the STAR steps ask for 8; they will be capped to --server_cpus and run slower)"
  else pass "cores: $cores (give the run its share with --server_cpus)"; fi
  mem_kb=$(awk '/^MemTotal:/ { print $2 }' "${PREFLIGHT_MEMINFO:-/proc/meminfo}" 2>/dev/null)
  if [ -z "$mem_kb" ]; then
    warn "memory: unknown (no ${PREFLIGHT_MEMINFO:-/proc/meminfo})"
  else
    mem_gb=$((mem_kb / 1048576))
    if [ "$mem_gb" -ge "${min_gb:-48}" ]; then pass "memory: $mem_gb GB (the STAR steps ask for ${min_gb:-48} GB; give the run its share with --server_memory)"
    else fail "memory: $mem_gb GB, below the ${min_gb:-48} GB the STAR index build and STARsolo ask for"; fi
  fi
  for r in singularity apptainer; do
    if command -v "$r" >/dev/null; then pass "$r: $("$r" --version 2>&1 | head -n 1) (-profile server,$r)"; found=1; fi
  done
  if command -v docker >/dev/null; then
    if v=$(docker info --format '{{.ServerVersion}}' 2>/dev/null) && [ -n "$v" ]; then pass "docker: server $v (-profile server,docker)"; found=1
    else warn "docker: installed, but the daemon is not reachable as $(id -un 2>/dev/null || echo this user)"; fi
  fi
  [ "$found" = 1 ] || fail "no container runtime: Singularity, Apptainer or a usable Docker is needed"
  for s in bsub sbatch qsub; do
    if command -v "$s" >/dev/null; then warn "$s is on PATH: if this machine is shared through a scheduler, use its profile rather than -profile server"; fi
  done
}

if [ "$local_mode" = 1 ]; then machine_checks; else cluster_checks; fi

# --- Tools --------------------------------------------------------------------
section "Tools"
if command -v git >/dev/null; then pass "git: $(git --version)"; else fail "git: not on PATH (the pre-registration check needs it)"; fi
if command -v tmux >/dev/null; then pass "tmux: present"; else warn "tmux: not on PATH (the runbook keeps the interactive session in tmux)"; fi

if command -v java >/dev/null; then
  jv=$(java -version 2>&1 | awk -F'"' '/version/ { print $2; exit }')
  jmajor=${jv%%.*}
  if [ "$jmajor" = 1 ]; then jmajor=$(cut -d. -f2 <<< "$jv"); fi   # 1.8.0_292 -> 8
  if [[ "$jmajor" =~ ^[0-9]+$ ]] && [ "$jmajor" -ge 17 ]; then pass "java: $jv"
  else fail "java: '${jv:-unknown}' (Nextflow needs Java 17 or newer)"; fi
else
  fail "java: not on PATH"
fi

if command -v nextflow >/dev/null; then
  # NXF_OFFLINE: never let the launcher download anything.
  out=$(NXF_OFFLINE=true nextflow -version 2>&1)
  nv=$(grep -oP 'version \K[0-9]+\.[0-9]+\.[0-9]+' <<< "$out" | head -n 1)
  if [ -z "$nv" ]; then
    fail "nextflow: version not reported (offline check; run 'nextflow -version' once):"
    head -n 3 <<< "$out" | indent
  elif ! version_ge "$nv" "$nxf_min"; then
    fail "nextflow: $nv (this pipeline needs >= $nxf_min)"
  elif [ "$nv" != "$nxf_ci" ]; then
    warn "nextflow: $nv (works; CI tests $nxf_ci)"
  else
    pass "nextflow: $nv"
  fi
else
  fail "nextflow: not on PATH (runbook section 2: micromamba activate nf)"
fi

if command -v nf-test >/dev/null; then
  tv=$(nf-test version 2>&1 | grep -oP 'nf-test \K[0-9]+\.[0-9]+\.[0-9]+' | head -n 1)
  if [ "$tv" = "$nft_ci" ]; then pass "nf-test: $tv"; else warn "nf-test: '${tv:-unknown}' (CI tests $nft_ci)"; fi
else
  fail "nf-test: not on PATH (runbook section 2)"
fi

# --- Storage ------------------------------------------------------------------
section "Storage"
# The image cache variable follows the container runtime: Singularity on
# the cluster; on a server whichever of Singularity or Apptainer is there
# (Docker keeps its own images).
cache_var=NXF_SINGULARITY_CACHEDIR
if [ "$local_mode" = 1 ] && ! command -v singularity >/dev/null; then
  if command -v apptainer >/dev/null; then cache_var=NXF_APPTAINER_CACHEDIR; else cache_var=; fi
fi
cache=${cache_var:+${!cache_var:-}}
if [ -z "$cache_var" ]; then
  info "image cache: none needed (Docker keeps its own images)"
elif [ -z "$cache" ]; then
  fail "$cache_var: not set (runbook section 2)"
elif [ ! -d "$cache" ]; then
  fail "$cache_var=$cache: no such directory"
elif [ ! -w "$cache" ]; then
  fail "$cache_var=$cache: not writable"
else
  g=$(free_gb "$cache")
  if [ -n "$g" ] && [ "$g" -lt 10 ]; then warn "$cache_var=$cache: writable, only ${g} GB free (the images need a few GB)"
  else pass "$cache_var=$cache: writable, ${g:-?} GB free"; fi
fi

if [ ! -d "$work_dir" ]; then
  fail "work directory $work_dir: no such directory"
elif [ ! -w "$work_dir" ]; then
  fail "work directory $work_dir: not writable"
else
  g=$(free_gb "$work_dir")
  if [ -z "$g" ]; then warn "work directory $work_dir: free space unknown (df failed)"
  elif [ "$g" -lt 10 ]; then fail "work directory $work_dir: ${g} GB free (the smoke test alone needs ~10 GB)"
  elif [ "$g" -lt 500 ]; then warn "work directory $work_dir: ${g} GB free (enough for the smoke test; 8 real samples need a few hundred GB: BAMs, FASTQs of about the same size, a 30 GB STAR index)"
  else pass "work directory $work_dir: ${g} GB free"; fi
  case "$work_dir/" in
    "${HOME:-/nonexistent}"/*) warn "work directory is under \$HOME (runbook: use project or scratch storage)" ;;
  esac
fi

mount=$(df -P "$work_dir" 2>/dev/null | awk 'NR == 2 { print $6 }')
if command -v mmlsquota >/dev/null; then
  info "quota (mmlsquota):"; mmlsquota -u "${USER:-$(id -un)}" --block-size auto 2>&1 | indent
elif command -v lfs >/dev/null; then
  info "quota (lfs quota on ${mount:-?}):"; lfs quota -h -u "${USER:-$(id -un)}" "${mount:-$work_dir}" 2>&1 | indent
elif command -v quota >/dev/null; then
  info "quota (quota -s):"; quota -s 2>&1 | indent
else
  warn "quota: no mmlsquota, lfs or quota command; check your quota for ${mount:-$work_dir} with HPCF"
fi

# --- Network: anonymous manifest requests, nothing is pulled -----------------
section "Container registries (anonymous; from this node)"
accept='application/vnd.oci.image.index.v1+json, application/vnd.docker.distribution.manifest.list.v2+json, application/vnd.oci.image.manifest.v1+json, application/vnd.docker.distribution.manifest.v2+json'
for img in $images; do
  registry=${img%%/*}
  rest=${img#*/}
  if [[ "$rest" == *@* ]]; then repo=${rest%%@*}; ref=${rest#*@}; repo=${repo%:*}
  else repo=${rest%:*}; ref=${rest##*:}; fi
  auth=()
  if [ "$registry" = ghcr.io ]; then
    token=$(curl -sS --max-time 20 "https://ghcr.io/token?scope=repository:${repo}:pull" 2>/dev/null | sed -n 's/.*"token" *: *"\([^"]*\)".*/\1/p')
    [ -n "$token" ] && auth=(-H "Authorization: Bearer $token")
  fi
  code=$(curl -sS -I -o /dev/null -w '%{http_code}' --max-time 20 ${auth[@]+"${auth[@]}"} -H "Accept: $accept" \
           "https://${registry}/v2/${repo}/manifests/${ref}" 2>/dev/null)
  case "$code" in
    200) pass "$img: reachable" ;;
    000|'') fail "$img: no connection to $registry (outbound HTTPS blocked from this node?)" ;;
    401|403) fail "$img: HTTP $code (not anonymously pullable: is the package public?)" ;;
    *) fail "$img: HTTP $code" ;;
  esac
done

# --- Download hosts: one ENA query, then HEAD requests; nothing is downloaded --
section "Download hosts (scripts/fetch_inputs.sh; from this node)"
fi_sh="$REPO/scripts/fetch_inputs.sh"
# head_check <label> <url>: one HEAD request, following redirects.
head_check() {
  local code
  if [ -z "$2" ]; then fail "$1: URL not found in scripts/fetch_inputs.sh"; return; fi
  code=$(curl -sS -I -L -o /dev/null -w '%{http_code}' --max-time 30 "$2" 2>/dev/null)
  case "$code" in
    200) pass "$1: reachable ($2)" ;;
    000|'') local h=${2#https://}; fail "$1: no connection to ${h%%/*} (outbound HTTPS blocked from this node?)" ;;
    *) fail "$1: HTTP $code for $2" ;;
  esac
}
ena_api=$(grep -oP '^ENA_API=\K\S+' "$fi_sh")
pilot=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
                    $col["run_order"] == 1 { print $col["srr"]; exit }' "$REPO/docs/results/phase1/samples.tsv")
bam_url=''
if [ -z "$ena_api" ] || [ -z "$pilot" ]; then
  fail "ENA portal API: could not read ENA_API from scripts/fetch_inputs.sh or the pilot run from docs/results/phase1/samples.tsv"
elif ! rec=$(curl -sS --max-time 30 "$ena_api?accession=$pilot&result=read_run&fields=run_accession,submitted_ftp,submitted_bytes&format=tsv" 2>&1); then
  h=${ena_api#https://}
  fail "ENA portal API: no connection to ${h%%/*} (outbound HTTPS blocked from this node?)"
else
  IFS=$'\t' read -r acc ftps _ <<< "$(sed -n 2p <<< "$rec")" || true
  for f in ${ftps//;/ }; do case "$f" in *.bam) bam_url="https://$f" ;; esac; done
  if [ "${acc:-}" = "$pilot" ] && [ -n "$bam_url" ]; then
    pass "ENA portal API: lists the BAM of the pilot run $pilot"
  else
    fail "ENA portal API: no BAM listed for the pilot run $pilot; response: $(head -c 200 <<< "$rec" | tr '\n\t' '  ')"
  fi
fi
# The ENA file server, through the pilot's own BAM when ENA named it.
head_check "ENA file server" "${bam_url:-https://ftp.sra.ebi.ac.uk/vol1/}"
head_check "GEO (published matrix folder)" "$(grep -oP '^\s*base=\Khttps://\S+' "$fi_sh")/"
head_check "10x reference GRCh38-3.0.0" "$(grep -oP 'GRCh38-3\.0\.0\)\s+url=\K\S+' "$fi_sh")"
wl=$(grep -oP '^\s*wl=\K\S+' "$fi_sh")
# shellcheck disable=SC2016  # the literal text $wl in fetch_inputs.sh
head_check "STARsolo whitelists" "${wl:+$wl/$(grep -oP 'fetch "\$wl/\K[^"]+' "$fi_sh" | head -n 1)}"

# --- Job submission from this node (opt-in) ---------------------------------
if [ "$submit_test" = 1 ]; then
  section "Job submission from this node (--submit-test)"
  timeout_s=${PREFLIGHT_SUBMIT_TIMEOUT:-900}
  poll_s=${PREFLIGHT_POLL:-15}
  if ! command -v bsub >/dev/null; then
    fail "bsub: not on PATH, cannot submit"
  else
    out=$(bsub -q short -W 1 -n 1 -J mb_preflight -o /dev/null -e /dev/null true 2>&1)
    jobid=$(grep -oP 'Job <\K[0-9]+(?=>)' <<< "$out" | head -n 1)
    if [ -z "$jobid" ]; then
      fail "bsub from this node failed:"
      indent <<< "$out"
    else
      pass "bsub from this node: job $jobid submitted (queue short)"
      waited=0
      while :; do
        state=$(bjobs "$jobid" 2>/dev/null | awk -v j="$jobid" '$1 == j { print $3 }')
        case "$state" in
          DONE) pass "job $jobid ran and exited 0 (seen after ${waited} s)"; break ;;
          EXIT) fail "job $jobid exited with an error (see bjobs -l $jobid, bhist -l $jobid)"; break ;;
        esac
        if [ "$waited" -ge "$timeout_s" ]; then
          bkill "$jobid" >/dev/null 2>&1 || true
          warn "job $jobid still ${state:-not listed} after ${timeout_s} s, killed: submission works, but it did not run (queue busy? re-run later)"
          break
        fi
        sleep "$poll_s"
        waited=$((waited + poll_s))
      done
    fi
  fi
fi

printf '\nsummary: %d pass, %d warn, %d fail\n' "$npass" "$nwarn" "$nfail"
[ "$nfail" -eq 0 ]
