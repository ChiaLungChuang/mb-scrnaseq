#!/usr/bin/env bash
# Read-only checks to run on an HPCF login node before the smoke test
# (docs/RUNBOOK.md, section 3). Prints one report to paste back.
#
# Usage: scripts/hpcf_preflight.sh [work_dir]
#   work_dir  where Nextflow's work/ directory will live (default: the
#             current directory)
#
# Checks: LSF commands and the queues conf/stjude.config routes to; the
# Singularity module it loads; git, Java, Nextflow and nf-test versions;
# NXF_SINGULARITY_CACHEDIR; free space and quota for the work directory;
# anonymous access to every container image the pipeline uses (manifest
# requests to ghcr.io and quay.io only).
#
# It changes nothing: no file is written, no image is pulled, no job is
# submitted, and the module is loaded only inside a subshell. Queue names,
# the module, the images and the versions are read from this repository's
# own config, so the checks follow the config.
#
# Exit status: 1 if any check FAILs, otherwise 0.
set -uo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
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

# --- LSF ---------------------------------------------------------------------
section "LSF"
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

# --- Tools --------------------------------------------------------------------
section "Tools"
if command -v git >/dev/null; then pass "git: $(git --version)"; else fail "git: not on PATH (the pre-registration check needs it)"; fi
if command -v tmux >/dev/null; then pass "tmux: present"; else warn "tmux: not on PATH (the runbook runs Nextflow inside tmux)"; fi

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
cache=${NXF_SINGULARITY_CACHEDIR:-}
if [ -z "$cache" ]; then
  fail "NXF_SINGULARITY_CACHEDIR: not set (runbook section 2)"
elif [ ! -d "$cache" ]; then
  fail "NXF_SINGULARITY_CACHEDIR=$cache: no such directory"
elif [ ! -w "$cache" ]; then
  fail "NXF_SINGULARITY_CACHEDIR=$cache: not writable"
else
  g=$(free_gb "$cache")
  if [ -n "$g" ] && [ "$g" -lt 10 ]; then warn "NXF_SINGULARITY_CACHEDIR=$cache: writable, only ${g} GB free (the images need a few GB)"
  else pass "NXF_SINGULARITY_CACHEDIR=$cache: writable, ${g:-?} GB free"; fi
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
section "Container registries (anonymous; from this login node)"
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

printf '\nsummary: %d pass, %d warn, %d fail\n' "$npass" "$nwarn" "$nfail"
[ "$nfail" -eq 0 ]
