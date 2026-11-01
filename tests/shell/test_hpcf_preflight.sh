#!/usr/bin/env bash
# Tests for scripts/hpcf_preflight.sh. LSF, Java, Nextflow, nf-test, the
# module system, curl and df are stand-ins; the script runs with `env -i`
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
mkdir -p "$SBX" "$SING" "$TMP/home" "$TMP/cache" "$WORK"
for t in bash sed awk grep sort head tail cut date dirname git id hostname cat tr; do
  ln -s "$(command -v "$t")" "$SBX/$t"
done

stub() { # stub <dir> <name> <body>: a logging stand-in command
  printf '#!/bin/bash\necho "%s $*" >> "$CALLS"\n%s\n' "$2" "$3" > "$1/$2"
  chmod +x "$1/$2"
}
stub "$SBX" bsub 'echo "Job <1> is submitted"'
stub "$SBX" bjobs 'echo "No unfinished job found"'
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
stub "$SBX" curl '
url=${!#}
case "$url" in
  *ghcr.io/token*) [ "$STUB_NET" = down ] && exit 7; echo "{\"token\":\"anon\"}"; exit 0 ;;
esac
case "$STUB_NET" in
  down) printf 000; exit 7 ;;
  private) case "$url" in *ghcr.io*) printf 401 ;; *) printf 200 ;; esac ;;
  *) printf 200 ;;
esac'
stub "$SING" singularity 'echo "singularity-ce version 4.1.1"'
cat > "$TMP/modinit.sh" <<EOF
module() {
  echo "module \$*" >> "\$CALLS"
  if [ "\$1" = load ] && [ "\$2" = "\${STUB_MODULE:-singularity/4.1.1}" ]; then PATH="$SING:\$PATH"; return 0; fi
  echo "Lmod has detected the following error: The following module(s) are unknown: \"\$2\"" >&2
  return 1
}
EOF

# pf [VAR=value ...]: run the preflight with good defaults, overridden by the arguments.
pf() {
  run env -i PATH="$SBX" HOME="$TMP/home" USER=tester CALLS="$CALLS" \
    MODULES_INIT="$TMP/modinit.sh" NXF_SINGULARITY_CACHEDIR="$TMP/cache" \
    STUB_QUEUES="short=Open:Active standard=Open:Active large_mem=Open:Active" \
    STUB_JAVA=21.0.5 STUB_NXF=26.10.4 STUB_NFT=0.9.5 STUB_DF_KB=2000000000 STUB_NET=ok \
    "$@" bash "$PF" "$WORK"
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
if ! grep -qE '^(bsub|bjobs) |^singularity (pull|build|exec|run)|^nextflow (run|pull|clone|self-update)' "$CALLS"; then
  ok "changes nothing: submits no job, pulls no image, runs no pipeline"
else bad "changes nothing: submits no job, pulls no image, runs no pipeline" "$(cat "$CALLS")"; fi

for q in short standard large_mem; do expect_line PASS "queue $q" "queue $q (read from conf/stjude.config): checked"; done
expect_line PASS "module load singularity/4.1.1" "singularity module from conf/stjude.config: loads"
expect_line PASS "ghcr.io/chialungchuang/mb-scrnaseq-r:" "pinned R image (nextflow.config): checked"
expect_line PASS "quay.io/biocontainers/star:" "STAR image (modules/): checked"
expect_line PASS "quay.io/biocontainers/sra-tools:" "sra-tools image (scripts/fetch_inputs.sh): checked"
expect_out "Block Limits" "quota tool output is shown"

pf STUB_QUEUES="short=Open:Active standard=Open:Active"
expect_rc 1 "a queue the profile uses is missing: exit 1"
expect_line FAIL "queue large_mem" "a queue the profile uses is missing: FAIL names it"
pf STUB_QUEUES="short=Open:Active standard=Closed:Adm large_mem=Open:Active"
expect_line WARN "queue standard" "a closed queue: WARN"

rm "$SBX/bsub"; pf; stub "$SBX" bsub 'echo "Job <1> is submitted"'
expect_line FAIL "bsub" "bsub missing: FAIL"

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

finish
