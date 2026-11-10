#!/usr/bin/env bash
# Every R step stays within its CPUs: each module that runs in the R image
# exports the thread limits as the task's cpus before R starts, and each
# R script it runs logs them (tests/main.nf.test checks the logged values
# in a real run). This static check also covers a value hard-coded where
# ${task.cpus} belongs, which a 1-CPU test run cannot tell apart.
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"

echo "R steps: thread limits"
# shellcheck disable=SC2016  # the module text is matched literally
mods=$(grep -l 'container "${params.r_container}"' "$REPO"/modules/local/*/main.nf | sort)
n=$(wc -l <<< "$mods")
if [ "$n" -ge 4 ]; then ok "R modules found: $n"; else bad "R modules found: $n" "expected at least 4: $mods"; fi
for m in $mods; do
  name=$(basename "$(dirname "$m")")
  script=$(sed -n '/^    script:/,/^    stub:/p' "$m")
  missing=()
  for v in OMP_NUM_THREADS OMP_THREAD_LIMIT OPENBLAS_NUM_THREADS MKL_NUM_THREADS R_DATATABLE_NUM_THREADS; do
    grep -qE "export .*\\b$v=\\\$\\{task\\.cpus\\}" <<< "$script" || missing+=("$v")
  done
  if [ ${#missing[@]} -eq 0 ]; then ok "$name: exports every thread limit as \${task.cpus}"
  else bad "$name: exports every thread limit as \${task.cpus}" "missing or not \${task.cpus}: ${missing[*]}"; fi
  while read -r r; do
    if grep -qx 'log_threads()' "$REPO/bin/$r"; then ok "$name: bin/$r logs its threads"
    else bad "$name: bin/$r logs its threads" "no log_threads() line in bin/$r"; fi
  done < <(grep -oE '[a-z_]+\.R ' <<< "$script" | sort -u)
done
finish
