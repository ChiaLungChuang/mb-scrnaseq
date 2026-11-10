#!/usr/bin/env bash
# The standalone-server profile (conf/server.config): --server_cpus and
# --server_memory cap every task (process.resourceLimits) and what runs at
# once (the local executor's cpus and memory), and the profile combines
# with each container profile. A tiny pipeline under conf/base.config +
# conf/server.config reports what each task was given.
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
command -v nextflow >/dev/null || { echo "nextflow is needed for these tests"; exit 1; }

cat > "$TMP/main.nf" <<'EOF'
params.part = 'caps'
process HIGH {
    label 'process_high'
    output: stdout
    script: "echo HIGH ${task.cpus} ${task.memory}"
}
process RSTEP {
    label 'process_r'
    output: stdout
    script: "echo RSTEP ${task.cpus} ${task.memory}"
}
process LOW {
    label 'process_low'
    output: stdout
    script: "echo LOW ${task.cpus} ${task.memory}"
}
// Two tasks that fit the share's memory but not its cores together...
process CPU_PAIR {
    cpus 2
    memory 100.MB
    input: val i
    output: stdout
    script: "start=\$(date +%s.%N); sleep 2; echo PAIR \$start \$(date +%s.%N)"
}
// ...and two that fit its cores but not its memory together.
process MEM_PAIR {
    cpus 1
    memory 2.GB
    input: val i
    output: stdout
    script: "start=\$(date +%s.%N); sleep 2; echo PAIR \$start \$(date +%s.%N)"
}
workflow {
    if (params.part == 'caps') {
        HIGH().view()
        RSTEP().view()
        LOW().view()
    } else if (params.part == 'cpu_pair') {
        CPU_PAIR(channel.of(1, 2)).view()
    } else {
        MEM_PAIR(channel.of(1, 2)).view()
    }
}
EOF

run_server() { # run_server [nextflow args...]: the tiny pipeline under the server profile's config
  (cd "$TMP" && NXF_OFFLINE=true nextflow run main.nf -c "$REPO/conf/base.config" -c "$REPO/conf/server.config" "$@" 2>&1)
}
line() { grep -m1 "^$1 " <<< "$out" || true; }
one_after_other() { # one_after_other <name>: the two PAIR lines in $out do not overlap in time
  local s1 e1 s2 e2
  read -r s1 e1 s2 e2 <<< "$(grep '^PAIR ' <<< "$out" | awk '{ print $2, $3 }' | sort -n | tr '\n' ' ')"
  if [ -n "${s2:-}" ] && awk -v e1="$e1" -v s2="$s2" 'BEGIN { exit !(s2 >= e1) }'; then ok "$1"
  else bad "$1" "starts/ends: ${s1:-?} ${e1:-?} ${s2:-?} ${e2:-?}; output: $out"; fi
}

echo "server profile: config"

run run_server --server_cpus 2 --server_memory '3.GB'
expect_rc 0 "runs with --server_cpus 2 --server_memory 3.GB"
if [ "$(line HIGH)" = "HIGH 2 3 GB" ]; then ok "an 8-CPU, 48 GB step gets 2 CPUs and 3 GB"; else bad "an 8-CPU, 48 GB step gets 2 CPUs and 3 GB" "got: '$(line HIGH)'"; fi
if [ "$(line RSTEP)" = "RSTEP 1 3 GB" ]; then ok "a 12 GB R step gets 3 GB"; else bad "a 12 GB R step gets 3 GB" "got: '$(line RSTEP)'"; fi
if [ "$(line LOW)" = "LOW 1 3 GB" ]; then ok "a 4 GB step gets 3 GB"; else bad "a 4 GB step gets 3 GB" "got: '$(line LOW)'"; fi
run run_server --server_cpus 2 --server_memory '3.GB' --part cpu_pair
one_after_other "two 2-CPU tasks on a 2-CPU share run one after the other"
run run_server --server_cpus 2 --server_memory '3.GB' --part mem_pair
one_after_other "two 2 GB tasks on a 3 GB share run one after the other"

for c in singularity apptainer docker; do
  out=$(cd "$REPO" && NXF_OFFLINE=true nextflow config -profile "server,$c" -o flat 2>&1)
  if grep -qx "$c.enabled = true" <<< "$out" && grep -qx "process.executor = 'local'" <<< "$out"; then
    ok "-profile server,$c parses: local executor, $c on"
  else
    bad "-profile server,$c parses: local executor, $c on" "$(grep -E "enabled|executor|ERROR|rror" <<< "$out" | head -5)"
  fi
done
out=$(cd "$REPO" && NXF_OFFLINE=true nextflow config -profile server,apptainer -o flat 2>&1)
if ! grep -qx "singularity.enabled = true" <<< "$out" && ! grep -qx "docker.enabled = true" <<< "$out"; then
  ok "-profile server,apptainer leaves Singularity and Docker off"
else
  bad "-profile server,apptainer leaves Singularity and Docker off" "$(grep enabled <<< "$out")"
fi

finish
