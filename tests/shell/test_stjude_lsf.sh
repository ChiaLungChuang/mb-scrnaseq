#!/usr/bin/env bash
# The LSF directives the St. Jude profile makes Nextflow write, checked
# against HPCF's documented policy: memory is requested per core (LSF
# multiplies rusage[mem] by the slots), short has a 30-minute hard limit,
# large_mem is for more than 1 TB, standard is the default. Nextflow runs
# a tiny pipeline under conf/base.config + conf/stjude.config with a
# stand-in bsub that only captures each job script: nothing is executed.
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
command -v nextflow >/dev/null || { echo "nextflow is needed for these tests"; exit 1; }

cat > "$TMP/main.nf" <<'EOF'
process HIGH {
    label 'process_high'
    script:
    "true"
}
process LOW {
    label 'process_low'
    script:
    "true"
}
process SHORTJOB {
    cpus 1
    memory 2.GB
    time 20.min
    script:
    "true"
}
process BIG {
    cpus 4
    memory 600.GB
    time 2.h
    script:
    "true"
}
workflow {
    HIGH()
    LOW()
    SHORTJOB()
    BIG()
}
EOF
mkdir -p "$TMP/bin"
printf '#!/bin/bash\nexit 0\n' > "$TMP/bin/bjobs"

# capture <name> [env...]: run the pipeline until all 4 jobs are submitted;
# leaves the captured scripts in $TMP/<name>/.
capture() {
  local cap="$TMP/$1" pid i
  shift
  mkdir -p "$cap"
  # shellcheck disable=SC2016  # $(...) and $((...)) belong to the stand-in
  printf '#!/bin/bash\nn=$(ls %s | wc -l); cat > %s/job.$n.sh; echo "Job <$((100+n))> is submitted to queue <q>."\n' "$cap" "$cap" > "$TMP/bin/bsub"
  chmod +x "$TMP/bin/bsub" "$TMP/bin/bjobs"
  (cd "$TMP" && exec env -u LSF_ENVDIR "$@" PATH="$TMP/bin:$PATH" NXF_OFFLINE=true \
     nextflow run main.nf -c "$REPO/conf/base.config" -c "$REPO/conf/stjude.config" > "$cap.log" 2>&1) &
  pid=$!
  for ((i = 0; i < 180; i++)); do
    [ "$(find "$cap" -name 'job.*.sh' | wc -l)" -ge 4 ] && break
    sleep 1
  done
  kill "$pid" 2>/dev/null || true
  wait "$pid" 2>/dev/null || true
}
directives() { # directives <capture dir> <job name>: the #BSUB lines of that job
  grep -l "^#BSUB -J $2\$" "$1"/job.*.sh | head -n 1 | xargs grep '^#BSUB' | grep -v '^#BSUB -o '
}
expect_directives() { # expect_directives <capture dir> <job> <expected lines> <name>
  local got
  got=$(directives "$1" "$2" 2>/dev/null || true)
  if [ "$got" = "$3" ]; then ok "$4"; else bad "$4" "got: $(printf '%q' "$got")"; fi
}

echo "St. Jude profile: LSF directives"

HIGH_EXPECTED='#BSUB -q standard
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00
#BSUB -M 6144MB
#BSUB -R "select[mem>=49152MB] rusage[mem=6144MB]"
#BSUB -J HIGH'

capture plain
expect_directives "$TMP/plain" HIGH "$HIGH_EXPECTED" "8 CPUs, 48 GB: rusage 6144 MB per slot (x 8 slots = 48 GB), one host, standard queue"
expect_directives "$TMP/plain" LOW '#BSUB -q standard
#BSUB -W 01:00
#BSUB -M 4096MB
#BSUB -R "select[mem>=4096MB] rusage[mem=4096MB]"
#BSUB -J LOW' "1 CPU, 4 GB: the whole request on the one slot"
expect_directives "$TMP/plain" SHORTJOB '#BSUB -q short
#BSUB -W 00:20
#BSUB -M 2048MB
#BSUB -R "select[mem>=2048MB] rusage[mem=2048MB]"
#BSUB -J SHORTJOB' "20 minutes: short queue (30-minute hard limit)"
expect_directives "$TMP/plain" BIG '#BSUB -q standard
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -W 02:00
#BSUB -M 153600MB
#BSUB -R "select[mem>=614400MB] rusage[mem=153600MB]"
#BSUB -J BIG' "600 GB: standard, not large_mem (HPCF recommends large_mem above 1 TB)"

# The cluster's lsf.conf must not change the request. Nextflow 26.10.4
# ignores lsf.conf for these two settings (unset means false); the profile
# sets both, so a Nextflow that honours lsf.conf gives the same directives.
mkdir -p "$TMP/lsfenv"
printf 'LSB_JOB_MEMLIMIT=Y\nRESOURCE_RESERVE_PER_TASK=N\nLSF_UNIT_FOR_LIMITS=KB\n' > "$TMP/lsfenv/lsf.conf"
capture lsfconf LSF_ENVDIR="$TMP/lsfenv"
expect_directives "$TMP/lsfconf" HIGH "$HIGH_EXPECTED" "an lsf.conf saying per-job limits and per-job reservation does not change them"

finish
