#!/usr/bin/env bash
# Run the R code of PREP_PUBLISHED, SAMPLE_QC and CONCORDANCE on the
# full-size synthetic inputs from make_scale_data.R and report each step's
# peak memory (maximum resident set size) and wall time, as
# /usr/bin/time -v would. See tests/scale/README.md.
#
# Usage: tests/scale/rehearse_r_steps.sh DATA_DIR [THREADS]
#   DATA_DIR  the --out directory of make_scale_data.R
#   THREADS   data.table threads (default 1: process_r tasks get one CPU)
# Rscript is taken from PATH (RSCRIPT=... to override). Writes the steps'
# outputs and logs to DATA_DIR/rehearsal/.
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
data=${1:?usage: rehearse_r_steps.sh DATA_DIR [THREADS]}
threads=${2:-1}
rscript=${RSCRIPT:-Rscript}
out="$data/rehearsal"
mkdir -p "$out"
export OMP_THREAD_LIMIT=$threads R_DATATABLE_NUM_THREADS=$threads

# measure <name> <command...>: run it, append "name peak_MB wall_s cpu_s rc" to $out/measurements.tsv.
measure() {
  local name=$1; shift
  python3 - "$name" "$out" "$@" <<'EOF'
import resource, subprocess, sys, time
name, out, cmd = sys.argv[1], sys.argv[2], sys.argv[3:]
t = time.time()
with open(f"{out}/{name}.log", "w") as log:
    rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
wall = time.time() - t
ru = resource.getrusage(resource.RUSAGE_CHILDREN)   # ru_maxrss: largest child, KB on Linux
row = f"{name}\t{ru.ru_maxrss / 1024:.0f}\t{wall:.1f}\t{ru.ru_utime + ru.ru_stime:.1f}\t{rc}"
print(row)
with open(f"{out}/measurements.tsv", "a") as f:
    f.write(row + "\n")
sys.exit(rc)
EOF
}

printf 'step\tpeak_rss_MB\twall_s\tcpu_s\texit\n' > "$out/measurements.tsv"
echo "R: $("$rscript" -e 'cat(R.version.string, "Matrix", format(packageVersion("Matrix")), "data.table", format(packageVersion("data.table")))' 2>/dev/null); threads: $threads"
cd "$out"
measure PREP_PUBLISHED "$rscript" "$REPO/bin/prep_published.R" --lib "$REPO/lib/R" \
  --counts "$data/published_counts.csv.gz" --out-prefix published_counts --chunk-rows 2000
measure SAMPLE_QC "$rscript" "$REPO/bin/sample_qc.R" --lib "$REPO/lib/R" --sample S1 \
  --filtered-dir "$data/filtered" --nmads 3 --out-prefix S1
measure CONCORDANCE "$rscript" "$REPO/bin/concordance.R" --lib "$REPO/lib/R" --sample S1 \
  --label "${LABEL:-1167}" --raw-dir "$data/raw" --filtered-dir "$data/filtered" \
  --published published_counts.rds --criteria "$REPO/preregistration/acceptance_criteria.yaml" \
  --out-prefix S1
cat "$out/measurements.tsv"
