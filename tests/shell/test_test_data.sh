#!/usr/bin/env bash
# The committed synthetic datasets must be exactly what their generators
# produce, so a change to a generator cannot leave stale test data behind.
# BAMs are compared as SAM text without @PG lines (compression and @PG
# lines differ across samtools versions).
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
command -v samtools >/dev/null || { echo "samtools is needed for these tests"; exit 1; }

compare_dir() { # compare_dir <generator> <committed dir> <name>
  local gen=$1 committed=$2 name=$3 f diffs=''
  python3 "$gen" --outdir "$TMP/$name" > /dev/null
  for f in "$TMP/$name"/*; do
    f=${f##*/}
    if [[ "$f" == *.bam ]]; then
      cmp -s <(samtools view -h --no-PG "$TMP/$name/$f" | grep -v '^@PG') \
        <(samtools view -h --no-PG "$committed/$f" | grep -v '^@PG') || diffs+=" $f"
    else
      cmp -s "$TMP/$name/$f" "$committed/$f" || diffs+=" $f"
    fi
  done
  if [ -z "$diffs" ]; then ok "$name: committed files match the generator"
  else bad "$name: committed files match the generator" "differ:$diffs (re-run the generator and commit)"; fi
}

echo "synthetic test data"
compare_dir "$REPO/tests/data/make_test_data.py" "$REPO/tests/data" basic
compare_dir "$REPO/tests/data/make_accuracy_data.py" "$REPO/tests/data/accuracy" accuracy

finish
