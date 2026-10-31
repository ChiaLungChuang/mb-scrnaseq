#!/usr/bin/env bash
# Tests for bin/check_reads.sh: read pairing and chemistry detection on
# hand-built FASTQ folders laid out the way bamtofastq writes them.
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
CR="$REPO/bin/check_reads.sh"

fq() { # fq <path> <read length> [records]: a gzipped FASTQ of equal-length reads
  local path=$1 len=$2 n=${3:-20} seq qual i
  mkdir -p "$(dirname "$path")"
  seq=$(printf '%*s' "$len" '' | tr ' ' 'A')
  qual=$(printf '%*s' "$len" '' | tr ' ' 'F')
  for ((i = 1; i <= n; i++)); do printf '@r%d\n%s\n+\n%s\n' "$i" "$seq" "$qual"; done | gzip > "$path"
}
lane() { # lane <dir> <lane> <chunk> <r1 len> <cdna len> <cdna tag>
  fq "$1/bamtofastq_S1_L00$2_R1_00$3.fastq.gz" "$4"
  fq "$1/bamtofastq_S1_L00$2_$6_00$3.fastq.gz" "$5"
  if [ "$6" = R3 ]; then fq "$1/bamtofastq_S1_L00$2_R2_00$3.fastq.gz" 8; fi   # index read in the R3 layout
}

echo "check_reads.sh"

d=$TMP/v3_r2/fastq/sampleA_0_1_HFLOW; lane "$d" 1 1 28 90 R2
run "$CR" "$TMP/v3_r2/fastq" v3 "$TMP/out.tsv"
expect_rc 0 "v3, cDNA in R2: accepted"
expect_out "cDNA read is R2; 1 file pair(s)" "v3, cDNA in R2: reports layout"
expect_file_eq "$TMP/out.tsv" "$(printf '%s\t%s' "$d/bamtofastq_S1_L001_R1_001.fastq.gz" "$d/bamtofastq_S1_L001_R2_001.fastq.gz")" "v3, cDNA in R2: pairs R1 with R2"

d=$TMP/v3_r3/fastq/sampleB_0_1_HFLOW; lane "$d" 1 1 28 90 R3
run "$CR" "$TMP/v3_r3/fastq" v3 "$TMP/out.tsv"
expect_rc 0 "cDNA in R3 (seqnames layout): accepted"
expect_file_eq "$TMP/out.tsv" "$(printf '%s\t%s' "$d/bamtofastq_S1_L001_R1_001.fastq.gz" "$d/bamtofastq_S1_L001_R3_001.fastq.gz")" "cDNA in R3: pairs R1 with R3, not the R2 index read"

d=$TMP/v2/fastq/x; lane "$d" 1 1 26 98 R2
run "$CR" "$TMP/v2/fastq" v2 "$TMP/out.tsv"
expect_rc 0 "v2 (R1 = 26 nt): accepted"
expect_out "chemistry v2 confirmed (R1 = 26 nt)" "v2: reports chemistry"

d=$TMP/multi/fastq/x
lane "$d" 1 1 28 90 R2; lane "$d" 1 2 28 90 R2; lane "$d" 2 1 28 90 R2
run "$CR" "$TMP/multi/fastq" v3 "$TMP/out.tsv"
expect_out "3 file pair(s)" "two lanes, two chunks: three pairs"
bad_pairs=$(awk -F'\t' '{ a = $1; sub(/_R1_/, "_R2_", a); if (a != $2) n++ } END { print n + 0 }' "$TMP/out.tsv")
if [ "$bad_pairs" -eq 0 ]; then ok "every R1 paired with its own lane and chunk"; else bad "every R1 paired with its own lane and chunk" "$(cat "$TMP/out.tsv")"; fi

d=$TMP/rgname/fastq/lib_R1_x_0_1_HFLOW; lane "$d" 1 1 28 90 R2
run "$CR" "$TMP/rgname/fastq" v3 "$TMP/out.tsv"
expect_rc 0 "folder name containing _R1_: accepted"
expect_file_eq "$TMP/out.tsv" "$(printf '%s\t%s' "$d/bamtofastq_S1_L001_R1_001.fastq.gz" "$d/bamtofastq_S1_L001_R2_001.fastq.gz")" "folder name containing _R1_: only the file name is rewritten"

ln -s "$TMP/v3_r2/fastq" "$TMP/linked"
run "$CR" "$TMP/linked" v3 "$TMP/out.tsv"
expect_rc 0 "input folder given as a symlink (as Nextflow stages it): accepted"

run "$CR" "$TMP/v3_r2/fastq" v2 "$TMP/out.tsv"
expect_fail "chemistry v2 needs 26" "v3 reads declared v2: refused"
run "$CR" "$TMP/v2/fastq" v3 "$TMP/out.tsv"
expect_fail "chemistry v3 needs 28" "v2 reads declared v3: refused"

d=$TMP/mixed/fastq/x; fq "$d/bamtofastq_S1_L001_R1_001.fastq.gz" 28 10; fq "$d/bamtofastq_S1_L001_R2_001.fastq.gz" 90
{ gzip -cd "$d/bamtofastq_S1_L001_R1_001.fastq.gz"; printf '@x\n%s\n+\n%s\n' AAAAAAAAAAAAAAAAAAAAAAAAAA FFFFFFFFFFFFFFFFFFFFFFFFFF; } | gzip > "$TMP/m.gz"
mv "$TMP/m.gz" "$d/bamtofastq_S1_L001_R1_001.fastq.gz"
run "$CR" "$TMP/mixed/fastq" v3 "$TMP/out.tsv"
expect_fail "barcode-read length(s) 26 28" "mixed barcode-read lengths: refused"

d=$TMP/nopartner/fastq/x; fq "$d/bamtofastq_S1_L001_R1_001.fastq.gz" 28
run "$CR" "$TMP/nopartner/fastq" v3 "$TMP/out.tsv"
expect_fail "no R2 partner" "R1 without a cDNA file: refused"

d=$TMP/swapped/fastq/x; lane "$d" 1 1 28 28 R2
run "$CR" "$TMP/swapped/fastq" v3 "$TMP/out.tsv"
expect_fail "files look swapped" "cDNA read no longer than the barcode read: refused"

d=$TMP/emptyfq/fastq/x; mkdir -p "$d"; printf '' | gzip > "$d/bamtofastq_S1_L001_R1_001.fastq.gz"; fq "$d/bamtofastq_S1_L001_R2_001.fastq.gz" 90
run "$CR" "$TMP/emptyfq/fastq" v3 "$TMP/out.tsv"
expect_fail "could not read sequences" "empty R1 file: refused"

mkdir -p "$TMP/none/fastq"
run "$CR" "$TMP/none/fastq" v3 "$TMP/out.tsv"
expect_fail "no *_R1_*.fastq.gz files" "folder without FASTQs: refused"

run "$CR" "$TMP/v3_r2/fastq" v4 "$TMP/out.tsv"
expect_rc 2 "unknown chemistry: usage error (exit 2)"

finish
