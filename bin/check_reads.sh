#!/usr/bin/env bash
# Pair up bamtofastq output and verify the declared 10x chemistry.
#
# Usage: check_reads.sh <fastq_dir> <declared_chemistry v2|v3> <out_tsv>
#
# bamtofastq names the cDNA read R2 or R3 depending on the BAM header (with a
# "10x_bam_to_fastq_seqnames" line the cDNA read is R3 and R2 is an index
# read). This script works out which, pairs every R1 (barcode + UMI) file
# with its cDNA file, and checks the R1 length against the declared
# chemistry: 26 nt for v2 (16 CB + 10 UMI), 28 nt for v3 (16 + 12).
# Any mismatch stops the run: quantifying with the wrong UMI length gives
# plausible-looking but wrong counts.
#
# Writes <out_tsv> with columns: r1<TAB>cdna (one row per lane chunk).
# find -L: Nextflow stages the input folder as a symlink.
set -euo pipefail

dir=${1:?fastq dir}
declared=${2:?declared chemistry}
out=${3:?output tsv}

case "$declared" in
  v2) expected=26 ;;
  v3) expected=28 ;;
  *) echo "ERROR: chemistry must be v2 or v3, got '$declared'" >&2; exit 2 ;;
esac

r1_files=$(find -L "$dir" -name '*_R1_*.fastq.gz' | sort)
if [ -z "$r1_files" ]; then
  echo "ERROR: no *_R1_*.fastq.gz files under $dir" >&2
  exit 1
fi

if [ -n "$(find -L "$dir" -name '*_R3_*.fastq.gz' | head -n 1)" ]; then
  cdna_tag=R3
else
  cdna_tag=R2
fi

: > "$out"
while IFS= read -r r1; do
  # Substitute in the file name only: the folder name comes from the BAM's
  # read group and may itself contain "_R1_".
  r1_name=${r1##*/}
  cdna="${r1%/*}/${r1_name/_R1_/_${cdna_tag}_}"
  if [ ! -s "$cdna" ]; then
    echo "ERROR: no ${cdna_tag} partner for $r1" >&2
    exit 1
  fi
  # Read lengths from the first 1000 records; '|| true' tolerates the
  # SIGPIPE gzip gets when head stops early. An empty result is caught below.
  r1_lens=$( { gzip -cd "$r1" 2>/dev/null || true; } | head -n 4000 | awk 'NR % 4 == 2 { print length($0) }' | sort -u)
  cdna_len=$( { gzip -cd "$cdna" 2>/dev/null || true; } | head -n 4000 | awk 'NR % 4 == 2 { print length($0) }' | sort -n | tail -n 1)
  if [ -z "$r1_lens" ] || [ -z "$cdna_len" ]; then
    echo "ERROR: could not read sequences from $r1 or $cdna" >&2
    exit 1
  fi
  if [ "$(printf '%s\n' "$r1_lens" | wc -l)" -ne 1 ] || [ "$r1_lens" -ne "$expected" ]; then
    echo "ERROR: $r1 has barcode-read length(s) $(tr '\n' ' ' <<< "$r1_lens")but chemistry $declared needs $expected." >&2
    echo "       Check the chemistry column of the samplesheet for this sample." >&2
    exit 1
  fi
  if [ "$cdna_len" -le "$expected" ]; then
    echo "ERROR: cDNA read $cdna (max length $cdna_len) is not longer than the barcode read; files look swapped." >&2
    exit 1
  fi
  printf '%s\t%s\n' "$r1" "$cdna" >> "$out"
done <<< "$r1_files"

echo "chemistry $declared confirmed (R1 = $expected nt); cDNA read is $cdna_tag; $(wc -l < "$out") file pair(s)"
