#!/usr/bin/env bash
# Tests for scripts/fetch_inputs.sh. Network and containers are replaced by
# stand-ins on PATH: `curl` serves local files (and a fake ENA record), and
# `singularity exec IMAGE cmd...` runs cmd directly (samtools must be on
# PATH; `prefetch` is simulated). So these test the script's checks and
# bookkeeping, not the real servers.
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
FI="$REPO/scripts/fetch_inputs.sh"
command -v samtools >/dev/null || { echo "samtools is needed for these tests"; exit 1; }

# ---- stand-ins ---------------------------------------------------------------
mkdir -p "$TMP/bin" "$TMP/served"
cat > "$TMP/bin/curl" <<'EOF'
#!/usr/bin/env bash
# Fake curl: ENA API queries print $MOCK_ENA (a file); any other URL is served
# from $MOCK_SERVED/<file name>. Every URL is logged to $MOCK_LOG.
out='' url=''
while [ $# -gt 0 ]; do
  case "$1" in -o) out=$2; shift 2 ;; -C|--retry) shift 2 ;; -*) shift ;; *) url=$1; shift ;; esac
done
echo "$url" >> "$MOCK_LOG"
if [[ "$url" == *ena/portal/api* ]]; then
  [ -z "${MOCK_ENA_FAIL:-}" ] || exit 22
  cat "$MOCK_ENA"
else
  src="$MOCK_SERVED/${url##*/}"
  [ -f "$src" ] || exit 22
  if [ -n "$out" ]; then cp "$src" "$out"; else cat "$src"; fi
fi
EOF
cat > "$TMP/bin/singularity" <<'EOF'
#!/usr/bin/env bash
# Fake singularity: "exec IMAGE cmd..." runs cmd; prefetch is simulated the
# way SRA names TenX downloads (<SRR>/<name>.bam.1, plus an index).
shift 2
if [ "$1" = prefetch ]; then
  dest='' srr=''
  while [ $# -gt 0 ]; do case "$1" in -O) dest=$2; shift 2 ;; --type|--max-size) shift 2 ;; prefetch) shift ;; *) srr=$1; shift ;; esac; done
  mkdir -p "$dest/$srr"
  cp "$MOCK_PREFETCH_BAM" "$dest/$srr/possorted_genome_bam.bam.1"
  echo index > "$dest/$srr/possorted_genome_bam.bam.bai"
  exit 0
fi
exec "$@"
EOF
chmod +x "$TMP/bin/curl" "$TMP/bin/singularity"
export PATH="$TMP/bin:$PATH" MOCK_SERVED="$TMP/served" MOCK_LOG="$TMP/curl.log" MOCK_ENA="$TMP/ena.tsv"

# ---- fixtures ----------------------------------------------------------------
src_bam="$REPO/tests/data/sampleB.bam"
samtools view -H "$src_bam" > "$TMP/h.sam"
sed -i 's/SM:[^[:space:]]*/SM:1_943/' "$TMP/h.sam"
printf '@PG\tID:STAR\tPN:STAR\tVN:STAR_2.5.1b\n@CO\tuser command line: STAR --genomeDir /x/refdata-cellranger-GRCh38-3.0.0/star\n' >> "$TMP/h.sam"
samtools reheader "$TMP/h.sam" "$src_bam" > "$TMP/served/1_943_possorted_genome_bam.bam"
good="$TMP/served/1_943_possorted_genome_bam.bam"
grep -v '10x_bam_to_fastq' "$TMP/h.sam" > "$TMP/h2.sam"
samtools reheader "$TMP/h2.sam" "$src_bam" > "$TMP/served/no10x_possorted_genome_bam.bam"
head -c 30000 "$good" > "$TMP/served/trunc_possorted_genome_bam.bam"

md5() { md5sum "$1" | cut -d' ' -f1; }
size() { stat -c %s "$1"; }
ena() { # ena <file> [<md5> [<size>]]: fake ENA record for a BAM plus its index
  local f=$1 m=${2:-$(md5 "$1")} s=${3:-$(size "$1")}
  printf 'run_accession\tsubmitted_ftp\tsubmitted_md5\tsubmitted_bytes\n' > "$MOCK_ENA"
  printf 'SRR1\tftp.example/vol1/%s;ftp.example/vol1/%s.bai\t%s;0\t%s;1\n' "${f##*/}" "${f##*/}" "$m" "$s" >> "$MOCK_ENA"
}
printf '# runs\n\nSRR1\n' > "$TMP/runs.txt"
fresh() { rm -rf "$TMP/out"; : > "$MOCK_LOG"; }

echo "fetch_inputs.sh bams"

fresh; ena "$good"; run bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_rc 0 "good BAM: accepted"
expect_out "read group sample 1_943" "good BAM: reports the read-group sample"
expect_file_eq "$TMP/out/bams_summary.tsv" "$(printf 'srr\tbam\tread_group_sample\treference\tstar_version\nSRR1\t%s\t1_943\trefdata-cellranger-GRCh38-3.0.0\tSTAR_2.5.1b' "$TMP/out/SRR1/1_943_possorted_genome_bam.bam")" "good BAM: summary row has sample, reference and STAR version"
if cmp -s "$good" "$TMP/out/SRR1/1_943_possorted_genome_bam.bam"; then ok "good BAM: downloaded file identical to the source"; else bad "good BAM: downloaded file identical to the source"; fi

: > "$MOCK_LOG"; run bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_out "exists:" "re-run: existing download is reused"
if ! grep -q 1_943 "$MOCK_LOG"; then ok "re-run: BAM not fetched again"; else bad "re-run: BAM not fetched again" "$(cat "$MOCK_LOG")"; fi

fresh; ena "$good" 00000000000000000000000000000000; run bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_fail "MD5 does not match" "wrong MD5: refused"

fresh; ena "$good" "" 123; run bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_fail "ENA says 123" "wrong size: refused"

fresh; ena "$TMP/served/trunc_possorted_genome_bam.bam"; run bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_fail "failed samtools quickcheck" "truncated BAM whose MD5 matches: caught by quickcheck"

fresh; ena "$TMP/served/no10x_possorted_genome_bam.bam"; run bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_fail "lacks the 10x_bam_to_fastq header lines" "BAM without the 10x header lines: refused"

fresh; ena "$good"; sed -i '2s#\t#\tftp.example/vol1/second.bam;#' "$MOCK_ENA"; run bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_fail "more than one submitted BAM" "two BAMs listed for one run: refused"

fresh; ena "$good"; run env DRY_RUN=1 bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_rc 0 "dry run: succeeds"
expect_out "1_943_possorted_genome_bam.bam" "dry run: lists the file"
if [ ! -e "$TMP/out/SRR1" ] && ! grep -q 1_943 "$MOCK_LOG"; then ok "dry run: downloads nothing"; else bad "dry run: downloads nothing" "$(cat "$MOCK_LOG")"; fi

fresh; export MOCK_PREFETCH_BAM="$good"; run env MOCK_ENA_FAIL=1 bash "$FI" bams "$TMP/out" "$TMP/runs.txt"
expect_rc 0 "ENA unavailable: falls back to SRA prefetch"
if [ -s "$TMP/out/SRR1/possorted_genome_bam.bam" ] && [ "$(cat "$TMP/out/SRR1/possorted_genome_bam.bam.bai")" = index ]; then
  ok "prefetch fallback: .bam.1 renamed to .bam, index left alone"
else bad "prefetch fallback: .bam.1 renamed to .bam, index left alone" "$(ls -la "$TMP/out/SRR1")"; fi

fresh; printf 'SRR1\nERR42\n' > "$TMP/badruns.txt"; ena "$good"; run bash "$FI" bams "$TMP/out" "$TMP/badruns.txt"
expect_fail "not an SRR accession: 'ERR42'" "non-SRR accession: refused"

echo "fetch_inputs.sh geo"
fresh
bc() { printf '%s' "$1" | head -c 16; }
{ printf '"",'; printf '"MB1_%s","MB1_%s-1",MB22_%s\n' "$(bc ACGTACGTACGTACGTAA)" "$(bc CCCCGGGGAAAATTTT)" "$(bc TTTTAAAACCCCGGGG)"
  printf 'GENE1,1,0,3\nGENE2,0,5,0\n'; } | gzip > "$TMP/served/GSE155446_human_raw_counts.csv.gz"
printf 'cell,label\nx,y\n' | gzip > "$TMP/served/GSE155446_human_cell_metadata.csv.gz"
run bash "$FI" geo "$TMP/out"
expect_rc 0 "GEO files: accepted"
expect_file_eq "$TMP/out/published_labels.tsv" "$(printf 'MB1\t2\nMB22\t1')" "labels and cells per label from the column names (quotes and -1 handled)"
fresh; printf 'not gzip' > "$TMP/served/GSE155446_human_cell_metadata.csv.gz"
run bash "$FI" geo "$TMP/out"
expect_fail "not a valid gzip file" "corrupt GEO download: refused"

echo "fetch_inputs.sh whitelists"
fresh; seq 1 737280 | gzip > "$TMP/served/737K-august-2016.txt.gz"; seq 1 6794880 | gzip > "$TMP/served/3M-february-2018.txt.gz"
run bash "$FI" whitelists "$TMP/out"
expect_out "whitelists OK (v2: 737280, v3: 6794880 barcodes)" "complete whitelists: accepted"
fresh; seq 1 737279 | gzip > "$TMP/served/737K-august-2016.txt.gz"
run bash "$FI" whitelists "$TMP/out"
expect_fail "expected 737280" "short v2 whitelist: refused"

echo "fetch_inputs.sh reference"
fresh; mkdir -p "$TMP/ref/refdata-cellranger-GRCh38-3.0.0/fasta" "$TMP/ref/refdata-cellranger-GRCh38-3.0.0/genes"
echo '>chr1' > "$TMP/ref/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"; echo 'gtf' > "$TMP/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
tar -czf "$TMP/served/refdata-cellranger-GRCh38-3.0.0.tar.gz" -C "$TMP/ref" refdata-cellranger-GRCh38-3.0.0
run bash "$FI" reference "$TMP/out" GRCh38-3.0.0
expect_out "genome_gtf:   $TMP/out/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf" "reference: reports the FASTA and GTF paths"
run bash "$FI" reference "$TMP/out" GRCh37
expect_fail "reference version must be" "unknown reference version: refused"

finish
