#!/usr/bin/env bash
# Download and verify the inputs for a real run. Run on a login node (or a
# data-transfer node if HPCF provides one); compute nodes may not have
# internet access.
#
# Usage:
#   scripts/fetch_inputs.sh geo        <dest_dir>
#   scripts/fetch_inputs.sh reference  <dest_dir> <GRCh38-2020-A|GRCh38-3.0.0>
#   scripts/fetch_inputs.sh whitelists <dest_dir>
#   scripts/fetch_inputs.sh bams       <dest_dir> <runs.txt>
#
# Every download is checked before it is accepted: whitelist line counts,
# gzip integrity, and the 10x header lines bamtofastq needs in each BAM.
# Tool containers are run with Singularity (Docker is not available on the
# St. Jude cluster): `module load singularity/4.1.1` first.
set -euo pipefail

SRA_IMAGE=docker://quay.io/biocontainers/sra-tools:3.4.1--h4304569_1
SAMTOOLS_IMAGE=docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1

die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found on PATH${2:+ ($2)}"; }

fetch() {  # fetch <url> <dest>; resumes partial downloads
  local url=$1 dest=$2
  if [ -s "$dest" ]; then echo "exists: $dest"; return; fi
  echo "downloading $url"
  curl -fL --retry 3 -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}

cmd=${1:-}
dest=${2:-}
[ -n "$cmd" ] && [ -n "$dest" ] || { sed -n '6,10p' "$0"; exit 2; }
mkdir -p "$dest"

case "$cmd" in
  geo)
    # Supplementary files of GSE155446: the published count matrix and cell metadata.
    base=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155446/suppl
    for f in GSE155446_human_raw_counts.csv.gz GSE155446_human_cell_metadata.csv.gz; do
      fetch "$base/$f" "$dest/$f"
      gzip -t "$dest/$f" || die "$f is not a valid gzip file"
    done
    echo "GEO files OK. Look at column NAMES only (not values) to fill published_label:"
    printf '  %s\n' "zcat $dest/GSE155446_human_raw_counts.csv.gz | head -n 1 | tr ',' '\\n' | sed -n '2,6p'"
    ;;

  reference)
    ver=${3:-}
    case "$ver" in
      GRCh38-2020-A) url=https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz ;;
      GRCh38-3.0.0)  url=https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz ;;
      *) die "reference version must be GRCh38-2020-A or GRCh38-3.0.0 (match the Cell Ranger reference named in GEO's data-processing notes)" ;;
    esac
    tarball="$dest/$(basename "$url")"
    fetch "$url" "$tarball"
    tar -xzf "$tarball" -C "$dest"
    dir=$(find "$dest" -maxdepth 1 -type d -name "refdata-*${ver#GRCh38-}*" | head -n 1)
    [ -s "$dir/fasta/genome.fa" ] && [ -s "$dir/genes/genes.gtf" ] || die "unexpected layout in $tarball"
    echo "reference OK:"
    echo "  genome_fasta: $dir/fasta/genome.fa"
    echo "  genome_gtf:   $dir/genes/genes.gtf"
    ;;

  whitelists)
    # 10x barcode whitelists, from the Cell Ranger source on GitHub. Line counts
    # are checked so a wrong or truncated file cannot slip through.
    gh=https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes
    fetch "$gh/737K-august-2016.txt" "$dest/737K-august-2016.txt"
    fetch "$gh/3M-february-2018.txt.gz" "$dest/3M-february-2018.txt.gz"
    n2=$(wc -l < "$dest/737K-august-2016.txt")
    n3=$(gzip -cd "$dest/3M-february-2018.txt.gz" | wc -l)
    [ "$n2" -eq 737280 ] || die "v2 whitelist has $n2 lines, expected 737280"
    [ "$n3" -eq 6794880 ] || die "v3 whitelist has $n3 lines, expected 6794880"
    echo "whitelists OK (v2: $n2, v3: $n3 barcodes)"
    ;;

  bams)
    runs=${3:-}
    [ -s "$runs" ] || die "give a file with one SRR accession per line"
    need singularity "module load singularity/4.1.1"
    while read -r srr _; do
      case "$srr" in ''|\#*) continue ;; esac
      [[ "$srr" =~ ^SRR[0-9]+$ ]] || die "not an SRR accession: '$srr'"
      echo "== $srr"
      # --type TenX fetches the original Cell Ranger BAM, not SRA's re-encoded reads.
      singularity exec "$SRA_IMAGE" prefetch --type TenX --max-size 200G -O "$dest" "$srr"
      bam=$(find "$dest/$srr" -name '*.bam*' -type f | head -n 1)
      [ -n "$bam" ] || die "no BAM downloaded for $srr"
      case "$bam" in *.bam) ;; *) mv "$bam" "${bam%.bam*}.bam"; bam="${bam%.bam*}.bam" ;; esac
      singularity exec "$SAMTOOLS_IMAGE" samtools view -H "$bam" | grep -q '10x_bam_to_fastq' \
        || die "$bam lacks the 10x_bam_to_fastq header lines; bamtofastq cannot convert it"
      echo "$srr OK: $bam"
    done < "$runs"
    echo "BAMs OK. Put their paths in the samplesheet's bam column."
    ;;

  *) die "unknown command '$cmd' (geo | reference | whitelists | bams)" ;;
esac
