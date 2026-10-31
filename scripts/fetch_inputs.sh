#!/usr/bin/env bash
# Download and verify the inputs for a real run. Run on a login node (or a
# data-transfer node if HPCF provides one); compute nodes may not have
# internet access.
#
# Usage:
#   scripts/fetch_inputs.sh geo        <dest_dir>
#   scripts/fetch_inputs.sh reference  <dest_dir> <GRCh38-3.0.0|GRCh38-2020-A>
#   scripts/fetch_inputs.sh whitelists <dest_dir>
#   scripts/fetch_inputs.sh bams       <dest_dir> <runs.txt>     (DRY_RUN=1 to list without downloading)
#
# Every download is checked before it is accepted: gzip integrity and
# whitelist line counts; for BAMs, size and MD5 against ENA's record,
# `samtools quickcheck` (catches truncation) and the 10x header lines
# bamtofastq needs. Tool containers are run with Singularity (Docker is not
# available on the St. Jude cluster): `module load singularity/4.1.1` first.
set -euo pipefail

SRA_IMAGE=docker://quay.io/biocontainers/sra-tools:3.4.1--h4304569_1
SAMTOOLS_IMAGE=docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1
ENA_API=https://www.ebi.ac.uk/ena/portal/api/filereport

die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found on PATH${2:+ ($2)}"; }

fetch() {  # fetch <url> <dest>; resumes partial downloads
  local url=$1 dest=$2
  if [ -s "$dest" ]; then echo "exists: $dest"; return; fi
  echo "downloading $url"
  curl -fL --retry 3 -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}

samtools_() { singularity exec "$SAMTOOLS_IMAGE" samtools "$@" </dev/null; }

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
    # Sample labels from the column NAMES of the counts file (its first line),
    # with the number of cells per label. No count values are read.
    { gzip -cd "$dest/GSE155446_human_raw_counts.csv.gz" || true; } | head -n 1 \
      | tr ',' '\n' | tail -n +2 | sed -E 's/^"//; s/"$//; s/_[ACGTN]{16}(-[0-9]+)?$//' \
      | sort | uniq -c | awk '{ print $2 "\t" $1 }' > "$dest/published_labels.tsv"
    [ -s "$dest/published_labels.tsv" ] || die "could not read cell IDs from the counts header"
    echo "GEO files OK. Labels and cells per label (from column names only): $dest/published_labels.tsv"
    ;;

  reference)
    ver=${3:-}
    case "$ver" in
      GRCh38-3.0.0)  url=https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz ;;
      GRCh38-2020-A) url=https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz ;;
      *) die "reference version must be GRCh38-3.0.0 or GRCh38-2020-A (the one named in the BAM headers; see 'bams')" ;;
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
    # 10x barcode whitelists, from the STARsolo distribution (the copies in
    # Cell Ranger's GitHub repository are no longer published). Line counts
    # are checked so a wrong or truncated file cannot slip through.
    wl=https://labshare.cshl.edu/shares/dobin/dobin/STARsolo/10X/whitelists
    fetch "$wl/737K-august-2016.txt.gz" "$dest/737K-august-2016.txt.gz"
    fetch "$wl/3M-february-2018.txt.gz" "$dest/3M-february-2018.txt.gz"
    n2=$(gzip -cd "$dest/737K-august-2016.txt.gz" | wc -l)
    n3=$(gzip -cd "$dest/3M-february-2018.txt.gz" | wc -l)
    [ "$n2" -eq 737280 ] || die "v2 whitelist has $n2 lines, expected 737280"
    [ "$n3" -eq 6794880 ] || die "v3 whitelist has $n3 lines, expected 6794880"
    echo "whitelists OK (v2: $n2, v3: $n3 barcodes)"
    ;;

  bams)
    # The original Cell Ranger BAMs, as submitted, from ENA (which mirrors
    # SRA and records each file's size and MD5). Falls back to SRA's
    # `prefetch --type TenX` if ENA lists no BAM for a run; that path has no
    # MD5 to check against, so it relies on samtools quickcheck alone.
    runs=${3:-}
    [ -s "$runs" ] || die "give a file with one SRR accession per line"
    need curl
    if [ "${DRY_RUN:-0}" != 1 ]; then
      need md5sum
      need singularity "module load singularity/4.1.1"
    fi
    summary="$dest/bams_summary.tsv"
    [ "${DRY_RUN:-0}" = 1 ] || printf 'srr\tbam\tread_group_sample\treference\tstar_version\n' > "$summary"
    # Read the run list on fd 3 so no command in the loop can consume it.
    while read -r srr _ <&3; do
      case "$srr" in ''|\#*) continue ;; esac
      [[ "$srr" =~ ^SRR[0-9]+$ ]] || die "not an SRR accession: '$srr'"
      echo "== $srr"
      rec=$(curl -fsS --retry 3 "$ENA_API?accession=$srr&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes&format=tsv" | tail -n +2 || true)
      IFS=$'\t' read -r _ ftps md5s sizes <<< "$rec" || true
      IFS=';' read -r -a ftp_a <<< "${ftps:-}"
      IFS=';' read -r -a md5_a <<< "${md5s:-}"
      IFS=';' read -r -a size_a <<< "${sizes:-}"
      url='' md5='' size=''
      for i in "${!ftp_a[@]}"; do
        case "${ftp_a[$i]}" in
          *.bam) [ -z "$url" ] || die "$srr has more than one submitted BAM; pick one by hand"
                 url="https://${ftp_a[$i]}"; md5=${md5_a[$i]:-}; size=${size_a[$i]:-} ;;
        esac
      done
      if [ -n "$url" ]; then
        name=${url##*/}
        if [ "${DRY_RUN:-0}" = 1 ]; then
          awk -v n="$name" -v s="$size" -v m="$md5" 'BEGIN { printf "  %s  %.1f GB  md5 %s\n", n, s / 1e9, m }'
          continue
        fi
        mkdir -p "$dest/$srr"
        bam="$dest/$srr/$name"
        fetch "$url" "$bam"
        [ -z "$size" ] || [ "$(stat -c %s "$bam")" = "$size" ] || die "$bam is $(stat -c %s "$bam") bytes, ENA says $size (delete it and re-run)"
        if [ -n "$md5" ]; then
          echo "checking MD5"
          [ "$(md5sum "$bam" | cut -d' ' -f1)" = "$md5" ] || die "$bam MD5 does not match ENA's record (delete it and re-run)"
        fi
      else
        [ "${DRY_RUN:-0}" = 1 ] && { echo "  ENA lists no BAM; would fall back to prefetch --type TenX"; continue; }
        echo "ENA lists no BAM for $srr; trying SRA prefetch --type TenX (no MD5 available)"
        singularity exec "$SRA_IMAGE" prefetch --type TenX --max-size 200G -O "$dest" "$srr" </dev/null
        bam=$(find "$dest/$srr" -type f -name '*.bam*' ! -name '*.bai' ! -name '*.part' | head -n 1)
        [ -n "$bam" ] || die "no BAM downloaded for $srr"
        case "$bam" in *.bam) ;; *) mv "$bam" "${bam%.bam*}.bam"; bam="${bam%.bam*}.bam" ;; esac
      fi
      samtools_ quickcheck -u "$bam" || die "$bam failed samtools quickcheck (truncated or corrupt)"
      hdr=$(samtools_ view -H "$bam")
      grep -q '10x_bam_to_fastq' <<< "$hdr" \
        || die "$bam lacks the 10x_bam_to_fastq header lines; bamtofastq cannot convert it"
      # Informational fields; absent ones print as '?' rather than stopping the script.
      sm=$( { grep '^@RG' <<< "$hdr" | grep -o 'SM:[^[:space:]]*' | cut -d: -f2 | sort -u | paste -sd, -; } || true)
      ref=$( { grep -o 'refdata-[A-Za-z0-9._-]*' <<< "$hdr" | sort -u | paste -sd, -; } || true)
      star=$( { grep -o 'STAR_[0-9][0-9.a-z]*' <<< "$hdr" | sort -u | paste -sd, -; } || true)
      printf '%s\t%s\t%s\t%s\t%s\n' "$srr" "$bam" "${sm:-?}" "${ref:-?}" "${star:-?}" >> "$summary"
      echo "$srr OK: $bam (read group sample ${sm:-?}; reference ${ref:-?}; ${star:-STAR version not found})"
    done 3< "$runs"
    [ "${DRY_RUN:-0}" = 1 ] && exit 0
    echo "BAMs OK. Summary: $summary"
    echo "Map each BAM to its published label by the read-group sample / file name, not by the SRA run title."
    ;;

  *) die "unknown command '$cmd' (geo | reference | whitelists | bams)" ;;
esac
