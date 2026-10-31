#!/usr/bin/env python3
"""Generate a tiny, fully synthetic 10x-like dataset for pipeline tests.

Everything here is invented: a two-contig genome, a handful of genes, and
reads for a few dozen cell barcodes per sample. The point is not biology;
it is that the ground truth is known exactly, so every pipeline step can be
checked against it.

Two samples are written because Cell Ranger BAMs come in two header layouts
that make bamtofastq name the cDNA read differently:
  sampleA  no "10x_bam_to_fastq_seqnames" header line -> cDNA read is R2
  sampleB  with "seqnames:R1,R3,I1,R2"                 -> cDNA read is R3
The pipeline has to handle both.

Outputs (written next to this script unless --outdir is given):
  genome.fa              two contigs: chrT1 (nuclear) and chrM (mito)
  genes.gtf              single-exon genes on the + strand
  whitelist_v3.txt       every barcode used, as a real whitelist would list
  sampleA.bam, sampleB.bam
                         Cell Ranger style BAMs (unmapped records carrying
                         CR/CY/UR/UY tags plus the header lines bamtofastq
                         reads)
  published_counts.csv   genes x cells UMI counts for the true cells of both
                         samples, laid out like a GEO supplementary matrix
                         (column names <label>_<barcode>, rows gene symbols)
  truth.json             per-cell expected UMI counts and which cells were
                         made deliberately low quality

Usage: python3 make_test_data.py [--outdir DIR] [--samtools PATH]
Standard library only, plus samtools for SAM -> BAM.
"""
import argparse
import json
import random
import shutil
import subprocess
from pathlib import Path

SEED = 20261030
CB_LEN, UMI_LEN = 16, 12          # 10x 3' v3 geometry
R2_LEN = 90
FLOWCELL = "HTESTFLOW"
SAMPLES = [
    # (sample id, label used in the "published" matrix, write seqnames line?)
    ("sampleA", "MB001", False),
    ("sampleB", "MB002", True),
]
ACGT = "ACGT"


def rand_seq(rng, n):
    return "".join(rng.choice(ACGT) for _ in range(n))


def build_genome(rng):
    contigs = {"chrT1": rand_seq(rng, 20000), "chrM": rand_seq(rng, 3000)}
    genes, pos = [], 500
    for i in range(1, 11):        # 10 nuclear genes
        length = rng.randint(600, 1200)
        genes.append(dict(chrom="chrT1", start=pos, end=pos + length - 1,
                          gene_id=f"ENSGTEST{i:05d}", gene_name=f"GENE{i}"))
        pos += length + rng.randint(300, 600)
    mpos = 200
    for j, name in enumerate(("MT-CO1", "MT-ND1"), start=1):
        genes.append(dict(chrom="chrM", start=mpos, end=mpos + 899,
                          gene_id=f"ENSGTESTM{j:04d}", gene_name=name))
        mpos += 1100
    return contigs, genes


def write_fasta(contigs, path):
    with open(path, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def write_gtf(genes, path):
    with open(path, "w") as fh:
        for g in genes:
            base = f'gene_id "{g["gene_id"]}"; gene_name "{g["gene_name"]}"; gene_type "protein_coding";'
            tx = base + f' transcript_id "{g["gene_id"].replace("ENSG", "ENST")}";'
            for feature, attrs in (("gene", base), ("transcript", tx), ("exon", tx)):
                fh.write("\t".join([g["chrom"], "synthetic", feature, str(g["start"]),
                                    str(g["end"]), ".", "+", ".", attrs]) + "\n")


def make_sample(rng, sample, barcodes, contigs, genes, with_seqnames, samtools, out):
    cells, ambient = barcodes[:20], barcodes[20:45]
    low_quality = set(cells[:3])   # high mito, low complexity: QC should flag these
    nuclear = [g for g in genes if g["chrom"] == "chrT1"]
    mito = [g for g in genes if g["chrom"] == "chrM"]

    molecules, truth = [], {}
    for cb in cells:
        if cb in low_quality:
            n_mol, pool, p_mito = rng.randint(60, 90), nuclear[:3], 0.6
        else:
            n_mol, pool, p_mito = rng.randint(250, 400), nuclear, 0.05
        umis, counts = set(), {}
        while len(umis) < n_mol:
            umi = rand_seq(rng, UMI_LEN)
            if umi in umis:
                continue
            umis.add(umi)
            gene = rng.choice(mito) if rng.random() < p_mito else rng.choice(pool)
            molecules.append((cb, umi, gene))
            counts[gene["gene_name"]] = counts.get(gene["gene_name"], 0) + 1
        truth[cb] = {"total_umi": n_mol, "low_quality": cb in low_quality, "counts": counts}
    for cb in ambient:             # empty droplets: a handful of ambient molecules
        for _ in range(rng.randint(2, 6)):
            molecules.append((cb, rand_seq(rng, UMI_LEN), rng.choice(nuclear)))

    rg_id = f"{sample}:0:1:{FLOWCELL}:1"
    header = [
        "@HD\tVN:1.6\tSO:unsorted",
        f"@RG\tID:{rg_id}\tSM:{sample}\tLB:0.1\tPU:{rg_id}\tPL:ILLUMINA",
        "@CO\t10x_bam_to_fastq:R1(CR:CY,UR:UY)",
        "@CO\t10x_bam_to_fastq:R2(SEQ:QUAL)",
    ]
    if with_seqnames:
        header.append("@CO\t10x_bam_to_fastq_seqnames:R1,R3,I1,R2")
    records, read_n = [], 0
    for cb, umi, gene in molecules:
        tx = contigs[gene["chrom"]][gene["start"] - 1:gene["end"]]
        for _ in range(rng.randint(1, 3)):          # PCR duplicates share a UMI
            read_n += 1
            start = rng.randint(max(0, len(tx) - 400), len(tx) - R2_LEN)   # 3' biased
            tags = [f"CR:Z:{cb}", f"CY:Z:{'F' * CB_LEN}", f"CB:Z:{cb}-1",
                    f"UR:Z:{umi}", f"UY:Z:{'F' * UMI_LEN}", f"RG:Z:{rg_id}"]
            records.append("\t".join([f"{FLOWCELL}:1:1101:{read_n}:{read_n}", "4", "*", "0", "0",
                                      "*", "*", "0", "0", tx[start:start + R2_LEN], "F" * R2_LEN] + tags))
    rng.shuffle(records)
    sam = out / f"{sample}.sam"
    sam.write_text("\n".join(header + records) + "\n")
    subprocess.run([samtools, "view", "--no-PG", "-b", "-o", str(out / f"{sample}.bam"), str(sam)], check=True)
    sam.unlink()
    return cells, truth, read_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(Path(__file__).parent))
    ap.add_argument("--samtools", default=shutil.which("samtools") or "samtools")
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    contigs, genes = build_genome(rng)
    write_fasta(contigs, out / "genome.fa")
    write_gtf(genes, out / "genes.gtf")

    whitelist = set()
    while len(whitelist) < 45 * len(SAMPLES):
        whitelist.add(rand_seq(rng, CB_LEN))
    whitelist = sorted(whitelist)
    (out / "whitelist_v3.txt").write_text("\n".join(whitelist) + "\n")

    published_cols, truth_all = [], {}
    for k, (sample, label, with_seqnames) in enumerate(SAMPLES):
        barcodes = whitelist[45 * k:45 * (k + 1)]
        cells, truth, n_reads = make_sample(rng, sample, barcodes, contigs, genes,
                                            with_seqnames, args.samtools, out)
        published_cols += [(label, cb, truth[cb]) for cb in cells]
        truth_all[sample] = {"published_label": label, "chemistry": "v3",
                             "cdna_read": "R3" if with_seqnames else "R2",
                             "n_reads": n_reads, "cells": truth}
        print(f"{sample}: {n_reads} reads, {len(cells)} cells, cDNA read {truth_all[sample]['cdna_read']}")

    with open(out / "published_counts.csv", "w") as fh:
        fh.write("," + ",".join(f"{label}_{cb}" for label, cb, _ in published_cols) + "\n")
        for g in genes:
            name = g["gene_name"]
            fh.write(name + "," + ",".join(str(t["counts"].get(name, 0)) for _, _, t in published_cols) + "\n")
    (out / "truth.json").write_text(json.dumps(truth_all, indent=1))


if __name__ == "__main__":
    main()
