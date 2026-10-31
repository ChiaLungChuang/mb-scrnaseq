#!/usr/bin/env python3
"""Generate the accuracy dataset: synthetic 10x reads whose correct counts
are known for every cell and gene under Cell Ranger's counting rules.

The basic test dataset (make_test_data.py) has only clean, unspliced,
uniquely mapping reads. Real data is not like that, and each of the
following is a place where a re-quantification can silently differ from
Cell Ranger. Every cell here gets molecules of each kind:

  kind          what the reads are                         counted as
  normal        sense, exonic, 1-4 PCR duplicates          1 UMI for the gene
  spliced       span the exon-exon junction of GSPL        1 UMI for GSPL
                (35-55 nt on the first exon, so they only
                map if the aligner splices them)
  minus_strand  sense reads of a gene on the - strand      1 UMI for GMINUS
  cb_error      every read's barcode has 1 substitution    1 UMI for the true cell
                (a single whitelist barcode is 1 away)     (barcode corrected)
  umi_error     3 reads with the UMI, 1 read with a        1 UMI, not 2
                1-substitution copy of it                  (UMI collapsed)
  multigene     one UMI: 3 reads on gene X, 1 on gene Y    1 UMI for X, 0 for Y
  antisense     reverse complement of a gene's sequence    0
  intronic      sense reads inside GSPL's intron           0
  intergenic    reads between genes                        0
  multimapper   reads inside a 500-nt block that two       0
                genes (GDUP1, GDUP2) share exactly

Two genes (TWIN on two loci) share a gene symbol, as ~20 human genes do;
the published matrix names them TWIN and TWIN.1, as Seurat's Read10X does.

Two samples, covering both chemistries and both bamtofastq layouts:
  accA  v3 (16 + 12 nt barcode read, 90-nt cDNA), cDNA read R2
  accB  v2 (16 + 10 nt barcode read, 98-nt cDNA), cDNA read R3

All whitelist barcodes are at least 3 substitutions apart, and all UMIs in
a cell at least 3 apart, so every correction above has exactly one right
answer.

Outputs (in --outdir, default: accuracy/ next to this script): genome.fa,
genes.gtf, whitelist_v2.txt, whitelist_v3.txt, accA.bam, accB.bam,
samplesheet.csv, published_counts.csv, truth.json (expected counts per
sample, barcode and gene ID, the true cells, and every molecule's kind).

Usage: python3 make_accuracy_data.py [--outdir DIR] [--samtools PATH]
"""
import argparse
import json
import random
import shutil
import subprocess
from pathlib import Path

SEED = 20261031
ACGT = "ACGT"
FLOWCELL = "HACCFLOW"
COMP = str.maketrans("ACGT", "TGCA")
SAMPLES = [
    # id, label in the published matrix, chemistry, UMI length, cDNA length, seqnames line (cDNA = R3)
    ("accA", "ACC1", "v3", 12, 90, False),
    ("accB", "ACC2", "v2", 10, 98, True),
]
N_CELLS, N_AMBIENT, N_DECOYS = 12, 8, 300
SPECIAL = {"cb_error": 3, "umi_error": 3, "multigene": 2, "antisense": 3,
           "intronic": 3, "intergenic": 2, "multimapper": 3}


def rand_seq(rng, n):
    return "".join(rng.choice(ACGT) for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def spaced_set(rng, n, length, min_dist, avoid=()):
    """n random sequences, pairwise (and to `avoid`) at least min_dist apart."""
    out = []
    while len(out) < n:
        s = rand_seq(rng, length)
        if len(set(s)) > 1 and all(hamming(s, t) >= min_dist for t in list(avoid) + out):
            out.append(s)
    return out


def substitute(rng, s):
    i = rng.randrange(len(s))
    return s[:i] + rng.choice([b for b in ACGT if b != s[i]]) + s[i + 1:]


def build_genome(rng):
    """Genes, their transcript sequences, and the regions the zero-count reads come from."""
    genes, seqs = [], {}
    layout, pos = [], 500
    singles = [f"G0{i}" for i in range(1, 7)] + ["TWIN", "TWIN"]
    for name in singles:
        length = rng.randint(700, 1000)
        layout.append(dict(name=name, strand="+", exons=[(pos, pos + length - 1)]))
        pos += length + rng.randint(400, 700)
    layout.append(dict(name="GSPL", strand="+", exons=[(pos, pos + 299), (pos + 1500, pos + 1899)]))
    intron = (pos + 300, pos + 1499)
    pos += 1900 + 500
    layout.append(dict(name="GMINUS", strand="-", exons=[(pos, pos + 799)]))
    pos += 800 + 400
    intergenic = (pos, pos + 1999)
    pos += 2000 + 500
    chrA1 = list(rand_seq(rng, pos))

    # chrA2: GDUP1 = unique 400 + shared 500; GDUP2 = shared 500 + unique 400.
    shared = rand_seq(rng, 500)
    chrA2 = list(rand_seq(rng, 300 + 900 + 800 + 900 + 300))
    chrA2[700:1200] = shared
    chrA2[2000:2500] = shared
    layout.append(dict(name="GDUP1", strand="+", exons=[(301, 1200)], chrom="chrA2"))
    layout.append(dict(name="GDUP2", strand="+", exons=[(2001, 2900)], chrom="chrA2"))
    chrM = list(rand_seq(rng, 3000))
    layout.append(dict(name="MT-CO1", strand="+", exons=[(200, 1099)], chrom="chrM"))
    layout.append(dict(name="MT-ND1", strand="+", exons=[(1300, 2199)], chrom="chrM"))
    contigs = {"chrA1": "".join(chrA1), "chrA2": "".join(chrA2), "chrM": "".join(chrM)}

    for k, g in enumerate(layout, start=1):
        g.setdefault("chrom", "chrA1")
        g["gene_id"] = f"ENSGACC{k:05d}"
        g["start"], g["end"] = g["exons"][0][0], g["exons"][-1][1]
        tx = "".join(contigs[g["chrom"]][s - 1:e] for s, e in g["exons"])
        seqs[g["gene_id"]] = revcomp(tx) if g["strand"] == "-" else tx
        genes.append(g)
    regions = {
        "intron": contigs["chrA1"][intron[0] - 1:intron[1]],
        "intergenic": contigs["chrA1"][intergenic[0] + 99:intergenic[1] - 100],
        "shared": shared,
    }
    return contigs, genes, seqs, regions


def write_fasta(contigs, path):
    with open(path, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def write_gtf(genes, path):
    with open(path, "w") as fh:
        for g in genes:
            base = f'gene_id "{g["gene_id"]}"; gene_name "{g["name"]}"; gene_type "protein_coding";'
            tx = base + f' transcript_id "{g["gene_id"].replace("ENSG", "ENST")}";'
            rows = [("gene", g["start"], g["end"], base), ("transcript", g["start"], g["end"], tx)]
            rows += [("exon", s, e, tx) for s, e in g["exons"]]
            for feature, s, e, attrs in rows:
                fh.write("\t".join([g["chrom"], "synthetic", feature, str(s), str(e), ".",
                                    g["strand"], ".", attrs]) + "\n")


def read_from(rng, genes_by_name, seqs, regions, kind, gene, rlen):
    """One cDNA read sequence for a molecule of the given kind."""
    if kind == "antisense":
        tx = seqs[gene["gene_id"]]
        s = rng.randint(len(tx) - 400, len(tx) - rlen)
        return revcomp(tx[s:s + rlen])
    if kind in ("intronic", "intergenic", "multimapper"):
        reg = regions[{"intronic": "intron", "intergenic": "intergenic", "multimapper": "shared"}[kind]]
        s = rng.randint(0, len(reg) - rlen)
        return reg[s:s + rlen]
    tx = seqs[gene["gene_id"]]
    name = gene["name"]
    if name == "GSPL":
        if kind == "spliced":
            s = 300 - rng.randint(35, 55)          # crosses the junction at transcript position 300
        else:
            s = rng.randint(300, len(tx) - rlen)   # within exon 2
    elif name == "GDUP1":
        s = rng.randint(0, 400 - rlen)             # unique part only
    elif name == "GDUP2":
        s = rng.randint(500, 900 - rlen)
    else:
        s = rng.randint(max(0, len(tx) - 400), len(tx) - rlen)   # 3' biased
    return tx[s:s + rlen]


def make_sample(rng, spec, barcodes, genes, seqs, regions, samtools, out):
    sample, label, chem, umi_len, rlen, with_seqnames = spec
    cells, ambient = barcodes[:N_CELLS], barcodes[N_CELLS:N_CELLS + N_AMBIENT]
    by_name = {}
    for g in genes:
        by_name.setdefault(g["name"], []).append(g)
    countable = [g for g in genes if g["name"] not in ("MT-CO1", "MT-ND1")]
    mito = by_name["MT-CO1"] + by_name["MT-ND1"]
    plus_singles = [g for g in genes if g["name"].startswith("G0")]

    molecules = []   # (cb, umi, kind, gene, counted gene_id or None, reads: [(cr, ur, gene)])
    for cb in cells + ambient:
        is_cell = cb in cells
        n_normal = rng.randint(120, 200) if is_cell else rng.randint(2, 5)
        kinds = ["normal"] * n_normal
        if is_cell:
            for k, n in SPECIAL.items():
                kinds += [k] * n
            kinds += ["spliced"] * 6 + ["minus_strand"] * 4
        umis = spaced_set(rng, len(kinds), umi_len, 3)
        for kind, umi in zip(kinds, umis):
            if kind in ("normal", "cb_error", "umi_error", "multigene"):
                gene = rng.choice(mito) if rng.random() < 0.05 else rng.choice(countable)
            elif kind == "spliced":
                gene = by_name["GSPL"][0]
            elif kind == "minus_strand":
                gene = by_name["GMINUS"][0]
            elif kind == "antisense":
                gene = rng.choice(plus_singles)
            else:
                gene = None
            n_reads = 1 if not is_cell else rng.randint(1, 4)
            reads = []
            if kind == "cb_error":
                bad_cb = substitute(rng, cb)
                reads = [(bad_cb, umi, gene, kind)] * rng.randint(2, 3)
            elif kind == "umi_error":
                reads = [(cb, umi, gene, kind)] * 3 + [(cb, substitute(rng, umi), gene, kind)]
            elif kind == "multigene":
                other = rng.choice([g for g in countable if g["gene_id"] != gene["gene_id"]])
                reads = [(cb, umi, gene, "normal")] * 3 + [(cb, umi, other, "normal")]
            elif kind == "spliced":
                reads = [(cb, umi, gene, "spliced")] * n_reads
            else:
                reads = [(cb, umi, gene, kind)] * n_reads
            counted = gene["gene_id"] if kind in ("normal", "cb_error", "umi_error", "multigene",
                                                    "spliced", "minus_strand") else None
            molecules.append(dict(cb=cb, umi=umi, kind=kind,
                                  gene=gene["gene_id"] if gene else None,
                                  other_gene=reads[-1][2]["gene_id"] if kind == "multigene" else None,
                                  counted=counted, reads=reads))

    rg_id = f"{sample}:0:1:{FLOWCELL}:1"
    header = [
        "@HD\tVN:1.6\tSO:unsorted",
        f"@RG\tID:{rg_id}\tSM:{sample}\tLB:0.1\tPU:{rg_id}\tPL:ILLUMINA",
        "@CO\t10x_bam_to_fastq:R1(CR:CY,UR:UY)",
        "@CO\t10x_bam_to_fastq:R2(SEQ:QUAL)",
    ]
    if with_seqnames:
        header.append("@CO\t10x_bam_to_fastq_seqnames:R1,R3,I1,R2")
    records, n = [], 0
    for m in molecules:
        for cr, ur, gene, kind in m["reads"]:
            n += 1
            seq = read_from(rng, by_name, seqs, regions, kind, gene, rlen)
            tags = [f"CR:Z:{cr}", f"CY:Z:{'F' * len(cr)}", f"UR:Z:{ur}", f"UY:Z:{'F' * len(ur)}", f"RG:Z:{rg_id}"]
            records.append("\t".join([f"{FLOWCELL}:1:1101:{n}:{n}", "4", "*", "0", "0", "*", "*", "0", "0",
                                      seq, "F" * rlen] + tags))
    rng.shuffle(records)
    sam = out / f"{sample}.sam"
    sam.write_text("\n".join(header + records) + "\n")
    subprocess.run([samtools, "view", "--no-PG", "-b", "-o", str(out / f"{sample}.bam"), str(sam)], check=True)
    sam.unlink()

    expected = {}
    for m in molecules:
        if m["counted"]:
            e = expected.setdefault(m["cb"], {})
            e[m["counted"]] = e.get(m["counted"], 0) + 1
    # Keep the special molecules for diagnosing a mismatch; normal ones are
    # fully described by `expected`.
    special = [{k: v for k, v in m.items() if k != "reads"} for m in molecules if m["kind"] != "normal"]
    return dict(published_label=label, chemistry=chem, cdna_read="R3" if with_seqnames else "R2",
                n_reads=n, cells=cells, expected=expected, special_molecules=special)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(Path(__file__).parent / "accuracy"))
    ap.add_argument("--samtools", default=shutil.which("samtools") or "samtools")
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    contigs, genes, seqs, regions = build_genome(rng)
    write_fasta(contigs, out / "genome.fa")
    write_gtf(genes, out / "genes.gtf")

    truth, sheet = {}, ["sample,bam,chemistry,published_label"]
    for spec in SAMPLES:
        sample, label, chem = spec[:3]
        used = spaced_set(rng, N_CELLS + N_AMBIENT, 16, 3)
        decoys = spaced_set(rng, N_DECOYS, 16, 3, avoid=used)
        (out / f"whitelist_{chem}.txt").write_text("\n".join(sorted(used + decoys)) + "\n")
        truth[sample] = make_sample(rng, spec, used, genes, seqs, regions, args.samtools, out)
        sheet.append(f"{sample},{sample}.bam,{chem},{label}")
        print(f"{sample}: {chem}, {truth[sample]['n_reads']} reads, cDNA read {truth[sample]['cdna_read']}")
    (out / "samplesheet.csv").write_text("\n".join(sheet) + "\n")

    # Published matrix: true cells only, gene symbols made unique in GTF order
    # (TWIN, TWIN.1), as Seurat's Read10X writes them.
    names, seen = [], {}
    for g in genes:
        k = seen.get(g["name"], 0)
        names.append(g["name"] if k == 0 else f"{g['name']}.{k}")
        seen[g["name"]] = k + 1
    cols = [(t["published_label"], cb, t["expected"].get(cb, {})) for t in truth.values() for cb in t["cells"]]
    with open(out / "published_counts.csv", "w") as fh:
        fh.write("," + ",".join(f"{label}_{cb}" for label, cb, _ in cols) + "\n")
        for g, name in zip(genes, names):
            fh.write(name + "," + ",".join(str(e.get(g["gene_id"], 0)) for _, _, e in cols) + "\n")
    (out / "truth.json").write_text(json.dumps(truth, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
