#!/usr/bin/env python3
"""Generate the phase 2 known-answer dataset (plan M5 step 1; DRAFT).

A small, fully synthetic single-cell count matrix with planted structure,
so each phase 2 step can be tested on whether it recovers what was put in:

  cell types     malignant cells plus four normal lineages (myeloid,
                 lymphoid, endothelial, oligodendrocyte), each with marker
                 genes ten times higher than elsewhere
  patients       13 synthetic patients: 4 SHH, 4 Group 3, 4 Group 4, 1 WNT
                 (the real cohort's shape: WNT cannot be learned and tested)
  patient effect every patient shifts 30 background genes in all its cells
                 (a batch effect integration should remove)
  copy number    each patient's malignant cells carry 1-2 arm-level gains
                 (x1.5) or losses (x0.5) along gene order; one patient
                 carries none
  subgroup       malignant cells of each subgroup express that subgroup's
                 genes five times higher (absent in the negative control)
  low quality    two cells per patient with few counts and high mito

The negative control (--null) is the same dataset with the subgroup signal
removed. Every random number is drawn in the same order in both, and the
Poisson sampler uses exactly one uniform per count, so the two differ only
through the subgroup genes' expected counts.

Outputs (in --outdir; tests/data/phase2 and tests/data/phase2_null):
  matrix.mtx.gz  genes x cells integer counts, Matrix Market coordinate
                 (gzip with mtime 0)
  genes.tsv      gene_id, gene_name, chrom, start, end, arm, role
  cells.tsv      cell_id (<patient>_<barcode>, as in the published matrix),
                 patient, subgroup, cell_type, malignant, low_quality
  truth.json     everything that was planted

Usage: python3 make_phase2_data.py [--outdir DIR] [--null]
Standard library only. Only Random.random() is used for randomness (its
stream is stable across Python versions); normals come from Box-Muller.
"""
import argparse
import gzip
import json
import math
import random
from pathlib import Path

SEED = 20261101
N_CHROM, GENES_PER_CHROM = 4, 60          # arms: first and second half of each chromosome
CELLS_PER_PATIENT = 60
COMPOSITION = [("malignant", 40), ("myeloid", 8), ("lymphoid", 6), ("endothelial", 3), ("oligodendrocyte", 3)]
TYPE_SCALE = {"malignant": 1.3, "myeloid": 1.0, "lymphoid": 0.6, "endothelial": 0.8, "oligodendrocyte": 0.9}
PATIENTS = [("P%02d" % (i + 1), sg) for i, sg in enumerate(["SHH"] * 4 + ["G3"] * 4 + ["G4"] * 4 + ["WNT"])]
NO_CNV_PATIENT = "P12"
N_LOW_QUALITY = 2

# Marker symbols follow the proposals in docs/phase2_design.md (D3 lineage
# markers, D7 Northcott 2012 subgroup genes); see the PR for what depends
# on those decisions.
MARKERS = {
    "myeloid": ["CD14", "CD68", "LYZ", "CSF1R", "C1QA"],
    "lymphoid": ["CD3E", "CD2", "NKG7"],
    "endothelial": ["PECAM1", "CLDN5", "VWF"],
    "oligodendrocyte": ["MBP", "PLP1", "OLIG2"],
}
IMMUNE = ["PTPRC"]                              # myeloid and lymphoid
SUBGROUP_MARKERS = {
    "SHH": ["PDLIM3", "EYA1", "HHIP", "ATOH1", "SFRP1"],
    "G3": ["IMPG2", "GABRA5", "EYS", "NRL", "MAB21L2", "NPR3"],
    "G4": ["KCNA1", "EOMES", "KHDRBS2", "RBM24", "UNC5D", "OAS1"],
    "WNT": ["WIF1", "TNC", "GAD1", "DKK2", "EMX2"],
}
N_SUBGROUP_EXTRA = 8                            # unnamed signal genes per subgroup, beyond the markers
N_MALIGNANT_PROGRAM = 6
MITO = ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CYB"]
N_PATIENT_GENES = 30

MARKER_FOLD, SUBGROUP_FOLD, MALIGNANT_FOLD = 10.0, 5.0, 8.0
GAIN, LOSS = 1.5, 0.5
NOISE_SD = 0.25          # per cell x gene log-normal noise
PATIENT_LFC_SD = 0.6     # patient effect: log fold change ~ N(0, 0.6)
LAMBDA_CAP = 400.0
LQ_SIZE, LQ_MITO = 0.1, 25.0  # low-quality cells: library size and mito expression


class Rng:
    """Random numbers from Random.random() only, so the stream is stable."""

    def __init__(self, seed):
        self.r = random.Random(seed)

    def u(self):
        return self.r.random()

    def normal(self):
        u1 = 1.0 - self.u()                     # (0, 1]
        return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * self.u())

    def below(self, n):
        return min(int(self.u() * n), n - 1)

    def sample(self, items, k):
        pool, out = list(items), []
        for _ in range(k):
            out.append(pool.pop(self.below(len(pool))))
        return out

    def poisson(self, lam):
        """Inversion with a single uniform, whatever lambda is."""
        u = self.u()
        k, p = 0, math.exp(-lam)
        cdf = p
        while u > cdf and k < 10000:
            k += 1
            p *= lam / k
            cdf += p
        return k


def build_genes(rng):
    named = (MITO + sum(MARKERS.values(), []) + IMMUNE + sum(SUBGROUP_MARKERS.values(), []))
    extra = {sg: ["%s_SIG%02d" % (sg, i + 1) for i in range(N_SUBGROUP_EXTRA)] for sg in SUBGROUP_MARKERS}
    program = ["MALIG%02d" % (i + 1) for i in range(N_MALIGNANT_PROGRAM)]
    nuclear_named = [g for g in named if g not in MITO] + sum(extra.values(), []) + program
    n_slots = N_CHROM * GENES_PER_CHROM
    order = rng.sample(range(n_slots), n_slots)           # random placement along the genome
    names = nuclear_named + ["GENE%03d" % (i + 1) for i in range(n_slots - len(nuclear_named))]
    role = {}
    for t, gs in MARKERS.items():
        role.update({g: "marker:" + t for g in gs})
    role.update({g: "marker:immune" for g in IMMUNE})
    for sg, gs in SUBGROUP_MARKERS.items():
        role.update({g: "subgroup:" + sg for g in gs + extra[sg]})
    role.update({g: "malignant_program" for g in program})
    genes = []
    for slot, name in zip(order, names):
        c, i = divmod(slot, GENES_PER_CHROM)
        genes.append(dict(slot=slot, name=name, chrom="chr%d" % (c + 1), start=100000 * (i + 1),
                          arm="p" if i < GENES_PER_CHROM // 2 else "q", role=role.get(name, "background")))
    genes.sort(key=lambda g: g["slot"])
    for i, name in enumerate(MITO):
        genes.append(dict(slot=n_slots + i, name=name, chrom="chrM", start=1000 * (i + 1), arm="-", role="mito"))
    for i, g in enumerate(genes):
        g["id"] = "SYNG%05d" % (i + 1)
        g["end"] = g["start"] + 999
    return genes, extra, program


def generate(null):
    rng = Rng(SEED)
    genes, extra, program = build_genes(rng)
    idx = {g["name"]: i for i, g in enumerate(genes)}
    # Baseline expected counts per gene (before cell size and multipliers).
    base = []
    for g in genes:
        mu = math.exp(0.2 + rng.normal())
        if g["role"].startswith("marker") or g["role"].startswith("subgroup") or g["role"] == "malignant_program":
            mu = 0.15
        if g["role"] == "mito":
            mu = 1.5
        base.append(min(max(mu, 0.05), 8.0))
    background = [i for i, g in enumerate(genes) if g["role"] == "background"]
    arms = [("chr%d" % (c + 1), a) for c in range(N_CHROM) for a in "pq"]

    patients, cells = {}, []
    for pid, sg in PATIENTS:
        lfc = {genes[i]["name"]: round(PATIENT_LFC_SD * rng.normal(), 4) for i in rng.sample(background, N_PATIENT_GENES)}
        n_events = 0 if pid == NO_CNV_PATIENT else 1 + rng.below(2)
        events = []
        for chrom, arm in rng.sample(arms, n_events):
            mult = GAIN if rng.u() < 0.5 else LOSS
            events.append(dict(chrom=chrom, arm=arm, multiplier=mult,
                               genes=[g["name"] for g in genes if g["chrom"] == chrom and g["arm"] == arm]))
        patients[pid] = dict(subgroup=sg, cnv=events, patient_effect_log=lfc)
        types = sum([[t] * n for t, n in COMPOSITION], [])
        low = set(rng.sample(range(CELLS_PER_PATIENT), N_LOW_QUALITY))
        used = set()
        for k, t in enumerate(types):
            while True:
                bc = "".join("ACGT"[rng.below(4)] for _ in range(16))
                if bc not in used:
                    used.add(bc)
                    break
            cells.append(dict(id="%s_%s" % (pid, bc), patient=pid, subgroup=sg, type=t,
                              malignant=t == "malignant", low_quality=k in low,
                              size=math.exp(0.35 * rng.normal()) * TYPE_SCALE[t]))

    entries = []                                   # (gene index, cell index, count)
    for j, c in enumerate(cells):
        p = patients[c["patient"]]
        mult = [1.0] * len(genes)
        for name, l in p["patient_effect_log"].items():
            mult[idx[name]] *= math.exp(l)
        if c["type"] in MARKERS:
            for name in MARKERS[c["type"]]:
                mult[idx[name]] *= MARKER_FOLD
        if c["type"] in ("myeloid", "lymphoid"):
            for name in IMMUNE:
                mult[idx[name]] *= MARKER_FOLD
        if c["malignant"]:
            for name in program:
                mult[idx[name]] *= MALIGNANT_FOLD
            if not null:
                for name in SUBGROUP_MARKERS[c["subgroup"]] + extra[c["subgroup"]]:
                    mult[idx[name]] *= SUBGROUP_FOLD
            for e in p["cnv"]:
                for name in e["genes"]:
                    mult[idx[name]] *= e["multiplier"]
        size = c["size"] * (LQ_SIZE if c["low_quality"] else 1.0)
        for i, g in enumerate(genes):
            m = mult[i] * (LQ_MITO if c["low_quality"] and g["role"] == "mito" else 1.0)
            lam = min(size * base[i] * m * math.exp(NOISE_SD * rng.normal()), LAMBDA_CAP)
            k = rng.poisson(lam)
            if k:
                entries.append((i + 1, j + 1, k))
    truth = dict(
        dataset="null" if null else "signal", seed=SEED, draft=True,
        folds=dict(marker=MARKER_FOLD, subgroup=0.0 if null else SUBGROUP_FOLD, malignant_program=MALIGNANT_FOLD,
                   gain=GAIN, loss=LOSS, low_quality_size=LQ_SIZE, low_quality_mito=LQ_MITO),
        markers=dict(MARKERS, immune=IMMUNE),
        subgroup_genes={sg: SUBGROUP_MARKERS[sg] + extra[sg] for sg in SUBGROUP_MARKERS},
        malignant_program=program, no_cnv_patient=NO_CNV_PATIENT, patients=patients,
        low_quality_cells=[c["id"] for c in cells if c["low_quality"]],
        n_genes=len(genes), n_cells=len(cells))
    return genes, cells, entries, truth


def write(outdir, genes, cells, entries, truth):
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "matrix.mtx.gz", "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
            gz.write(b"%%MatrixMarket matrix coordinate integer general\n")
            gz.write(("%d %d %d\n" % (len(genes), len(cells), len(entries))).encode())
            for e in entries:
                gz.write(("%d %d %d\n" % e).encode())
    with open(outdir / "genes.tsv", "w") as f:
        f.write("gene_id\tgene_name\tchrom\tstart\tend\tarm\trole\n")
        for g in genes:
            f.write("\t".join(str(g[k]) for k in ("id", "name", "chrom", "start", "end", "arm", "role")) + "\n")
    with open(outdir / "cells.tsv", "w") as f:
        f.write("cell_id\tpatient\tsubgroup\tcell_type\tmalignant\tlow_quality\n")
        for c in cells:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (c["id"], c["patient"], c["subgroup"], c["type"],
                                                  "TRUE" if c["malignant"] else "FALSE",
                                                  "TRUE" if c["low_quality"] else "FALSE"))
    with open(outdir / "truth.json", "w") as f:
        json.dump(truth, f, indent=1, sort_keys=True)
        f.write("\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--null", action="store_true", help="the negative control: no subgroup signal")
    ap.add_argument("--outdir", type=Path, default=None)
    a = ap.parse_args()
    here = Path(__file__).resolve().parent
    outdir = a.outdir or here / ("phase2_null" if a.null else "phase2")
    write(outdir, *generate(a.null))
    print("wrote", outdir)


if __name__ == "__main__":
    main()
