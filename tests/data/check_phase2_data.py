#!/usr/bin/env python3
"""Check that a phase 2 dataset contains what its truth.json says was planted.

Usage: python3 check_phase2_data.py DIR [DIR ...]

For each directory written by make_phase2_data.py (signal or null), on
counts normalised per cell (count / total x 10,000):
  - matrix, genes.tsv, cells.tsv and truth.json agree in size; cell IDs are
    <library>_<16-nt barcode>; every library maps to the patient truth.json
    gives it, and exactly one patient has two libraries, the second with a
    label that is not the patient's;
  - each lineage's markers are >= 4x higher in that lineage than elsewhere;
  - the malignant program is >= 4x higher in malignant than normal cells;
  - each planted arm gain (loss) raises (lowers) the arm in a library's
    malignant cells relative to its normal cells, both measured against the
    arms without an event (summed raw counts; ratio > 1.25, < 0.75), in
    every library of the patient; the patient with no copy-number change
    stays within 0.8-1.25 on every arm;
  - signal: each subgroup's genes are >= 2.5x higher in its malignant cells
    than in other subgroups' (summed ratio; per patient, genes on that
    patient's copy-number arms are left out); null: within 0.7-1.45x;
  - low-quality cells have < half their library's median total and > 2x
    the median mito fraction;
  - patient effects of |log fold| >= 0.8 go the planted way in >= 90% of
    such genes overall and >= 70% within each patient, and within each
    library of the two-library patient (the shared signal that makes
    library-wise validation leak);
  - the second library's own effect goes the planted way against the
    patient's first library in >= 80% of its genes.
Exits 1 if any check fails. Standard library only.
"""
import gzip
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    lines = [l.rstrip("\n").split("\t") for l in open(path)]
    return [dict(zip(lines[0], r)) for r in lines[1:]]


def load(d):
    d = Path(d)
    genes = read_tsv(d / "genes.tsv")
    cells = read_tsv(d / "cells.tsv")
    truth = json.load(open(d / "truth.json"))
    with gzip.open(d / "matrix.mtx.gz", "rt") as f:
        header = f.readline()
        assert header.startswith("%%MatrixMarket matrix coordinate integer general"), header
        ng, nc, nnz = map(int, f.readline().split())
        counts = [defaultdict(int) for _ in range(nc)]          # per cell: gene index -> count
        n = 0
        for line in f:
            i, j, k = map(int, line.split())
            counts[j - 1][i - 1] = k
            n += 1
    return genes, cells, truth, counts, (ng, nc, nnz, n)


def check_dir(d):
    genes, cells, truth, counts, (ng, nc, nnz, n) = load(d)
    fails = []

    def check(ok, what):
        if not ok:
            fails.append(what)

    check(ng == len(genes) == truth["n_genes"], "gene counts disagree")
    check(nc == len(cells) == truth["n_cells"], "cell counts disagree")
    check(nnz == n, "matrix declares %d entries, has %d" % (nnz, n))
    check(all(re.fullmatch(r"P\d\d_[ACGT]{16}", c["cell_id"]) for c in cells), "cell IDs are not <library>_<barcode>")
    check(all(c["cell_id"].split("_")[0] == c["library"] for c in cells), "a cell ID prefix is not its library")
    check(len({c["cell_id"] for c in cells}) == len(cells), "duplicate cell IDs")

    patients, libraries = truth["patients"], truth["libraries"]
    check(all(libraries.get(c["library"], {}).get("patient") == c["patient"] for c in cells),
          "cells.tsv library -> patient disagrees with truth.json")
    check(all(libraries[lib]["patient"] == pid for pid, p in patients.items() for lib in p["libraries"]),
          "truth.json patients and libraries disagree")
    multi = {pid: p["libraries"] for pid, p in patients.items() if len(p["libraries"]) > 1}
    check(len(multi) == 1, "expected exactly one patient with two libraries, found %d" % len(multi))
    check(all(lib != pid for pid, libs in multi.items() for lib in libs[1:]),
          "a second library's label is the patient's (it must not reveal the patient)")
    check({c["library"] for c in cells} == set(libraries), "libraries in cells.tsv and truth.json differ")

    gidx = {g["gene_name"]: i for i, g in enumerate(genes)}
    totals = [sum(c.values()) for c in counts]
    norm = [{i: v * 1e4 / t for i, v in c.items()} if t else {} for c, t in zip(counts, totals)]
    ok_cells = [j for j, c in enumerate(cells) if c["low_quality"] == "FALSE"]

    def mean(gene_i, js):
        return sum(norm[j].get(gene_i, 0.0) for j in js) / max(len(js), 1)

    def ratio(gene_names, js_in, js_out):
        """Geometric mean over genes of mean(in) / mean(out), with a small floor."""
        logs = [math.log((mean(gidx[g], js_in) + 1e-3) / (mean(gidx[g], js_out) + 1e-3)) for g in gene_names]
        return math.exp(sum(logs) / len(logs))

    by_type = defaultdict(list)
    for j in ok_cells:
        by_type[cells[j]["cell_type"]].append(j)
    for t, gs in truth["markers"].items():
        types = ["myeloid", "lymphoid"] if t == "immune" else [t]
        js_in = [j for tt in types for j in by_type[tt]]
        js_out = [j for j in ok_cells if cells[j]["cell_type"] not in types]
        for g in gs:
            r = ratio([g], js_in, js_out)
            check(r >= 4, "marker %s of %s: only %.2fx" % (g, t, r))
    mal = [j for j in ok_cells if cells[j]["malignant"] == "TRUE"]
    normal = [j for j in ok_cells if cells[j]["malignant"] == "FALSE"]
    r = ratio(truth["malignant_program"], mal, normal)
    check(r >= 4, "malignant program: only %.2fx" % r)

    # Copy number: each arm's background genes in a library's malignant
    # cells relative to its normal cells, both measured against the arms with
    # no event in that patient (as inferCNV uses a normal reference). Raw
    # counts, summed, so per-cell normalisation cannot dilute a gain. Every
    # library of a patient must carry the patient's events.
    check(truth["no_cnv_patient"] in patients, "the no-CNV patient %s is not in the data" % truth["no_cnv_patient"])
    arm_genes = defaultdict(list)
    for i, g in enumerate(genes):
        if g["role"] == "background":
            arm_genes[(g["chrom"], g["arm"])].append(i)
    for pid, p in patients.items():
        events = {(e["chrom"], e["arm"]): e["multiplier"] for e in p["cnv"]}
        ref = [i for arm, gs in arm_genes.items() if arm not in events for i in gs]
        for lib in p["libraries"]:
            mal_l = [j for j in mal if cells[j]["library"] == lib]
            nor_l = [j for j in normal if cells[j]["library"] == lib]

            def rel(gs, js):
                return (sum(counts[j].get(i, 0) for j in js for i in gs) /
                        max(sum(counts[j].get(i, 0) for j in js for i in ref), 1))
            for arm, gs in arm_genes.items():
                r = rel(gs, mal_l) / max(rel(gs, nor_l), 1e-9)
                if arm in events:
                    if events[arm] > 1:
                        check(r > 1.25, "%s (library %s) gain on %s%s: ratio %.2f" % (pid, lib, arm[0], arm[1], r))
                    else:
                        check(r < 0.75, "%s (library %s) loss on %s%s: ratio %.2f" % (pid, lib, arm[0], arm[1], r))
                elif pid == truth["no_cnv_patient"]:
                    check(0.8 <= r <= 1.25, "%s (no CNV) %s%s: ratio %.2f" % (pid, arm[0], arm[1], r))
        check(pid != truth["no_cnv_patient"] or not p["cnv"], "the no-CNV patient has events")

    # Subgroup signal as a summed ratio (stable at the few counts of a single
    # WNT patient), per in-group patient leaving out the subgroup genes on
    # that patient's own copy-number arms (they measure dosage, not
    # subgroup), pooled over the subgroup's patients.
    gene_arm = {g["gene_name"]: (g["chrom"], g["arm"]) for g in genes}
    for sg, gs in truth["subgroup_genes"].items():
        js_out = [j for j in mal if cells[j]["subgroup"] != sg]
        num = den = 0.0
        for pid, p in patients.items():
            if p["subgroup"] != sg:
                continue
            arms = {(e["chrom"], e["arm"]) for e in p["cnv"]}
            use = [g for g in gs if gene_arm[g] not in arms]
            js_p = [j for j in mal if cells[j]["patient"] == pid]
            num += sum(mean(gidx[g], js_p) for g in use)
            den += sum(mean(gidx[g], js_out) for g in use)
        r = num / max(den, 1e-9)
        if truth["dataset"] == "signal":
            check(r >= 2.5, "subgroup %s genes: only %.2fx in signal data" % (sg, r))
        else:
            check(0.7 <= r <= 1.45, "subgroup %s genes: %.2fx in the negative control" % (sg, r))

    mito = [i for i, g in enumerate(genes) if g["role"] == "mito"]
    jidx = {c["cell_id"]: j for j, c in enumerate(cells)}
    for cid in truth["low_quality_cells"]:
        j = jidx[cid]
        peers = [k for k in ok_cells if cells[k]["library"] == cells[j]["library"]]
        med_total = sorted(totals[k] for k in peers)[len(peers) // 2]
        fr = lambda k: sum(counts[k].get(i, 0) for i in mito) / max(totals[k], 1)
        med_mito = sorted(fr(k) for k in peers)[len(peers) // 2]
        check(totals[j] < 0.5 * med_total, "low-quality %s: total %d vs median %d" % (cid, totals[j], med_total))
        check(fr(j) > 2 * med_mito, "low-quality %s: mito %.3f vs median %.3f" % (cid, fr(j), med_mito))

    def agreement(effects, js_in, js_out):
        return sum(1 for g, l in effects.items() if (ratio([g], js_in, js_out) > 1) == (l > 0))

    n_strong = n_agree = 0
    for pid, p in patients.items():
        strong = {g: l for g, l in p["patient_effect_log"].items() if abs(l) >= 0.8}
        js_out = [j for j in ok_cells if cells[j]["patient"] != pid]
        groups = [("all libraries", [j for j in ok_cells if cells[j]["patient"] == pid])]
        if len(p["libraries"]) > 1:
            groups += [("library " + lib, [j for j in ok_cells if cells[j]["library"] == lib]) for lib in p["libraries"]]
        for what, js_in in groups:
            agree = agreement(strong, js_in, js_out)
            check(not strong or agree >= 0.7 * len(strong),
                  "%s (%s): patient effect: %d of %d strong genes go the planted way" % (pid, what, agree, len(strong)))
            if what == "all libraries":
                n_strong += len(strong)
                n_agree += agree
    check(n_agree >= 0.9 * n_strong, "patient effects: %d of %d strong genes go the planted way" % (n_agree, n_strong))

    for pid, libs in multi.items():
        first = [j for j in ok_cells if cells[j]["library"] == libs[0]]
        for lib in libs[1:]:
            eff = libraries[lib]["library_effect_log"]
            check(len(eff) > 0, "second library %s has no library effect" % lib)
            js_in = [j for j in ok_cells if cells[j]["library"] == lib]
            agree = agreement(eff, js_in, first)
            check(agree >= 0.8 * len(eff), "library %s effect: %d of %d genes go the planted way against %s"
                  % (lib, agree, len(eff), libs[0]))
    return fails


def main():
    bad = 0
    for d in sys.argv[1:]:
        fails = check_dir(d)
        if fails:
            bad += 1
            print("FAIL %s:" % d)
            for f in fails:
                print("  - " + f)
        else:
            print("ok   %s: planted structure present" % d)
    sys.exit(1 if bad or not sys.argv[1:] else 0)


if __name__ == "__main__":
    main()
