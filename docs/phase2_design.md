# Phase 2 design (draft for decision)

Status: **draft, nothing decided.** This is plan M4 step 1. Every
"recommendation" below is a proposal; Chia-Lung chooses (M4 step 2). The
thresholds that will judge the results are not set here: they go into
`preregistration/phase2_criteria.yaml` (M4 step 3), drafted from the
published cell metadata's column *names* only and locked before any of its
values is read (CLAUDE.md rule 1).

Already settled and not reopened here (CLAUDE.md, 31 Oct 2026): the main
analysis is R (Seurat + Bioconductor); a Python module answers exactly two
questions, the scVI-versus-Harmony integration comparison and the subgroup
classifier; phase 1 is not ported; phase 2 analyses all 28 tumors from the
published matrix.

Every compute figure in this document is an **estimate, not a
measurement**. The M5 pilot on the cluster replaces them with `trace.tsv`
numbers, as M3 does for phase 1.

## Decisions needed

| # | Step | Options | Recommendation |
|---|---|---|---|
| D1 | QC | no re-filter / re-filter with phase 1 MAD rules / flag only | flag only; sensitivity analysis without the flagged cells |
| D2 | Normalisation | LogNormalize / SCTransform v2 | LogNormalize |
| D3 | Annotation | markers / SingleR (which reference) / both | markers for lineage, SingleR + Monaco for immune subtypes |
| D4 | Malignant cells | inferCNV / copykat / numbat; normal reference | inferCNV on all 28 with a pooled immune + endothelial reference; numbat on the 8 phase 1 tumors as an optional check |
| D5 | Tumor programs | none / score published programs / de novo NMF | score published programs; NMF only as labelled exploration |
| D6 | Integration comparison | cells, embeddings, metrics, decision rule shape, seeds | non-malignant cells only; none / Harmony / scVI; scib-metrics; 5 seeds |
| D7 | Subgroup classifier | pseudobulk / cell-level; tuning | cell-level logistic regression, fixed hyperparameters, leave-one-patient-out |
| D8 | R to Python hand-off | Matrix Market + TSV / `.h5ad` | Matrix Market + TSV with a column schema |
| D9 | Python layout and image | pip/uv lock / conda | `python/` package, uv lock, `python:3.12-slim` image |
| D10 | Wiring | `--phase` switch in `main.nf` / separate entry script | separate subworkflows, selected by a required `--phase` |

D9 is engineering, with no effect on results; I list it so it is visible,
but it does not need a decision.

## Inputs and blinding

- **Counts:** `GSE155446_human_raw_counts.csv.gz`, read with the phase 1
  `PREP_PUBLISHED` step (already tested: streams the dense CSV into a sparse
  matrix in chunks and detects its orientation). The plan's figure is about
  40,000 cells; the exact dimensions come from the header in M5, which
  blinding allows.
- **Tumor of each cell:** from the cell-ID prefix (names are allowed), as
  phase 1's `published_labels.tsv` does.
- **Subgroup of each tumor:** from the paper's sample table (Riemondy *et
  al.* 2022), not from the cell metadata file. The cohort is 1 WNT, 9 SHH,
  7 Group 3, 11 Group 4. Before M5 the mapping is committed as
  `docs/results/phase2/tumor_subgroups.tsv` with the paper table it came
  from.
- **Gene positions** (needed for copy-number inference): from the
  GRCh38-3.0.0 GTF that phase 1 uses, matched by symbol in GTF order, the
  same way phase 1 matches genes (symbols made unique as `RGS5.1` included).
- **Cell metadata** (`GSE155446_human_cell_metadata.csv.gz`): its values
  stay unread until `phase2_criteria.yaml` is locked and released. After
  that, its labels serve only as the reference that the criteria compare
  against: they are never inputs to annotation, malignant calls,
  integration or the classifier.

## R: main analysis

### D1. QC: re-filter or not

The authors already filtered the published matrix. Their thresholds are in
the paper's methods, not re-derived here.

| Option | For | Against |
|---|---|---|
| a. Use the published cells as they are | Same cells as the published labels, so agreement is measured on identical sets; nothing to tune | Doublets or low-quality cells that the authors kept stay in |
| b. Re-filter with the phase 1 rules (median ± 3 MAD per tumor on log UMI, log genes, % mito) | Same rule as phase 1 | MADs on an already-filtered distribution are narrower, so a second pass removes real cells. Malignant and cycling cells sit in the tails (high UMI, high mito) and would be removed preferentially. The removed cells then have no counterpart in the label comparison |
| c. Flag with the phase 1 rules, remove nothing (**recommended**) | Reuses tested code (`qc_metrics`, `qc_flag`); reports how many cells each tumor would lose; keeps (a)'s comparability | One more output to read |

Recommendation (c): the main analysis uses every published cell. In
addition, the annotation-agreement and classifier results are re-run
without the flagged cells as a sensitivity analysis and reported next to
the main result, whatever it shows. Doublet scoring (scDblFinder) could be
added in the same flag-only way; I do not propose it as a filter.

- **Known-answer test:** the phase 1 generator already plants low-quality
  cells, and `qc_flag` is tested against scuttle and against median ± n MAD
  computed by hand. New: a pipeline test that the phase 2 cell count equals
  the input cell count (nothing silently dropped), and that the flags equal
  the planted low-quality cells.
- **Compute:** minutes; memory dominated by loading the matrix (below).

### D2. Normalisation

| Option | For | Against |
|---|---|---|
| a. `LogNormalize` (log1p of counts per 10,000) (**recommended**) | A closed formula: an independent implementation in base R and in Python gives the exact answer, so it can have an exact accuracy test like phase 1's. Identical to scanpy's `normalize_total` + `log1p`, so R and Python see the same values. Deterministic | Depth effects remain in highly expressed genes; variance not stabilised |
| b. `SCTransform` v2 (glmGamPoi) | Better variance stabilisation; usually cleaner clusters on cells with very different depths | Output depends on a fitted model, so there is no exact known answer to test against. Slower, more memory. Residuals are the wrong input for DE and for the Python classifier, so two representations would be in play |

Recommendation (a). If clusters look depth-driven in M5, SCTransform v2 can
be run as a sensitivity analysis on clustering only, reported as such.

- **Known-answer test:** `lib/R` normalisation equals
  `log1p(1e4 * x / colSums(x))` from a dense base-R implementation on random
  sparse matrices, including all-zero genes and a one-count cell; the
  Python module's normalisation of the handed-off counts equals R's to
  within 1e-12.
- **Compute:** seconds to minutes.

### Shared: features, PCA, clustering

These are not separate decisions, but they feed D3, D4 and D6, so they are
fixed in the design rather than tuned later: 2,000 highly variable genes
(Seurat `vst`), PCA with 30 components, a shared-nearest-neighbour graph
and Leiden clustering at a resolution set in the design (proposal: 0.8),
seeds fixed. Mitochondrial, ribosomal and cell-cycle genes stay in (they
are flagged, not removed, as in D1).

- **Known-answer test:** on the generator's planted cell types (M5 step 1),
  clustering recovers them (ARI against the planted types above a bound
  set in the test, e.g. 0.9 on well-separated synthetic types).
- **Compute:** 10–20 minutes on one core; peak memory perhaps 16–32 GB for
  about 40,000 cells × 33,000 genes. Loading the dense CSV is the
  unmeasured part (the README's known limitation for `CONCORDANCE`
  applies here too).

### D3. Annotation

| Option | For | Against |
|---|---|---|
| a. Marker-based: pre-specified marker lists, module scores per cluster, a written rule from scores to labels | Transparent; matches how medulloblastoma papers annotate; works for brain cell types that general references lack | Subjective unless the rule is written down before the data is seen; cluster-level labels hide mixed clusters |
| b. Reference-based, SingleR with a celldex reference (HumanPrimaryCellAtlas, BlueprintEncode or MonacoImmune) | Per-cell, reproducible, no hand curation | Bulk-derived references have no medulloblastoma or developing-cerebellum states: tumor cells get forced into the nearest normal type (neurons, astrocytes). The reference must be baked into the image, not downloaded at run time |
| c. Reference-based with a developing-cerebellum single-cell atlas | Closest normal counterpart of medulloblastoma cells | Large; licence and version to check; mapping tumor cells to developmental states is a research question in itself, not an annotation step |
| d. Two levels (**recommended**) | Each method used where it is reliable | Two methods to test and describe |

Recommendation (d):

1. **Lineage** (malignant-candidate, myeloid, lymphoid, endothelial,
   pericyte/fibroblast, astrocyte, oligodendrocyte lineage) by module
   scores of a marker list fixed in the design. Proposal, to be checked
   against the paper before locking: PTPRC; myeloid CD14, CD68, LYZ,
   CSF1R, C1QA; T/NK CD3E, CD2, NKG7; B CD79A, MS4A1; endothelial PECAM1,
   CLDN5, VWF; pericyte/fibroblast PDGFRB, RGS5, COL1A1; astrocyte GFAP,
   AQP4; oligodendrocyte lineage MBP, PLP1, OLIG2. The rule: a cluster takes
   the lineage whose mean score is highest by a margin set in the design;
   otherwise "unresolved".
2. **Immune subtypes** by SingleR with `celldex::MonacoImmuneData`, version
   pinned and stored in the R image, run on immune cells only.
3. **Malignant or not** is decided by copy number (D4), not by annotation.

Agreement with the published labels (after the lock) is a phase 2
criterion; its threshold is set in M4 step 3.

- **Known-answer test:** generator with planted cell types whose marker
  genes are known; the rule recovers them, and returns "unresolved" for a
  planted cluster that expresses no marker set. SingleR is checked on a
  reference built from the generator's own type profiles. A test fixes the
  rule's output for a given seed, so a change in it shows.
- **Compute:** module scores in minutes; SingleR on perhaps 5,000–15,000
  immune cells against Monaco in 5–15 minutes on 4 cores.

### D4. Copy-number inference: which cells are malignant

| Option | For | Against |
|---|---|---|
| a. inferCNV (Bioconductor) | The most used method; per-tumor heatmaps are the standard figure; works from expression only, so on all 28 tumors | Needs a normal reference; many parameters (cutoff 0.1 for 10x, denoising, HMM model); the HMM step is slow and memory-hungry; malignant calls come from clustering the CNV profiles, which needs a rule |
| b. copykat | Finds diploid cells itself; per-cell aneuploid/diploid call; cheaper than the inferCNV HMM | Weak on tumors with few copy-number changes, which some medulloblastomas are; GitHub-only package (pin a commit) |
| c. numbat | Haplotype-aware (allele counts plus expression), much more specific than expression alone | Needs BAMs, so only the 8 phase 1 tumors; needs cellsnp-lite pileups and phasing against a 1000 Genomes hg38 panel (several GB); the most expensive by far |

**Normal reference.** Options: per-tumor normal cells (most faithful, but
some tumors may have few); pooled normal cells from all tumors; or
immune-only. Recommendation: pooled immune and endothelial cells from all
tumors, chosen by the lineage rule (D3 step 1) before any CNV call, with
myeloid and lymphoid kept as separate reference groups (inferCNV's
`ref_group_names`). Glial cells are not in the reference, because
annotation can confuse them with tumor cells.

**Malignant call.** Per tumor: each cell's CNV score (mean squared
deviation of the smoothed profile) and correlation with the tumor's mean
profile, the two-axis rule of Tirosh *et al.* 2016 / Hovestadt *et al.*
2019, with cut-offs set in the design and recorded per tumor. Cells that
fall between the two clusters are "unresolved", not forced.

Recommendation: (a) on all 28 tumors (without the HMM, which the malignant
call does not need, unless arm-level calls are wanted for the write-up),
plus (c) on the 8 phase 1 tumors as an **optional** orthogonal check,
reported as agreement between the two methods' calls. Expected arm-level
events from the literature (for example 17p loss / 17q gain in Groups 3
and 4, chromosome 6 loss in WNT) are a sanity check only, never a
criterion: the subgroup is what the classifier predicts.

- **Known-answer test:** generator with genes ordered along a few
  synthetic chromosomes; malignant cells in each tumor carry planted gains
  (expression × 1.5) and losses (× 0.5) over chosen arms; normal cells are
  diploid. Tests: the planted malignant cells are called malignant and the
  normal ones not (sensitivity and specificity above bounds set in the
  test); a tumor with no planted change produces no malignant call rather
  than an arbitrary split; arm-level output, if used, finds the planted
  arms.
- **Compute:** inferCNV without HMM, perhaps 10–30 minutes and 8–32 GB per
  tumor, so 28 LSF jobs and about 10–30 CPU hours in total; with the i6 HMM,
  several times that. numbat: pileup about 1–3 hours on 8 cores per BAM,
  and numbat itself several hours on 8 cores per tumor, so of the order of
  hundreds of CPU hours for 8 tumors, plus storage for the phasing panel.

### D5. Tumor programs

| Option | For | Against |
|---|---|---|
| a. None | Not needed by any criterion or by the Python questions | The analysis says nothing about intra-tumor heterogeneity, the paper's main topic |
| b. Score published programs (**recommended**) | Gene lists fixed before the data is seen (Hovestadt *et al.* 2019 programs; the paper's own); no discovery, so nothing to tune; testable | Can only confirm known programs |
| c. De novo per-tumor NMF and meta-programs (Gavish *et al.* 2023 style: k = 4–9 per tumor, robust programs, clustering by top-50-gene Jaccard) | The standard discovery method; good for the write-up | Many parameters; results are hard to validate; easy to over-interpret |

Recommendation (b) in the pre-registered analysis. If (c) is done, it is
labelled exploratory in the report and has no criterion.

- **Known-answer test:** generator with planted programs (gene sets
  co-varying in subsets of malignant cells); scoring separates the planted
  cells; for (c), NMF recovers each planted program (top-gene Jaccard above
  a bound) and finds none in a control with no planted program.
- **Compute:** scoring in minutes. NMF with RcppML at k = 4–9 and 10
  restarts per tumor: perhaps 1–4 CPU hours in total.

## Python module

### D6. Integration comparison: Harmony versus scVI

**Cells.** Non-malignant cells only (D4), all tumors. Malignant cells are
patient-specific. Integrating them by patient would remove real
between-tumor differences, including the subgroup signal the classifier
needs. They are not integrated: they are analysed per tumor, and the
classifier uses their unintegrated log-normalised expression (D7). Cells
called "unresolved" in D4 are excluded and counted.

**Batch key.** The sequencing library, which is the tumor except where a
tumor has more than one library; the design records the mapping from
names. **Label key** for biological conservation: the published cell-type
labels (after the lock). Using our own annotation would be circular,
because it was made from one of the embeddings being compared.

**Embeddings compared** (same cells, same 2,000 genes):

1. none: PCA of log-normalised expression (baseline);
2. Harmony on that PCA, run in R, exported (D8);
3. scVI (scvi-tools) on raw counts of the same genes, batch = library.

**Metrics:** scib-metrics. Batch mixing: iLISI, graph connectivity, kBET
(if affordable), PCR comparison. Biological conservation: ARI and NMI of
Leiden clusters against the labels, cLISI, label silhouette. Every metric
has an independent check on toy inputs (ARI and NMI from the contingency
table by hand; LISI on a constructed neighbourhood with a known answer).

**Decision rule** (shape only; the numbers go into `phase2_criteria.yaml`):
prefer scVI over Harmony only if its batch-mixing score is higher by at
least δ_batch **and** its label ARI is lower by no more than δ_bio, both
judged on the median over seeds. Otherwise prefer Harmony, the simpler
method already in the R analysis. If neither beats "none" on batch mixing
by δ_batch, report that integration is not needed for these cells.

**Seeds.** scVI training and Harmony's k-means initialisation are random:
5 seeds each, all reported (median and range), not only the best.

| Option | For | Against |
|---|---|---|
| Hardware: CPU (**recommended**) | No queue wait; about 10,000–20,000 cells is small for scVI | Slower per seed |
| Hardware: GPU queue | Minutes per seed | Queue wait; CUDA build of torch makes the image much larger; results can differ slightly between CPU and GPU |

- **Known-answer test:** generator with planted cell types and a planted
  patient effect (gene-wise log-fold shifts per patient). Tests: "none"
  scores low on mixing; Harmony and scVI raise it; label ARI stays high.
  A control with no patient effect: integration must not lower label ARI
  below a bound. A confounded control (one cell type only in one patient)
  must keep that type separate, which is the over-correction failure.
  Seeds fixed in tests; a test that two runs with the same seed give
  identical embeddings.
- **Compute:** scVI with about 400 epochs on 10,000–20,000 cells × 2,000
  genes, 8 CPU cores: perhaps 15–30 minutes per seed, so 2–3 hours for 5
  seeds (a few minutes per seed on a GPU). Harmony: minutes per seed.
  scib-metrics on 11 embeddings: under an hour. Memory under 16 GB.

### D7. Subgroup classifier

**Question.** Can a tumor's subgroup be predicted from its malignant
cells' expression, for a tumor the model has not seen?

**Unit of validation: the patient.** Leave-one-patient-out (LOPO) over the
27 SHH / Group 3 / Group 4 tumors. No patient's cells are ever split
between training and test. WNT (one tumor) cannot be both learned and
tested: it is excluded from training and evaluation, and its marker-rule
scores and the 3-class model's prediction are reported separately,
without a claim.

| Option | For | Against |
|---|---|---|
| a. Pseudobulk: one mean profile per tumor, 27 rows | Simple; fast, so a large permutation null is cheap | Ignores within-tumor heterogeneity; 27 rows for 3 classes and thousands of genes |
| b. Cell-level model, predictions averaged per patient (**recommended**) | Uses every malignant cell; the question as worded ("from its malignant cells"); same validation unit | Cells of a patient are correlated, so the effective n is still 27; patients with many cells dominate unless capped |

Recommendation (b), with (a) reported as a secondary analysis.

- **Model:** multinomial logistic regression, L2 penalty, standardised
  features, balanced class weights. **Hyperparameters fixed in the design**
  (proposal: C = 1.0), not tuned. Tuning would need nested grouped CV and
  makes the permutation null 10–20 times more expensive; fixing them
  removes tuning as a route to leakage.
- **Features:** log-normalised expression of the 2,000 most variable genes
  **computed inside each training fold, on training patients only**, so the
  held-out patient never influences feature selection.
- **Cells per patient:** at most 300 malignant cells per patient, sampled
  with a fixed seed, so every patient weighs about the same. Patients with
  fewer cells contribute all of theirs; the minimum is reported.
- **Patient prediction:** the class with the highest mean log-probability
  over the patient's cells.
- **Primary metric:** patient-level balanced accuracy over the 27 LOPO
  predictions, with the confusion matrix. Expect errors between Group 3 and
  Group 4 (plan, "Known risks"); the point is a correctly validated
  estimate, not a high number.
- **Baselines.** Majority class (Group 4; balanced accuracy 1/3). A
  known-marker rule with no training: per tumor, mean z-score of each
  subgroup's marker genes in malignant cells, predicted subgroup = argmax.
  Proposal: the Northcott *et al.* 2012 NanoString genes (SHH: PDLIM3,
  EYA1, HHIP, ATOH1, SFRP1; Group 3: IMPG2, GABRA5, EGFL11/EYS, NRL,
  MAB21L2, NPR3; Group 4: KCNA1, EOMES, KHDRBS2, RBM24, UNC5D, OAS1), to be
  checked against that paper before locking.
- **Chance level:** permute subgroup labels **across patients** (never
  across cells), rerun the whole LOPO procedure, 1,000 permutations; the
  p-value is (1 + number of permutations at or above the real score) /
  1,001.
- **Leaky contrast:** the same model with cells split at random (5-fold,
  patients spread across folds), run once, reported next to the real
  result as cell-level and patient-level accuracy, to show how much leakage
  inflates accuracy.
- **Claim rule** (shape only; numbers in `phase2_criteria.yaml`, as in the
  plan's M4 step 3): say "subgroup is predictable from malignant-cell
  expression" only if patient-level balanced accuracy beats the
  permutation null at p < 0.05 **and** beats the marker baseline.

**Known-answer tests** (plan M5 step 1, written before the classifier):

- planted subgroup signal (subgroup-specific shifts in a gene set in
  malignant cells): LOPO balanced accuracy well above chance;
- negative control with patient-specific signal but **no** subgroup
  signal: LOPO accuracy stays at chance (the permutation p-value is not
  small). The same control with random cell-level splits scores far above
  chance, because the model memorises patients. That is the leakage test:
  if the real pipeline's splitter were swapped for a random one, the
  negative-control test would fail;
- a test that no patient ID occurs in both the training and test indices
  of any fold, and that feature selection receives training patients only.

**Compute:** at 27 × 300 cells × 2,000 genes, one fit takes seconds, one
LOPO run is 27 fits (about a minute), and 1,000 permutations are perhaps
10–20 CPU hours, embarrassingly parallel (an LSF job array: under an hour of
wall time). Tuned (nested CV) instead: 10–20 times that.

### D8. Language boundary: R to Python

| Option | For | Against |
|---|---|---|
| a. Matrix Market counts + `features.tsv` + `barcodes.tsv` + `cells.tsv` (metadata), with a JSON schema giving each metadata column's type (**recommended**) | Plain text, readable by both languages without extra packages; phase 1's `read_mex` already reads it; integers stay exact; easy to checksum | Larger files and slower I/O; types (factor levels, logicals, NA) must be stated explicitly, which the schema does |
| b. `.h5ad` written from R (zellkonverter, which brings a conda environment into the R image through basilisk; or anndataR) | One file; scanpy-native | Heavy dependency in the R image; factor, logical and NA conversions are where round-trip bugs live; harder to inspect |

Recommendation (a) for everything crossing the boundary in either
direction (counts and labels R→Python; Harmony embeddings R→Python as TSV
written with 17 significant digits; predictions and metrics Python→R as
TSV). The Python side may write `.h5ad` internally for its own use.

- **Round-trip test (runs in CI):** R writes a synthetic dataset with every
  awkward case (all-zero gene, all-zero cell, gene names with dots and
  dashes, made-unique symbols, NA, a factor with an unused level, logical
  and integer columns, non-ASCII text); Python reads it and writes it back;
  R reads that and checks the sparse matrix (indices and values),
  dimnames, order and every metadata column are identical, and that the
  SHA-256 of a canonical serialisation matches. Embeddings: max absolute
  difference 0 after the 17-digit round trip.

### D9. Python layout and image (engineering)

- `python/pyproject.toml`: package `mbsc` under `python/src/mbsc/`
  (`io.py` for the boundary, `integration.py`, `classifier.py`,
  `metrics.py`); thin entry points `bin/integration_compare.py` and
  `bin/subgroup_classifier.py`; tests in `python/tests/` (pytest,
  hypothesis); ruff configured in `pyproject.toml`.
- Dependencies pinned in a `uv` lock file (exact versions and hashes):
  Python 3.12, scvi-tools, torch (CPU wheels), anndata, scanpy,
  scikit-learn, scib-metrics, numpy, scipy, pandas, pytest, hypothesis,
  ruff. Alternative: a conda environment (conda-lock). Nothing here needs
  conda's non-Python packages, and uv is faster and smaller.
- Image `containers/py-analysis/Dockerfile` (`python:3.12-slim` + the
  lock), published as `ghcr.io/chialungchuang/mb-scrnaseq-py:<version>` on
  version tags and pinned by digest, like the R image (plan M0 step 6).
  CPU torch keeps it to about 1–2 GB; a CUDA build would be several GB.
- CI: a `python-tests` job (`pytest python/tests`, `ruff check python`)
  and the round-trip test, which needs both images (R writes, Python reads
  and writes back, R checks).

### D10. Wiring

| Option | For | Against |
|---|---|---|
| a. Two subworkflows (`PHASE1`, `PHASE2`) in `subworkflows/local/`, chosen by a required `--phase 1\|2` in `main.nf` (**recommended**) | One provenance, pre-registration gate and `onComplete` code path; the St. Jude profile and resource labels work unchanged; validated like `--cellranger_compat` | `main.nf` grows. Phase 1 code moves into a subworkflow, a refactor after phase 1 has run: its tests must pass unchanged, and phase 1 results stay reproducible from their tag (`v0.2.0`) |
| b. A separate entry script (`phase2.nf`) sharing `modules/` and `conf/` | Phase 1's `main.nf` is never touched after its run | Duplicates the gate, provenance and parameter checks; two entry points to document and test |

Recommendation (a), with phase 2 gated on `phase2_criteria.yaml` exactly
as phase 1 is gated on `acceptance_criteria.yaml`. If moving phase 1 into a
subworkflow after M3 feels too risky, (b) is the safe alternative.

## Compute summary (estimates)

| Step | Jobs | CPU hours | Peak memory per job |
|---|---|---|---|
| Load, QC flags, normalise, PCA, cluster | 1 | < 1 | 16–32 GB (dense CSV load unmeasured) |
| Annotation (scores + SingleR) | 1 | < 1 | < 16 GB |
| inferCNV, no HMM | 28 | 10–30 | 8–32 GB |
| numbat (optional, 8 tumors) | 8 | hundreds | 16–64 GB |
| Programs: scoring / NMF (optional) | 1 / 28 | < 1 / 1–4 | < 16 GB |
| Integration: Harmony + scVI × 5 seeds + metrics | ~11 | 3–5 | < 16 GB |
| Classifier: LOPO + 1,000 permutations + leaky run | array | 10–20 | < 8 GB |

Without numbat, phase 2 is roughly 25–60 CPU hours, well within a normal
LSF allocation; numbat dominates if chosen.

## Test data for M5 step 1

One generator, `tests/data/make_phase2_data.py`, committed with its output
and checked by `tests/shell/test_test_data.sh` like the phase 1
generators. Proposal: 12 synthetic patients (4 per subgroup), about 200
cells each, about 600 genes on 4 synthetic chromosomes, with:

- planted cell types (malignant plus 4 normal lineages) with known marker
  genes;
- planted patient effects on all cells;
- planted CNV gains and losses in each patient's malignant cells, and one
  patient with none;
- planted subgroup signal in malignant cells;
- a second, negative-control output with the same structure and no
  subgroup signal.

Each test above states which planted feature it recovers and is shown to
fail when the corresponding bug is planted (CLAUDE.md rule 3), for
example: splitting by cell instead of patient (negative control fails),
selecting features on all patients (leak test fails), using malignant
cells as the CNV reference (no malignant calls), and normalising per gene
instead of per cell (exact normalisation test fails).

## References

- Riemondy KA, *et al.* Neuro-Oncology 24(2):273–286, 2022 (the dataset).
- Hovestadt V, *et al.* Resolving medulloblastoma cellular architecture by
  single-cell genomics. Nature 572:74–79, 2019.
- Tirosh I, *et al.* Science 352:189–196, 2016 (CNV-based malignant calls).
- Gao R, *et al.* (copykat) Nature Biotechnology 39:599–608, 2021.
- Gao T, *et al.* (numbat) Nature Biotechnology 41:417–426, 2023.
- Gavish A, *et al.* Hallmarks of transcriptional intratumour
  heterogeneity across a thousand tumours. Nature 618:598–606, 2023.
- Aran D, *et al.* (SingleR) Nature Immunology 20:163–172, 2019.
- Choudhary S, Satija R. (SCTransform v2) Genome Biology 23:27, 2022.
- Korsunsky I, *et al.* (Harmony) Nature Methods 16:1289–1296, 2019.
- Lopez R, *et al.* (scVI) Nature Methods 15:1053–1058, 2018.
- Luecken MD, *et al.* (scIB) Nature Methods 19:41–50, 2022.
- Northcott PA, *et al.* Rapid, reliable, and reproducible molecular
  sub-grouping of clinical medulloblastoma samples. Acta Neuropathologica
  123:615–626, 2012.

I wrote these references from memory, without network access to check
them. Volume and page numbers in particular must be checked before this
document is cited anywhere.
