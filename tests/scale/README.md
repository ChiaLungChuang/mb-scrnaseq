# Full-size rehearsal of the R steps

The R steps (`PREP_PUBLISHED`, `SAMPLE_QC`, `CONCORDANCE`) are tested in CI
on matrices a few hundred cells wide. This rehearsal runs the same R code on
synthetic inputs at the real dimensions, to size their memory and time
requests before the first real run. It is not part of CI (the inputs are
several GB and take minutes to make) and nothing it makes is committed.
No real data is involved.

```bash
Rscript tests/scale/make_scale_data.R --out /scratch/scale   # ~4 min; 0.5 GB (2.7 GB while writing)
bash tests/scale/rehearse_r_steps.sh /scratch/scale          # ~4 min
# a denser case: twice the cells, 5x the ambient barcodes, ~3,500 genes per cell
Rscript tests/scale/make_scale_data.R --out /scratch/dense --cells 10000 --ambient 1500000 --umis 12000
bash tests/scale/rehearse_r_steps.sh /scratch/dense
```

`make_scale_data.R` writes:

- a STARsolo raw matrix of 33,538 genes x 6,794,880 barcodes (the v3
  whitelist) with the called cells plus ambient barcodes of 1-6 UMIs, and
  the matching filtered matrix;
- a dense published-counts CSV.gz of 26,841 genes x 39,946 cells, with the
  30 labels and per-label cell counts of the real matrix (from
  `docs/proposals/gse155446_runs.tsv`, i.e. from cell IDs only). The
  rehearsal label (1167, the largest, 3,073 cells) holds a subset of the
  called cells with the same counts, so concordance does its full work.

`rehearse_r_steps.sh` runs each step's `bin/` script as the pipeline does
and records peak resident memory and wall time (`getrusage` of the child,
the figures `/usr/bin/time -v` reports), with data.table on one thread, as
a one-CPU task has.

## Results (3 Nov 2026)

R 4.4.3, Matrix 1.7.6, data.table 1.18.6 (the pipeline's image has R 4.4.3
from a 2025 package snapshot), one thread, Intel Xeon 2.8 GHz.

| Inputs | typical | dense |
|---|---|---|
| called cells, median genes per cell | 5,000, 2,059 | 10,000, 3,520 |
| raw matrix entries (ambient barcodes) | 11.6 M (300,000) | 40.5 M (1,500,000) |
| published non-zeros (density) | 78.3 M (7.3%) | 131.6 M (12.3%) |

| Step | typical: peak, time | dense: peak, time | request |
|---|---|---|---|
| `PREP_PUBLISHED`, before | 5.4 GB, 193 s | 6.9 GB, 228 s | |
| `PREP_PUBLISHED`, after | 3.3 GB, 144 s | 4.8 GB, 179 s | 12 GB, 1 h (`process_r`; was 16 GB, 2 h) |
| `SAMPLE_QC` | 0.8 GB, 8-16 s | 1.9 GB, 13-14 s | 4 GB, 1 h (`process_low`; was 16 GB, 2 h) |
| `CONCORDANCE` | 3.4 GB, 39-54 s | 5.8 GB, 69-83 s | 12 GB, 1 h (`process_r`; was 16 GB, 2 h) |

Requests double on each retry (`conf/base.config`), so a first attempt that
runs out still finishes on the second.

"Before" is the published-matrix reader as it was; "after" is the same
reader without two avoidable copies (`lib/R/io.R`):

- Each chunk of 2,000 rows went through a dense double matrix (640 MB for
  39,946 cells) on its way to sparse. Chunks are now converted column by
  column from the parsed table (`chunk_to_sparse()`).
- The chunks were stacked with `do.call(rbind, ...)`, which held about 2.5
  times the result at its peak (6.0 GB of the 6.9 GB on the dense case).
  `rbind_csc()` writes each chunk's entries straight into the result.

The output is identical to the old reader's on the dense case (every slot,
the dimnames and the orientation).

Not measured here: the real data (the published matrix's density is
unknown until it is read, which waits for the criteria to be locked),
HPCF's CPUs and file system, and the image's exact package versions. Check
`peak_rss` and `realtime` in `pipeline_info/trace.tsv` after the first
real run.
