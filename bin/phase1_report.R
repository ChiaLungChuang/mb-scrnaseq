#!/usr/bin/env Rscript
# Write the phase 1 results report (docs/results/phase1/README.md) from a
# finished run's output directory: verdict per sample and criterion,
# pre-registration status, commit and hashes, QC flags and compute.
#
# Usage: phase1_report.R --lib DIR --outdir RUN_OUTDIR --out FILE [--trace FILE]
#   --trace  the run's trace if not RUN_OUTDIR/pipeline_info/trace.tsv
# It stops, writing nothing, if the run's outputs are incomplete or do not
# belong together.
argv <- commandArgs(trailingOnly = TRUE)
lib <- argv[which(argv == "--lib") + 1]
if (!length(lib) || is.na(lib)) stop("--lib DIR is required (the repository's lib/R)")
for (f in list.files(lib, pattern = "\\.R$", full.names = TRUE)) source(f)
a <- parse_args(argv, required = c("lib", "outdir", "out"))
md <- phase1_report(a$outdir, trace = a$trace)
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
writeLines(md, a$out)
message("wrote ", a$out)
