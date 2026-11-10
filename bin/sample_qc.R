#!/usr/bin/env Rscript
# Per-sample QC on STARsolo's filtered cells: metrics, MAD-based flags, and
# the thresholds that produced them. Flags only; nothing is removed here.
#
# Usage: sample_qc.R --lib DIR --sample ID --filtered-dir DIR --nmads N --out-prefix PREFIX
argv <- commandArgs(trailingOnly = TRUE)
lib <- argv[which(argv == "--lib") + 1]
for (f in list.files(lib, pattern = "\\.R$", full.names = TRUE)) source(f)
log_threads()
a <- parse_args(argv, required = c("lib", "sample", "filtered_dir", "nmads", "out_prefix"))

x <- read_mex(a$filtered_dir)
if (ncol(x$counts) < 3) stop("only ", ncol(x$counts), " cells passed STARsolo's filter for ", a$sample)
q <- qc_metrics(x$counts, x$features$gene_name)
res <- qc_flag(q, nmads = as.numeric(a$nmads))
res$cells$sample <- a$sample

con <- gzfile(paste0(a$out_prefix, ".qc_cells.tsv.gz"), "w")
utils::write.table(res$cells, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
write_json(c(list(sample = a$sample), res$thresholds), paste0(a$out_prefix, ".qc_thresholds.json"))
message(sprintf("%s: %d cells, %d flagged", a$sample, res$thresholds$n_cells, res$thresholds$n_flagged))
