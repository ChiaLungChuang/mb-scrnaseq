#!/usr/bin/env Rscript
# Convert the published (GEO) counts CSV into a sparse RDS, once.
#
# Usage: prep_published.R --lib DIR --counts FILE --out-prefix PREFIX [--chunk-rows N]
# Writes PREFIX.rds (genes x cells dgCMatrix) and PREFIX.summary.json
# (orientation, dimensions, cells per sample label).
argv <- commandArgs(trailingOnly = TRUE)
lib <- argv[which(argv == "--lib") + 1]
for (f in list.files(lib, pattern = "\\.R$", full.names = TRUE)) source(f)
a <- parse_args(argv, required = c("lib", "counts", "out_prefix"))
chunk <- as.integer(if (is.null(a$chunk_rows)) 2000 else a$chunk_rows)

m <- read_published_counts(a$counts, chunk_rows = chunk)
ids <- parse_cell_ids(colnames(m))
saveRDS(m, paste0(a$out_prefix, ".rds"))
write_json(list(
  source = basename(a$counts),
  orientation = attr(m, "orientation"),
  n_genes = nrow(m), n_cells = ncol(m),
  gene_key = gene_key_type(rownames(m)),
  cells_per_label = as.list(table(ids$label))
), paste0(a$out_prefix, ".summary.json"))
message(sprintf("published matrix: %d genes x %d cells (%s), %d labels",
                nrow(m), ncol(m), attr(m, "orientation"), length(unique(ids$label))))
