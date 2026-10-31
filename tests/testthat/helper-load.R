# Load the pipeline's R library for testing.
lib_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "lib", "R"))
for (f in list.files(lib_dir, pattern = "\\.R$", full.names = TRUE)) source(f)

# A small genes x cells count matrix with dimnames.
make_counts <- function(n_genes = 30, n_cells = 12, seed = 1, rate = 3) {
  set.seed(seed)
  m <- matrix(stats::rpois(n_genes * n_cells, rate), n_genes, n_cells)
  m[sample(length(m), length(m) %/% 3)] <- 0
  bc <- vapply(seq_len(n_cells), function(i) paste(sample(c("A", "C", "G", "T"), 16, TRUE), collapse = ""), "")
  dimnames(m) <- list(paste0("GENE", seq_len(n_genes)), bc)
  Matrix::Matrix(m, sparse = TRUE)
}

# Wrap a matrix as read_mex() would return it.
as_mex <- function(m, gene_ids = paste0("ENSG", sprintf("%011d", seq_len(nrow(m))))) {
  feats <- data.frame(gene_id = gene_ids, gene_name = rownames(m), stringsAsFactors = FALSE)
  m2 <- m
  rownames(m2) <- gene_ids
  list(counts = m2, features = feats, barcodes = colnames(m))
}

# Add all-zero columns, as STARsolo's raw matrix has for every whitelist
# barcode that no read carried.
pad_whitelist <- function(mex, extra_barcodes) {
  z <- Matrix::Matrix(0, nrow(mex$counts), length(extra_barcodes), sparse = TRUE,
                      dimnames = list(rownames(mex$counts), extra_barcodes))
  mex$counts <- cbind(mex$counts, z)
  mex$barcodes <- colnames(mex$counts)
  mex
}
random_barcodes <- function(n, seed = 99) {
  set.seed(seed)
  vapply(seq_len(n), function(i) paste(sample(c("A", "C", "G", "T"), 16, TRUE), collapse = ""), "")
}

default_criteria <- list(
  barcode_recovery_min = 0.99, umi_total_spearman_min = 0.95,
  umi_total_median_ratio = list(0.85, 1.15), pseudobulk_pearson_min = 0.98,
  per_cell_pearson_median_min = 0.90
)
