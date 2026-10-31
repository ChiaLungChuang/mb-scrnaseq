published_from <- function(m, label = "MB001") {
  p <- m
  colnames(p) <- paste0(label, "_", colnames(m))
  p
}

test_that("colwise_pearson matches cor() column by column", {
  x <- make_counts(seed = 2); y <- make_counts(seed = 3)
  colnames(y) <- colnames(x)
  r <- colwise_pearson(x, y)
  ref <- vapply(seq_len(ncol(x)), function(j) stats::cor(as.numeric(x[, j]), as.numeric(y[, j])), 0)
  expect_equal(unname(r), ref, tolerance = 1e-10)
})

test_that("colwise_pearson gives NA for a zero-variance column", {
  x <- make_counts(n_cells = 3); y <- x
  x[, 2] <- 0
  expect_true(is.na(colwise_pearson(x, y)[2]))
})

test_that("identical counts give perfect agreement on every metric", {
  m <- make_counts()
  res <- concordance_metrics(as_mex(m), published_from(m), "MB001")$metrics
  expect_equal(res$barcode_recovery, 1)
  expect_equal(res$gene_match_fraction, 1)
  expect_equal(res$umi_total_spearman, 1)
  expect_equal(res$umi_total_median_ratio, 1)
  expect_equal(res$pseudobulk_pearson, 1)
  expect_equal(res$per_cell_pearson_median, 1)
  expect_true(all(evaluate_criteria(res, default_criteria)$pass))
})

test_that("barcodes missing from our matrix lower recovery and fail C1", {
  m <- make_counts(n_cells = 20)
  ours <- as_mex(m[, 1:18])
  res <- concordance_metrics(ours, published_from(m), "MB001")$metrics
  expect_equal(res$barcode_recovery, 18 / 20)
  ev <- evaluate_criteria(res, default_criteria)
  expect_false(ev$pass[ev$criterion == "C1_barcode_recovery"])
})

test_that("whitelist columns with no UMI do not count as recovered", {
  m <- make_counts(n_cells = 20)
  ours <- as_mex(m)
  ours$counts[, 1:4] <- 0                      # 4 published cells got no reads here
  ours <- pad_whitelist(ours, random_barcodes(500))
  res <- concordance_metrics(ours, published_from(m), "MB001")$metrics
  expect_equal(res$barcode_recovery, 16 / 20)
  expect_equal(res$n_matched_cells, 16)
})

test_that("the wrong sample fails C1 and is recorded, not crashed on", {
  m1 <- make_counts(seed = 1); m2 <- make_counts(seed = 9)
  pub <- cbind(published_from(m1, "MB001"), published_from(m2, "MB002"))
  # Sample 1's reads; its raw matrix also has (empty) columns for sample 2's
  # barcodes, because both are on the whitelist.
  ours <- pad_whitelist(as_mex(m1), colnames(m2))
  out <- concordance_metrics(ours, pub, "MB002")
  expect_equal(out$metrics$barcode_recovery, 0)
  expect_true(is.na(out$metrics$pseudobulk_pearson))
  expect_equal(nrow(out$per_cell), 0)
  ev <- evaluate_criteria(out$metrics, default_criteria)
  expect_false(any(ev$pass))
})

test_that("a systematic UMI deficit is caught by the median ratio, not by correlation", {
  m <- make_counts(rate = 20)
  ours <- as_mex(round(m * 0.5))
  res <- concordance_metrics(ours, published_from(m), "MB001")$metrics
  expect_gt(res$umi_total_spearman, 0.95)
  expect_lt(res$umi_total_median_ratio, 0.6)
  ev <- evaluate_criteria(res, default_criteria)
  expect_false(ev$pass[ev$criterion == "C2b_umi_total_median_ratio"])
})

test_that("shuffled cell labels fail per-cell agreement", {
  m <- make_counts(n_cells = 15, seed = 5)
  scrambled <- m
  colnames(scrambled) <- sample(colnames(m))
  res <- concordance_metrics(as_mex(scrambled), published_from(m), "MB001")$metrics
  expect_lt(res$per_cell_pearson_median, 0.9)
})

test_that("only cells with the requested label are compared", {
  m1 <- make_counts(seed = 1); m2 <- make_counts(seed = 9)
  pub <- cbind(published_from(m1, "MB001"), published_from(m2, "MB002"))
  res <- concordance_metrics(as_mex(m1), pub, "MB001")$metrics
  expect_equal(res$n_published_cells, ncol(m1))
  expect_equal(res$barcode_recovery, 1)
  expect_error(concordance_metrics(as_mex(m1), pub, "MB003"), "label 'MB003' not found.*MB001, MB002")
})

test_that("Cell Ranger '-1' suffixes on our barcodes do not break matching", {
  m <- make_counts()
  ours <- as_mex(m)
  colnames(ours$counts) <- paste0(colnames(m), "-1")
  expect_equal(concordance_metrics(ours, published_from(m), "MB001")$metrics$barcode_recovery, 1)
})

test_that("published Ensembl IDs are matched on gene ID with versions stripped", {
  m <- make_counts(n_genes = 10)
  ids <- paste0("ENSG", sprintf("%011d", 1:10))
  pub <- published_from(m)
  rownames(pub) <- paste0(ids, ".7")
  res <- concordance_metrics(as_mex(m, gene_ids = ids), pub, "MB001")$metrics
  expect_equal(res$gene_key, "gene_id")
  expect_equal(res$gene_match_fraction, 1)
})

test_that("genes sharing a symbol are summed before comparison", {
  m <- make_counts(n_genes = 6)
  ours <- as_mex(m)
  ours$features$gene_name[6] <- ours$features$gene_name[5]   # two IDs, one symbol
  pub <- published_from(m)[1:5, ]
  pub[5, ] <- m[5, ] + m[6, ]
  res <- concordance_metrics(ours, pub, "MB001")$metrics
  expect_equal(res$duplicate_symbols, "sum")
  expect_equal(res$pseudobulk_pearson, 1)
})

test_that("symbols made unique Seurat-style (X, X.1) are matched gene by gene", {
  m <- make_counts(n_genes = 6)
  ours <- as_mex(m)
  ours$features$gene_name[5:6] <- "DUP"
  pub <- published_from(m)
  rownames(pub)[5:6] <- c("DUP", "DUP.1")
  res <- concordance_metrics(ours, pub, "MB001")$metrics
  expect_equal(res$duplicate_symbols, "make_unique")
  expect_equal(res$n_shared_genes, 6)
  expect_equal(res$per_cell_pearson_median, 1)
})

test_that("STARsolo's own cell calls are reported when supplied", {
  m <- make_counts(n_cells = 10)
  res <- concordance_metrics(as_mex(m), published_from(m), "MB001", called_barcodes = paste0(colnames(m)[1:8], "-1"))$metrics
  expect_equal(res$published_cells_called_by_starsolo, 0.8)
})

test_that("evaluate_criteria fails a missing or NA metric and needs every threshold", {
  res <- list(barcode_recovery = 1, umi_total_spearman = NA, umi_total_median_ratio = 1,
              pseudobulk_pearson = 1)   # per-cell metric absent
  ev <- evaluate_criteria(res, default_criteria)
  expect_false(ev$pass[ev$criterion == "C2a_umi_total_spearman"])
  expect_false(ev$pass[ev$criterion == "C4_per_cell_pearson_median"])
  bad <- default_criteria; bad$pseudobulk_pearson_min <- NULL
  expect_error(evaluate_criteria(res, bad), "missing the threshold for C3")
})

test_that("criteria boundaries are inclusive", {
  res <- list(barcode_recovery = 0.99, umi_total_spearman = 0.95, umi_total_median_ratio = 1.15,
              pseudobulk_pearson = 0.98, per_cell_pearson_median = 0.90)
  expect_true(all(evaluate_criteria(res, default_criteria)$pass))
})
