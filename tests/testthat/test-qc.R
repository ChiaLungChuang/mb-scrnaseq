test_that("qc_metrics computes totals, detected genes and mito percentage", {
  m <- Matrix::Matrix(matrix(c(10, 0, 5,  0, 0, 0,  2, 2, 0), 3, 3), sparse = TRUE)
  dimnames(m) <- list(c("GENE1", "MT-CO1", "mt-Nd1"), c("a", "b", "c"))
  q <- qc_metrics(m, rownames(m))
  expect_equal(q$total_umi, c(15, 0, 4))
  expect_equal(q$detected_genes, c(2, 0, 2))
  expect_equal(q$pct_mito, c(100 * 5 / 15, NA, 50))   # mito match is case-insensitive
})

test_that("qc_flag flags the deliberately poor cells and records thresholds", {
  set.seed(11)
  good <- data.frame(barcode = paste0("g", 1:50), total_umi = round(stats::rlnorm(50, log(3000), 0.2)),
                     detected_genes = round(stats::rlnorm(50, log(1200), 0.15)), pct_mito = stats::runif(50, 2, 6))
  bad <- data.frame(barcode = paste0("b", 1:3), total_umi = c(150, 3000, 3100),
                    detected_genes = c(90, 1100, 1250), pct_mito = c(4, 45, 60))
  res <- qc_flag(rbind(good, bad), nmads = 3)
  flagged <- res$cells$barcode[res$cells$flag_any]
  expect_true(all(c("b1", "b2", "b3") %in% flagged))
  expect_lte(sum(startsWith(flagged, "g")), 2)
  expect_true(res$cells$flag_low_umi[res$cells$barcode == "b1"])
  expect_true(res$cells$flag_high_mito[res$cells$barcode == "b3"])
  expect_true(is.numeric(res$thresholds$total_umi$lower))
  expect_equal(res$thresholds$n_cells, 53)
  expect_false(res$thresholds$pct_mito$mad_zero)
})

test_that("qc_flag reports a zero MAD instead of hiding it", {
  df <- data.frame(barcode = letters[1:6], total_umi = c(100, 100, 100, 100, 100, 120),
                   detected_genes = c(50, 52, 49, 51, 50, 48), pct_mito = c(1, 2, 3, 2, 1, 2))
  res <- qc_flag(df)
  expect_true(res$thresholds$total_umi$mad_zero)
})
