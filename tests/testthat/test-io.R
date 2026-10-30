test_that("parse_args builds a named list and enforces pairs and required keys", {
  a <- parse_args(c("--sample-id", "s1", "--nmads", "3"), required = c("sample_id"))
  expect_equal(a$sample_id, "s1")
  expect_equal(a$nmads, "3")
  expect_error(parse_args(c("--sample-id")), "pairs")
  expect_error(parse_args(c("sample", "x")), "expected '--key'")
  expect_error(parse_args(c("--a", "1"), required = "b"), "missing required argument")
  expect_error(parse_args(c("--a", "1", "--a", "2")), "duplicated")
})

test_that("parse_cell_ids splits label and barcode in either order", {
  bc <- "ACGTACGTACGTACGT"
  p <- parse_cell_ids(c(paste0("MB001_", bc), paste0(bc, "-1_MB002"), paste0("1166:", bc), paste0(bc, "-1")))
  expect_equal(p$label, c("MB001", "MB002", "1166", ""))
  expect_equal(unique(p$barcode), bc)
})

test_that("parse_cell_ids refuses IDs it would have to guess about", {
  expect_error(parse_cell_ids("MB001_cell7"), "0 barcode-like tokens")
  expect_error(parse_cell_ids("ACGTACGTACGTACGT_TTTTACGTACGTACGT"), "2 barcode-like tokens")
})

test_that("read_published_counts detects orientation, handles gzip and any chunk size", {
  m <- as.matrix(make_counts(n_genes = 7, n_cells = 5))
  colnames(m) <- paste0("MB001_", colnames(m))
  dir <- withr::local_tempdir()

  gxc <- file.path(dir, "gxc.csv")
  utils::write.csv(m, gxc, quote = FALSE)
  cxg <- file.path(dir, "cxg.csv.gz")
  con <- gzfile(cxg, "w"); utils::write.csv(t(m), con, quote = FALSE); close(con)

  for (k in c(1L, 2L, 3L, 100L)) {
    a <- read_published_counts(gxc, chunk_rows = k)
    expect_equal(attr(a, "orientation"), "genes_x_cells")
    expect_equal(unname(as.matrix(a)), unname(m), ignore_attr = TRUE)
    expect_equal(dimnames(a), dimnames(m))
    b <- read_published_counts(cxg, chunk_rows = k)
    expect_equal(attr(b, "orientation"), "cells_x_genes")
    expect_equal(dimnames(b), dimnames(m))
    expect_equal(unname(as.matrix(b)), unname(m), ignore_attr = TRUE)
  }
})

test_that("read_published_counts refuses a file with no recognisable cell IDs", {
  dir <- withr::local_tempdir()
  p <- file.path(dir, "x.csv")
  writeLines(c(",cellA,cellB", "GENE1,1,2", "GENE2,0,3"), p)
  expect_error(read_published_counts(p), "could not find cell IDs")
})

test_that("read_published_counts rejects non-numeric values", {
  dir <- withr::local_tempdir()
  p <- file.path(dir, "x.csv")
  writeLines(c(",MB1_ACGTACGTACGTACGT", "GENE1,1", "GENE2,oops"), p)
  expect_error(read_published_counts(p), "non-numeric")
})

test_that("read_mex round-trips a STARsolo-style directory", {
  m <- make_counts(n_genes = 6, n_cells = 4)
  dir <- withr::local_tempdir()
  Matrix::writeMM(m, file.path(dir, "matrix.mtx"))
  writeLines(colnames(m), file.path(dir, "barcodes.tsv"))
  writeLines(paste(paste0("ENSG", 1:6), rownames(m), "Gene Expression", sep = "\t"), file.path(dir, "features.tsv"))
  x <- read_mex(dir)
  expect_equal(dim(x$counts), dim(m))
  expect_equal(x$features$gene_name, rownames(m))
  expect_equal(unname(as.matrix(x$counts)), unname(as.matrix(m)))
  writeLines(colnames(m)[1:3], file.path(dir, "barcodes.tsv"))
  expect_error(read_mex(dir), "barcodes file has 3")
})
