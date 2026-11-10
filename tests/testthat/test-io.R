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

test_that("read_published_counts never gives fread a one-line string without a newline", {
  # data.table < 1.15 reads fread(text = x), when x is one string without a
  # newline, as a file name, or as a shell command if it contains a space.
  # The header line and one-line chunks are exactly that. This stand-in
  # refuses such input on any data.table version, so the test fails if the
  # reader hands fread a bare line again.
  real_fread <- data.table::fread
  strict_fread <- function(...) {
    text <- list(...)$text
    if (!is.null(text) && length(text) == 1L && !grepl("[\r\n]", text)) {
      stop("fread(text =) got one line without a newline: a file name or command to data.table < 1.15")
    }
    real_fread(...)
  }
  local_mocked_bindings(fread = strict_fread, .package = "data.table")

  m <- as.matrix(make_counts(n_genes = 7, n_cells = 5))
  colnames(m) <- paste0("MB001_", colnames(m))
  p <- file.path(withr::local_tempdir(), "gxc.csv")
  utils::write.csv(m, p, quote = FALSE)
  # 7 data rows: k = 1 gives seven one-line chunks, k = 2 and k = 3 end on a
  # one-line chunk, k = 100 reads all rows at once.
  for (k in c(1L, 2L, 3L, 100L)) {
    a <- read_published_counts(p, chunk_rows = k)
    expect_equal(unname(as.matrix(a)), unname(m), ignore_attr = TRUE)
    expect_equal(dimnames(a), dimnames(m))
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

test_that("read_published_counts refuses TRUE/FALSE and missing values, in any column of a chunk", {
  # fread reads a column of TRUE/FALSE as logical, which as.matrix() used to
  # turn into 1/0 without a word. A missing value must not become a zero.
  dir <- withr::local_tempdir()
  p <- file.path(dir, "x.csv")
  hdr <- ",MB1_ACGTACGTACGTACGT,MB1_TTTTACGTACGTACGT"
  for (body in list(c("GENE1,1,TRUE", "GENE2,0,FALSE"), c("GENE1,TRUE,2", "GENE2,FALSE,0"))) {
    writeLines(c(hdr, body), p)
    expect_error(read_published_counts(p), "non-numeric or missing", info = body[1])
  }
  for (body in list(c("GENE1,1,", "GENE2,0,3"), c("GENE1,,2", "GENE2,4,0"), c("GENE1,1,NA", "GENE2,0,3"))) {
    writeLines(c(hdr, body), p)
    for (k in c(1, 2)) expect_error(read_published_counts(p, chunk_rows = k), "non-numeric or missing", info = body[1])
  }
})

test_that("rbind_csc stacks row blocks exactly as rbind does", {
  set.seed(5)
  blk <- function(r, n, density) Matrix::rsparsematrix(r, n, density, rand.x = function(k) sample(1:9, k, TRUE))
  cases <- list(
    list(blk(3, 6, 0.4), blk(5, 6, 0.3), blk(1, 6, 0.5)),
    list(blk(4, 6, 0.4), blk(2, 6, 0), blk(3, 6, 0.6)),        # an all-zero block
    list(blk(1, 1, 1), blk(1, 1, 0), blk(1, 1, 1)),            # one column
    list(blk(7, 9, 0.2)),                                       # one block
    list(blk(2, 5, 0), blk(3, 5, 0))                            # nothing at all
  )
  for (k in 1:30) cases[[length(cases) + 1]] <- lapply(sample(1:6, sample(1:8, 1), TRUE), blk, n = 12, density = runif(1, 0, 0.5))
  for (b in cases) {
    want <- methods::as(do.call(rbind, b), "CsparseMatrix")
    got <- rbind_csc(b)
    expect_s4_class(got, "dgCMatrix")
    expect_true(methods::validObject(got))
    expect_identical(got@Dim, want@Dim)
    expect_identical(got@p, want@p)
    expect_identical(got@i, want@i)
    expect_identical(got@x, as.double(want@x))
  }
  expect_error(rbind_csc(list(blk(2, 3, 0.5), blk(2, 4, 0.5))), "number of columns")
})

test_that("read_published_counts works in a fresh R session, before Matrix is loaded", {
  # bin/prep_published.R sources lib/R into a new session; nothing has
  # loaded the Matrix namespace when the first chunk is converted.
  dir <- withr::local_tempdir()
  p <- file.path(dir, "x.csv")
  writeLines(c(",MB1_ACGTACGTACGTACGT,MB1_TTTTACGTACGTACGT", "GENE1,1,0", "GENE2,0,2"), p)
  code <- sprintf(paste0("stopifnot(!'Matrix' %%in%% loadedNamespaces()); ",
                         "for (f in list.files('%s', full.names = TRUE)) source(f); ",
                         "m <- read_published_counts('%s'); cat(class(m), sum(m), '\\n')"), lib_dir, p)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(code)),
                                  stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  expect_match(paste(out, collapse = "\n"), "dgCMatrix 3")
})

test_that("read_published_counts keeps non-integer values exactly", {
  dir <- withr::local_tempdir()
  p <- file.path(dir, "x.csv")
  writeLines(c(",MB1_ACGTACGTACGTACGT,MB1_TTTTACGTACGTACGT", "GENE1,1.5,0", "GENE2,0,2", "GENE3,1e3,0.25"), p)
  for (k in c(1, 2, 3)) {
    got <- read_published_counts(p, chunk_rows = k)
    expect_identical(unname(as.matrix(got)), matrix(c(1.5, 0, 1000, 0, 2, 0.25), 3))
    expect_s4_class(got, "dgCMatrix")
    expect_identical(Matrix::nnzero(got), 4L)   # zeros are not stored
  }
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

test_that("read_mtx reads STARsolo's header layout and agrees with readMM", {
  dir <- withr::local_tempdir()
  p <- file.path(dir, "matrix.mtx")
  # STARsolo writes a bare '%' comment line, integer values, and many empty columns.
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "%", "3 1000 4",
               "1 1 5", "3 1 1", "2 999 2", "1 1000 7"), p)
  m <- read_mtx(p)
  expect_equal(dim(m), c(3L, 1000L))
  expect_equal(sum(m), 15)
  expect_equal(m[2, 999], 2)
  expect_equal(unname(as.matrix(m)), unname(as.matrix(Matrix::readMM(p))))
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "3 5 0"), p)
  expect_equal(sum(read_mtx(p)), 0)
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "3 5 2", "1 1 1"), p)
  expect_error(read_mtx(p), "declares 2 entries but has 1")
  writeLines(c("%%MatrixMarket matrix array real general", "2 2"), p)
  expect_error(read_mtx(p), "not a general coordinate")
})

test_that("write_json writes NA as null, not as the string \"NA\"", {
  p <- withr::local_tempfile(fileext = ".json")
  write_json(list(a = NA_real_, b = 1), p)
  expect_null(jsonlite::read_json(p)$a)
  expect_false(grepl('"NA"', paste(readLines(p), collapse = "")))
})

test_that("log_threads reports the threads data.table will use and the limits in the environment", {
  # Fresh sessions, as the R steps start: data.table reads the limits then.
  run <- function(env) {
    code <- sprintf("for (f in list.files('%s', full.names = TRUE)) source(f); log_threads()", lib_dir)
    out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(code)),
                   env = env, stdout = TRUE, stderr = TRUE)
    grep("^threads:", out, value = TRUE)
  }
  expect_identical(run(c("OMP_NUM_THREADS=1", "OMP_THREAD_LIMIT=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1", "R_DATATABLE_NUM_THREADS=1")),
                   "threads: data.table 1, OMP_NUM_THREADS=1, OMP_THREAD_LIMIT=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1")
  # The data.table figure is data.table's own count, not copied from the variables.
  two <- min(2L, parallel::detectCores())
  expect_identical(run(c(sprintf("OMP_THREAD_LIMIT=%d", two), sprintf("R_DATATABLE_NUM_THREADS=%d", two), "OMP_NUM_THREADS=7")),
                   sprintf("threads: data.table %d, OMP_NUM_THREADS=7, OMP_THREAD_LIMIT=%d, OPENBLAS_NUM_THREADS=unset, MKL_NUM_THREADS=unset", two, two))
})
