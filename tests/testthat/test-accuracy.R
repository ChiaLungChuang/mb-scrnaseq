# Accuracy tests: each function in lib/R is checked against an independent
# implementation (written here in plain dense base R, or taken from the
# reference package) on many random inputs. The unit tests in the other
# files check chosen cases; these check that the numbers are right in
# general, including the awkward inputs real data produces: empty whitelist
# columns, "-1" barcode suffixes, duplicated gene symbols, unmatched genes,
# cells with no counts.

n_trials <- 25

rand_barcodes <- function(n) {
  unique(vapply(seq_len(n * 2), function(i) paste(sample(c("A", "C", "G", "T"), 16, TRUE), collapse = ""), ""))[seq_len(n)]
}

# A random "published" matrix: several labels, gene symbols, sparse counts.
rand_published <- function(n_genes, labels, cells_per_label) {
  bcs <- rand_barcodes(length(labels) * cells_per_label)
  ids <- paste0(rep(labels, each = cells_per_label), "_", bcs)
  m <- matrix(stats::rpois(n_genes * length(ids), stats::runif(1, 0.5, 4)), n_genes, length(ids))
  m[sample(length(m), length(m) %/% 2)] <- 0
  dimnames(m) <- list(paste0("SYM", seq_len(n_genes)), ids)
  m
}

# Independent reference for concordance_metrics(): dense matrices, base R
# only, written from the criteria's definitions rather than from the code.
naive_concordance <- function(ours_dense, ours_ids, ours_names, pub_dense, label) {
  lab <- sub("_[ACGTN]{16}$", "", colnames(pub_dense))
  bc <- sub("^.*_", "", colnames(pub_dense))
  pub <- pub_dense[, lab == label, drop = FALSE]
  colnames(pub) <- bc[lab == label]

  colnames(ours_dense) <- sub("-1$", "", colnames(ours_dense))
  ours_dense <- ours_dense[, colSums(ours_dense) > 0, drop = FALSE]
  uniq <- make.unique(ours_names)
  if (any(uniq[uniq != ours_names] %in% rownames(pub))) {
    rownames(ours_dense) <- uniq
  } else {
    ours_dense <- rowsum(ours_dense, ours_names, reorder = TRUE)
  }
  shared <- intersect(rownames(pub), rownames(ours_dense))
  matched <- intersect(colnames(pub), colnames(ours_dense))
  out <- list(n_published_cells = ncol(pub), n_matched_cells = length(matched),
              barcode_recovery = length(matched) / ncol(pub),
              n_shared_genes = length(shared), gene_match_fraction = length(shared) / nrow(pub))
  if (length(matched) < 3) return(out)
  a <- ours_dense[shared, matched, drop = FALSE]
  b <- pub[shared, matched, drop = FALSE]
  ta <- colSums(a); tb <- colSums(b)
  pa <- log1p(rowSums(a) / sum(a) * 1e6); pb <- log1p(rowSums(b) / sum(b) * 1e6)
  r <- vapply(seq_along(matched), function(j) {
    x <- log1p(a[, j]); y <- log1p(b[, j])
    if (stats::sd(x) == 0 || stats::sd(y) == 0) NA_real_ else stats::cor(x, y)
  }, 0)
  c(out, list(umi_total_spearman = stats::cor(ta, tb, method = "spearman"),
              umi_total_median_ratio = stats::median(ta[tb > 0] / tb[tb > 0]),
              pseudobulk_pearson = stats::cor(pa, pb),
              per_cell_pearson_median = stats::median(r, na.rm = TRUE),
              per_cell_pearson_na = sum(is.na(r))))
}

test_that("concordance_metrics agrees with an independent dense implementation on random data", {
  set.seed(20261031)
  for (trial in seq_len(n_trials)) {
    labels <- c("MB1", "MB2", "MB3")
    pub <- rand_published(n_genes = sample(15:40, 1), labels, cells_per_label = sample(5:20, 1))
    target <- sample(labels, 1)
    in_target <- startsWith(colnames(pub), paste0(target, "_"))
    target_bc <- sub("^.*_", "", colnames(pub)[in_target])

    # Our side: the target's cells with noise, some cells lost, extra genes
    # the published matrix lacks, some published genes missing, empty
    # whitelist columns, other samples' barcodes as empty columns, "-1".
    keep_cells <- target_bc[stats::runif(length(target_bc)) > stats::runif(1, 0, 0.3)]
    pub_rows <- sample(rownames(pub), round(nrow(pub) * stats::runif(1, 0.7, 1)))
    base <- pub[pub_rows, in_target, drop = FALSE][, match(keep_cells, target_bc), drop = FALSE]
    noisy <- matrix(stats::rpois(length(base), base * stats::runif(1, 0.7, 1.3)), nrow(base))
    extra <- matrix(stats::rpois(3 * ncol(base), 1), 3, ncol(base))
    ours_names <- c(pub_rows, paste0("EXTRA", 1:3))
    ours_dense <- rbind(noisy, extra)

    mode <- sample(c("none", "sum", "make_unique"), 1)
    if (mode != "none") {      # one published gene appears twice on our side
      dup <- pub_rows[1]
      ours_dense <- rbind(ours_dense, stats::rpois(ncol(ours_dense), 2))
      ours_names <- c(ours_names, dup)
      if (mode == "make_unique") {
        pub <- rbind(pub, stats::rpois(ncol(pub), 2))
        rownames(pub)[nrow(pub)] <- paste0(dup, ".1")
      }
    }
    # STARsolo's raw matrix has a column for every whitelist barcode, so lost
    # cells and other samples' cells are present as empty columns.
    lost_cells <- setdiff(target_bc, keep_cells)
    pad <- c(rand_barcodes(30), sub("^.*_", "", colnames(pub)[!in_target]), lost_cells)
    ours_dense <- cbind(ours_dense, matrix(0, nrow(ours_dense), length(pad)))
    colnames(ours_dense) <- paste0(c(keep_cells, pad), "-1")

    ours <- list(counts = Matrix::Matrix(ours_dense, sparse = TRUE,
                                         dimnames = list(paste0("ENSG", seq_along(ours_names)), colnames(ours_dense))),
                 features = data.frame(gene_id = paste0("ENSG", seq_along(ours_names)), gene_name = ours_names),
                 barcodes = colnames(ours_dense))
    got <- concordance_metrics(ours, Matrix::Matrix(pub, sparse = TRUE), target)$metrics
    want <- naive_concordance(ours_dense, colnames(ours_dense), ours_names, pub, target)
    for (k in names(want)) {
      expect_equal(got[[k]], want[[k]], tolerance = 1e-10,
                   info = sprintf("trial %d (%s duplicates), metric %s", trial, mode, k))
    }
    if (mode != "none") expect_equal(got$duplicate_symbols, if (mode == "sum") "sum" else "make_unique")
  }
})

test_that("evaluate_criteria agrees with a direct reading of the thresholds", {
  set.seed(7)
  crit <- default_criteria
  for (trial in seq_len(200)) {
    m <- list(barcode_recovery = stats::runif(1, 0.95, 1), umi_total_spearman = stats::runif(1, 0.9, 1),
              umi_total_median_ratio = stats::runif(1, 0.8, 1.2), pseudobulk_pearson = stats::runif(1, 0.96, 1),
              per_cell_pearson_median = stats::runif(1, 0.85, 1))
    want <- c(m$barcode_recovery >= 0.99, m$umi_total_spearman >= 0.95,
              m$umi_total_median_ratio >= 0.85 && m$umi_total_median_ratio <= 1.15,
              m$pseudobulk_pearson >= 0.98, m$per_cell_pearson_median >= 0.90)
    expect_identical(evaluate_criteria(m, crit)$pass, want)
  }
})

test_that("colwise_pearson equals cor() on random sparse matrices, including constant columns", {
  set.seed(3)
  for (trial in seq_len(n_trials)) {
    n <- sample(5:60, 1); k <- sample(2:30, 1)
    x <- matrix(stats::rpois(n * k, 2), n, k); y <- x + matrix(stats::rpois(n * k, 1), n, k)
    x[, 1] <- 0; y[, 2] <- 3                          # zero-variance columns -> NA
    got <- colwise_pearson(Matrix::Matrix(x, sparse = TRUE), Matrix::Matrix(y, sparse = TRUE))
    want <- suppressWarnings(vapply(seq_len(k), function(j) stats::cor(x[, j], y[, j]), 0))
    expect_equal(unname(got), want, tolerance = 1e-10)
  }
})

test_that("qc_metrics agrees with scuttle::perCellQCMetrics", {
  set.seed(5)
  for (trial in seq_len(n_trials)) {
    n_genes <- sample(20:80, 1); n_cells <- sample(5:40, 1)
    m <- matrix(stats::rpois(n_genes * n_cells, stats::runif(1, 0.2, 3)), n_genes, n_cells)
    genes <- paste0("G", seq_len(n_genes))
    mt <- sample(n_genes, 4)
    genes[mt] <- paste0(sample(c("MT-", "mt-"), 4, TRUE), "X", seq_along(mt))
    m[, 1] <- pmax(m[, 1], 1)                       # at least one cell with counts everywhere
    dimnames(m) <- list(genes, paste0("c", seq_len(n_cells)))
    sm <- Matrix::Matrix(m, sparse = TRUE)
    got <- qc_metrics(sm, genes)
    ref <- scuttle::perCellQCMetrics(sm, subsets = list(Mito = grepl("^mt-", genes, ignore.case = TRUE)))
    expect_equal(got$total_umi, as.numeric(ref$sum))
    expect_equal(got$detected_genes, as.numeric(ref$detected))
    ok <- ref$sum > 0
    expect_equal(got$pct_mito[ok], as.numeric(ref$subsets_Mito_percent[ok]), tolerance = 1e-12)
  }
})

test_that("qc_flag thresholds and flags equal median -/+ n MADs computed by hand", {
  set.seed(9)
  for (trial in seq_len(n_trials)) {
    n <- sample(30:200, 1)
    df <- data.frame(barcode = paste0("c", seq_len(n)),
                     total_umi = round(stats::rlnorm(n, log(2000), stats::runif(1, 0.1, 0.6))) + 1,
                     detected_genes = round(stats::rlnorm(n, log(800), stats::runif(1, 0.1, 0.5))) + 1,
                     pct_mito = stats::rgamma(n, 2, 0.5))
    nmads <- sample(c(2, 3, 4), 1)
    res <- qc_flag(df, nmads = nmads)
    lo <- function(x) exp(stats::median(log(x)) - nmads * stats::mad(log(x)))
    hi <- function(x) stats::median(x) + nmads * stats::mad(x)
    expect_equal(res$thresholds$total_umi$lower, lo(df$total_umi), tolerance = 1e-10)
    expect_equal(res$thresholds$detected_genes$lower, lo(df$detected_genes), tolerance = 1e-10)
    expect_equal(res$thresholds$pct_mito$higher, hi(df$pct_mito), tolerance = 1e-10)
    expect_identical(res$cells$flag_low_umi, df$total_umi < lo(df$total_umi))
    expect_identical(res$cells$flag_low_genes, df$detected_genes < lo(df$detected_genes))
    expect_identical(res$cells$flag_high_mito, df$pct_mito > hi(df$pct_mito))
    expect_identical(res$cells$flag_any, res$cells$flag_low_umi | res$cells$flag_low_genes | res$cells$flag_high_mito)
    expect_equal(res$thresholds$n_flagged, sum(res$cells$flag_any))
  }
})

test_that("read_published_counts reproduces random matrices exactly, any layout and chunk size", {
  set.seed(11)
  for (trial in seq_len(n_trials)) {
    pub <- rand_published(n_genes = sample(3:50, 1), c("L1", "L22"), cells_per_label = sample(2:15, 1))
    cells_x_genes <- runif(1) < 0.5
    gz <- runif(1) < 0.5
    p <- withr::local_tempfile(fileext = if (gz) ".csv.gz" else ".csv")
    tab <- if (cells_x_genes) t(pub) else pub
    con <- if (gz) gzfile(p, "w") else file(p, "w")
    writeLines(paste0(",", paste(colnames(tab), collapse = ",")), con)
    for (i in seq_len(nrow(tab))) writeLines(paste0(rownames(tab)[i], ",", paste(tab[i, ], collapse = ",")), con)
    close(con)
    got <- read_published_counts(p, chunk_rows = sample(c(1, 2, 7, 1000), 1))
    expect_equal(attr(got, "orientation"), if (cells_x_genes) "cells_x_genes" else "genes_x_cells")
    expect_equal(dimnames(got), dimnames(pub))
    expect_equal(unname(as.matrix(got)), unname(pub * 1.0))
  }
})

test_that("read_mtx equals Matrix::readMM on random sparse matrices", {
  set.seed(13)
  for (trial in seq_len(n_trials)) {
    nr <- sample(1:40, 1); nc <- sample(1:300, 1)
    m <- Matrix::rsparsematrix(nr, nc, density = stats::runif(1, 0.001, 0.3),
                               rand.x = function(n) as.numeric(sample(1:50, n, TRUE)))
    p <- withr::local_tempfile(fileext = ".mtx")
    Matrix::writeMM(m, p)
    got <- read_mtx(p)
    expect_equal(dim(got), dim(m))
    expect_equal(unname(as.matrix(got)), unname(as.matrix(Matrix::readMM(p))))
  }
})
