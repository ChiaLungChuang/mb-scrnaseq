# Agreement between our re-quantified counts and the published counts.
#
# The question these functions answer is narrow: given the same reads, does
# our STARsolo run reproduce the matrix the authors deposited? Disagreement
# here means a wrong chemistry, whitelist, reference or sample mapping, and
# has to be resolved before anything downstream is worth running.

#' Sum rows that share a key (e.g. several gene IDs with one symbol).
collapse_rows <- function(m, keys) {
  if (length(keys) != nrow(m)) stop("keys length ", length(keys), " != nrow ", nrow(m))
  f <- factor(keys)
  out <- Matrix::fac2sparse(f) %*% m
  rownames(out) <- levels(f)
  methods::as(out, "CsparseMatrix")
}

#' Decide whether published rows are Ensembl IDs or gene symbols.
gene_key_type <- function(published_genes, min_frac = 0.9) {
  if (mean(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", published_genes)) >= min_frac) "gene_id" else "gene_name"
}

#' Keys for matching our genes to published gene symbols.
#'
#' Seurat's Read10X() makes duplicated symbols unique in features.tsv order
#' (RGS5, RGS5.1), so a published matrix built that way keeps the genes
#' separate. If the published names show that pattern we do the same (STAR,
#' like Cell Ranger, lists genes in GTF order); otherwise genes that share a
#' symbol are summed. Decided from gene names only, never from counts.
symbol_keys <- function(our_names, published_genes) {
  uniq <- make.unique(our_names)
  renamed <- uniq[uniq != our_names]
  if (length(renamed) && any(renamed %in% published_genes)) {
    list(keys = uniq, mode = "make_unique")
  } else {
    list(keys = our_names, mode = "sum")
  }
}

#' Per-column Pearson correlation of two equally shaped sparse matrices.
#'
#' Exact, computed from sums so it never densifies the matrices. Columns with
#' zero variance give NA, as cor() would.
colwise_pearson <- function(x, y) {
  stopifnot(identical(dim(x), dim(y)))
  n <- nrow(x)
  sx <- Matrix::colSums(x); sy <- Matrix::colSums(y)
  sxx <- Matrix::colSums(x * x); syy <- Matrix::colSums(y * y)
  sxy <- Matrix::colSums(x * y)
  num <- sxy - sx * sy / n
  den <- sqrt((sxx - sx^2 / n) * (syy - sy^2 / n))
  r <- num / den
  r[!is.finite(r)] <- NA_real_
  r
}

log1p_sparse <- function(m) {
  m <- methods::as(m, "CsparseMatrix")
  m@x <- log1p(m@x)
  m
}

#' Compare one sample's STARsolo raw matrix with the published matrix.
#'
#' A published barcode counts as recovered only if it has at least one UMI
#' in our raw matrix. STARsolo's raw matrix has a column for every whitelist
#' barcode (6.8M for v3), so mere presence of a column says nothing. With
#' fewer than 3 recovered cells the other metrics are NA (and fail): a wrong
#' sample is a result to record, not a crash.
#'
#' @param ours     list from read_mex() on Solo.out/Gene/raw
#' @param published genes x cells matrix from read_published_counts()
#' @param label    how this sample is named in the published cell IDs
#' @param called_barcodes optional: barcodes STARsolo's own cell filter kept
#' @return list of metrics plus a per-cell data.frame
concordance_metrics <- function(ours, published, label, called_barcodes = NULL) {
  ids <- parse_cell_ids(colnames(published))
  keep <- ids$label == label
  if (!any(keep)) {
    labs <- sort(unique(ids$label))
    stop(sprintf("label '%s' not found in published cell IDs; labels present: %s",
                 label, paste(utils::head(labs, 40), collapse = ", ")))
  }
  pub <- published[, keep, drop = FALSE]
  colnames(pub) <- ids$barcode[keep]
  if (anyDuplicated(colnames(pub))) stop("duplicated barcodes within published label '", label, "'")

  key <- gene_key_type(rownames(published))
  ours_counts <- ours$counts
  colnames(ours_counts) <- sub("-[0-9]+$", "", colnames(ours_counts))
  # Drop empty whitelist columns first: recovery needs >= 1 UMI, and the
  # raw matrix shrinks from millions of columns to the barcodes seen.
  ours_counts <- ours_counts[, Matrix::colSums(ours_counts) > 0, drop = FALSE]
  if (key == "gene_id") {
    ours_keys <- sub("\\.[0-9]+$", "", ours$features$gene_id)
    symbol_mode <- NA_character_
  } else {
    sk <- symbol_keys(ours$features$gene_name, rownames(pub))
    ours_keys <- sk$keys
    symbol_mode <- sk$mode
  }
  ours_g <- collapse_rows(ours_counts, ours_keys)
  pub_keys <- if (key == "gene_id") sub("\\.[0-9]+$", "", rownames(pub)) else rownames(pub)
  pub_g <- collapse_rows(pub, pub_keys)

  shared_genes <- intersect(rownames(pub_g), rownames(ours_g))
  matched <- intersect(colnames(pub_g), colnames(ours_g))
  if (length(shared_genes) < 2) stop("fewer than 2 genes shared between our features and the published matrix")

  metrics <- list(
    label = label,
    gene_key = key,
    duplicate_symbols = symbol_mode,
    n_published_cells = ncol(pub_g),
    n_matched_cells = length(matched),
    barcode_recovery = length(matched) / ncol(pub_g),
    n_published_genes = nrow(pub_g),
    n_shared_genes = length(shared_genes),
    gene_match_fraction = length(shared_genes) / nrow(pub_g),
    umi_total_spearman = NA_real_,
    umi_total_median_ratio = NA_real_,
    pseudobulk_pearson = NA_real_,
    per_cell_pearson_median = NA_real_,
    per_cell_pearson_na = NA_integer_
  )
  per_cell <- data.frame(barcode = character(), umi_ours = numeric(), umi_published = numeric(),
                         pearson_log1p = numeric())

  if (length(matched) >= 3) {
    a <- ours_g[shared_genes, matched, drop = FALSE]
    b <- pub_g[shared_genes, matched, drop = FALSE]
    tot_a <- Matrix::colSums(a); tot_b <- Matrix::colSums(b)
    ratio <- tot_a[tot_b > 0] / tot_b[tot_b > 0]

    cpm <- function(v) v / sum(v) * 1e6
    pb_a <- log1p(cpm(Matrix::rowSums(a))); pb_b <- log1p(cpm(Matrix::rowSums(b)))
    per_cell_r <- colwise_pearson(log1p_sparse(a), log1p_sparse(b))

    metrics$umi_total_spearman <- suppressWarnings(stats::cor(tot_a, tot_b, method = "spearman"))
    metrics$umi_total_median_ratio <- stats::median(ratio)
    metrics$pseudobulk_pearson <- stats::cor(pb_a, pb_b)
    metrics$per_cell_pearson_median <- stats::median(per_cell_r, na.rm = TRUE)
    metrics$per_cell_pearson_na <- sum(is.na(per_cell_r))
    per_cell <- data.frame(barcode = matched, umi_ours = tot_a, umi_published = tot_b,
                           pearson_log1p = per_cell_r, row.names = NULL)
  }
  if (!is.null(called_barcodes)) {
    called <- sub("-[0-9]+$", "", called_barcodes)
    metrics$published_cells_called_by_starsolo <- mean(colnames(pub_g) %in% called)
    metrics$n_starsolo_called_cells <- length(called)
  }
  list(metrics = metrics, per_cell = per_cell)
}

#' Check metrics against the pre-registered thresholds.
#'
#' Thresholds come from the phase1_concordance block of the criteria file.
#' Returns one row per criterion with the observed value, the rule and
#' whether it passed. A missing or NA metric fails; it never passes by
#' default.
evaluate_criteria <- function(metrics, crit) {
  rules <- list(
    list(id = "C1_barcode_recovery", metric = "barcode_recovery", op = ">=", value = crit$barcode_recovery_min),
    list(id = "C2a_umi_total_spearman", metric = "umi_total_spearman", op = ">=", value = crit$umi_total_spearman_min),
    list(id = "C2b_umi_total_median_ratio", metric = "umi_total_median_ratio", op = "in", value = crit$umi_total_median_ratio),
    list(id = "C3_pseudobulk_pearson", metric = "pseudobulk_pearson", op = ">=", value = crit$pseudobulk_pearson_min),
    list(id = "C4_per_cell_pearson_median", metric = "per_cell_pearson_median", op = ">=", value = crit$per_cell_pearson_median_min)
  )
  rows <- lapply(rules, function(r) {
    if (is.null(r$value)) stop("criteria file is missing the threshold for ", r$id)
    obs <- metrics[[r$metric]]
    obs <- if (is.null(obs)) NA_real_ else as.numeric(obs)
    pass <- if (is.na(obs)) FALSE else switch(r$op,
      ">=" = obs >= r$value,
      "in" = obs >= r$value[[1]] && obs <= r$value[[2]])
    data.frame(criterion = r$id, metric = r$metric, observed = obs,
               rule = if (r$op == "in") sprintf("in [%s, %s]", r$value[[1]], r$value[[2]]) else sprintf(">= %s", r$value),
               pass = pass, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
