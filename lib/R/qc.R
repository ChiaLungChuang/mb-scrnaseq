# Per-sample cell QC.
#
# Cells are flagged, not removed: thresholds are adaptive (median +/- n MADs,
# computed per sample) and recorded alongside the flags, so the decision to
# drop cells is made later, in the open, rather than inside this step.

#' Per-barcode QC metrics.
qc_metrics <- function(counts, gene_names, mito_pattern = "^MT-") {
  if (length(gene_names) != nrow(counts)) stop("gene_names length != nrow(counts)")
  mito <- grepl(mito_pattern, gene_names, ignore.case = TRUE)
  total <- Matrix::colSums(counts)
  data.frame(
    barcode = colnames(counts),
    total_umi = as.numeric(total),
    detected_genes = as.numeric(Matrix::colSums(counts > 0)),
    pct_mito = ifelse(total > 0, 100 * as.numeric(Matrix::colSums(counts[mito, , drop = FALSE])) / total, NA_real_),
    stringsAsFactors = FALSE
  )
}

#' Flag outlier cells with scuttle::isOutlier and record the thresholds.
#'
#' Low total UMI and low detected genes are judged on the log scale (lower
#' tail); mitochondrial percentage on the linear scale (upper tail). When a
#' metric's MAD is zero the threshold collapses onto the median, which would
#' flag every cell that differs at all, so that case is reported rather than
#' silently applied.
qc_flag <- function(df, nmads = 3) {
  one <- function(x, type, log) {
    o <- scuttle::isOutlier(x, nmads = nmads, type = type, log = log)
    th <- attr(o, "thresholds")
    mad_val <- stats::mad(if (log) log(x[x > 0]) else x, na.rm = TRUE)
    list(flag = as.logical(o), lower = unname(th["lower"]), higher = unname(th["higher"]), mad_zero = isTRUE(mad_val == 0))
  }
  tot <- one(df$total_umi, "lower", TRUE)
  det <- one(df$detected_genes, "lower", TRUE)
  mt <- one(df$pct_mito, "higher", FALSE)
  df$flag_low_umi <- tot$flag
  df$flag_low_genes <- det$flag
  df$flag_high_mito <- mt$flag
  df$flag_any <- df$flag_low_umi | df$flag_low_genes | df$flag_high_mito
  thresholds <- list(
    nmads = nmads,
    total_umi = list(lower = tot$lower, mad_zero = tot$mad_zero),
    detected_genes = list(lower = det$lower, mad_zero = det$mad_zero),
    pct_mito = list(higher = mt$higher, mad_zero = mt$mad_zero),
    n_cells = nrow(df),
    n_flagged = sum(df$flag_any)
  )
  list(cells = df, thresholds = thresholds)
}
