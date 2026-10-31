# Input/output helpers shared by the bin/ scripts.
#
# Kept dependency-light (Matrix, data.table, jsonlite) so the unit tests run
# without a single-cell framework installed.

#' Parse "--key value" command-line arguments into a named list.
#'
#' Flags without a value are rejected rather than guessed, and every name in
#' `required` must be present. Hyphens in keys become underscores.
parse_args <- function(argv, required = character()) {
  if (length(argv) %% 2 != 0) {
    stop("arguments must come in '--key value' pairs; got: ", paste(argv, collapse = " "))
  }
  keys <- argv[c(TRUE, FALSE)]
  vals <- argv[c(FALSE, TRUE)]
  if (!all(startsWith(keys, "--"))) {
    stop("expected '--key' at positions 1, 3, 5...; got: ", paste(keys, collapse = " "))
  }
  keys <- gsub("-", "_", sub("^--", "", keys))
  if (anyDuplicated(keys)) stop("duplicated argument: ", keys[duplicated(keys)][1])
  out <- as.list(stats::setNames(vals, keys))
  missing <- setdiff(required, names(out))
  if (length(missing)) stop("missing required argument(s): --", paste(gsub("_", "-", missing), collapse = ", --"))
  out
}

#' Read a 10x/STARsolo Market Exchange directory.
#'
#' Accepts plain or gzipped matrix.mtx, barcodes.tsv and features.tsv (or
#' genes.tsv). Returns a list with a genes x barcodes dgCMatrix (rownames are
#' gene IDs), the feature table, and the barcodes with any "-1" suffix kept
#' as written.
read_mex <- function(dir) {
  pick <- function(stems) {
    for (s in stems) for (ext in c("", ".gz")) {
      p <- file.path(dir, paste0(s, ext))
      if (file.exists(p)) return(p)
    }
    stop("none of ", paste(stems, collapse = "/"), " found in ", dir)
  }
  m <- read_mtx(pick("matrix.mtx"))
  barcodes <- data.table::fread(pick("barcodes.tsv"), header = FALSE)[[1]]
  feats <- data.table::fread(pick(c("features.tsv", "genes.tsv")), header = FALSE)
  feats <- as.data.frame(feats)
  names(feats)[1:2] <- c("gene_id", "gene_name")
  if (nrow(m) != nrow(feats)) stop("matrix has ", nrow(m), " rows but features file has ", nrow(feats))
  if (ncol(m) != length(barcodes)) stop("matrix has ", ncol(m), " columns but barcodes file has ", length(barcodes))
  dimnames(m) <- list(feats$gene_id, barcodes)
  list(counts = m, features = feats[, 1:2], barcodes = barcodes)
}

#' Read a coordinate MatrixMarket file as a dgCMatrix.
#'
#' Plain files are read with data.table::fread, which is much faster and
#' leaner than Matrix::readMM on STARsolo's raw matrices (tens of millions
#' of entries). Gzipped files fall back to readMM.
read_mtx <- function(path) {
  if (grepl("\\.gz$", path)) return(methods::as(Matrix::readMM(path), "CsparseMatrix"))
  con <- file(path, "r")
  on.exit(close(con))
  first <- readLines(con, n = 1)
  if (!length(first) || !grepl("^%%MatrixMarket matrix coordinate (integer|real) general", first)) {
    stop("not a general coordinate MatrixMarket file: ", path)
  }
  n_skip <- 1L
  repeat {
    line <- readLines(con, n = 1)
    if (!length(line)) stop("no size line in ", path)
    n_skip <- n_skip + 1L
    if (!startsWith(line, "%")) break
  }
  dims <- suppressWarnings(as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]]))
  if (length(dims) != 3 || anyNA(dims)) stop("bad size line '", line, "' in ", path)
  if (dims[3] == 0) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = dims[1:2]))
  }
  tr <- data.table::fread(path, skip = n_skip, header = FALSE, select = 1:3, col.names = c("i", "j", "x"))
  if (nrow(tr) != dims[3]) stop(path, " declares ", dims[3], " entries but has ", nrow(tr))
  Matrix::sparseMatrix(i = tr$i, j = tr$j, x = as.numeric(tr$x), dims = dims[1:2])
}

#' Split published cell IDs into a sample label and a 16-nt barcode.
#'
#' IDs are split on "_" or ":"; exactly one token must look like a 10x
#' barcode (16 nt, optional "-N" suffix). Everything else, re-joined with
#' "_", is the label. An ID with zero or several barcode-like tokens is an
#' error, because guessing would silently mis-assign cells to samples.
parse_cell_ids <- function(ids) {
  bc_re <- "^[ACGTN]{16}(-[0-9]+)?$"
  parts <- strsplit(ids, "[_:]")
  res <- lapply(seq_along(parts), function(i) {
    tok <- parts[[i]]
    hit <- grepl(bc_re, tok)
    if (sum(hit) != 1) {
      stop(sprintf("cell ID '%s' has %d barcode-like tokens (need exactly 1)", ids[i], sum(hit)))
    }
    c(label = paste(tok[!hit], collapse = "_"), barcode = sub("-[0-9]+$", "", tok[hit]))
  })
  res <- do.call(rbind, res)
  data.frame(id = ids, label = res[, "label"], barcode = res[, "barcode"], stringsAsFactors = FALSE)
}

#' Is this vector of IDs mostly parseable as cell IDs?
looks_like_cell_ids <- function(ids, min_frac = 0.9) {
  ids <- ids[nzchar(ids)]
  if (!length(ids)) return(FALSE)
  ok <- vapply(ids, function(x) tryCatch({ parse_cell_ids(x); TRUE }, error = function(e) FALSE), logical(1))
  mean(ok) >= min_frac
}

#' Stream a published dense counts CSV (plain or .gz) into a sparse matrix.
#'
#' GEO supplementary matrices are often dense text far larger than memory
#' once parsed as doubles, so the file is read once, `chunk_rows` lines at a
#' time, and each chunk is converted to sparse before the next is read.
#'
#' Orientation is detected, not assumed: if the header (after the first
#' column) holds cell IDs the file is genes x cells; if the first column
#' holds cell IDs it is cells x genes. Neither -> error. The result is always
#' genes x cells.
read_published_counts <- function(path, chunk_rows = 2000L) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  header_line <- readLines(con, n = 1)
  if (!length(header_line)) stop("empty file: ", path)
  header <- data.table::fread(text = header_line, header = FALSE, colClasses = "character")
  header <- unname(unlist(header))
  header[is.na(header)] <- ""
  col_ids <- header[-1]

  row_ids <- character()
  blocks <- list()
  repeat {
    lines <- readLines(con, n = chunk_rows)
    if (!length(lines)) break
    dt <- data.table::fread(text = lines, header = FALSE, colClasses = list(character = 1))
    if (ncol(dt) != length(header)) {
      stop("row has ", ncol(dt), " fields but header has ", length(header), " (", path, ")")
    }
    row_ids <- c(row_ids, dt[[1]])
    vals <- as.matrix(dt[, -1, with = FALSE])
    suppressWarnings(storage.mode(vals) <- "double")   # NAs are caught on the next line
    if (anyNA(vals)) stop("non-numeric or missing count values in ", path)
    blocks[[length(blocks) + 1]] <- Matrix::Matrix(vals, sparse = TRUE)
  }
  if (!length(blocks)) stop("no data rows in ", path)
  m <- do.call(rbind, blocks)

  if (looks_like_cell_ids(col_ids)) {
    orientation <- "genes_x_cells"
    dimnames(m) <- list(row_ids, col_ids)
  } else if (looks_like_cell_ids(row_ids)) {
    orientation <- "cells_x_genes"
    dimnames(m) <- list(row_ids, col_ids)
    m <- Matrix::t(m)
  } else {
    stop("could not find cell IDs (label + 16-nt barcode) in either the header or the first column of ", path)
  }
  m <- methods::as(m, "CsparseMatrix")
  attr(m, "orientation") <- orientation
  m
}

#' Write a JSON file with stable formatting.
write_json <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
}
