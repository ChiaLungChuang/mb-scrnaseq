#!/usr/bin/env Rscript
# Full-size synthetic inputs for the R steps' resource rehearsal
# (tests/scale/README.md). Writes, under OUT:
#
#   raw/       STARsolo-style raw matrix: 33,538 genes x 6,794,880 barcodes
#              (the v3 whitelist size); CELLS called cells with realistic
#              depth, AMBIENT barcodes with a handful of UMIs each
#   filtered/  the CELLS called cells only
#   published_counts.csv.gz
#              dense genes x cells CSV, 26,841 genes x 39,946 cells, with
#              the 30 labels and per-label cell counts of the real matrix
#              (docs/proposals/gse155446_runs.tsv: names and counts of cell
#              IDs only). Label LABEL's cells are a subset of the called
#              cells, with the same counts, so concordance has real work.
#
# Nothing here comes from count values of the real data. The output is
# several GB and is never committed.
#
# Usage: make_scale_data.R --out DIR [--cells 5000] [--ambient 300000]
#                          [--umis 5000] [--label 1167] [--seed 1]
#                          [--barcodes 6794880]
#   --umis  median UMIs per cell: 5000 gives ~2,000 genes per cell (7% of
#           the published matrix non-zero), 12000 ~3,300 (a dense case)
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
argv <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d) { i <- match(paste0("--", k), argv); if (is.na(i)) d else argv[i + 1] }
out <- opt("out", NA); if (is.na(out)) stop("--out DIR is required")
n_cells <- as.integer(opt("cells", 5000)); n_ambient <- as.integer(opt("ambient", 300000))
median_umis <- as.numeric(opt("umis", 5000))
label <- opt("label", "1167"); set.seed(as.integer(opt("seed", 1)))
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
runs <- fread(file.path(here, "..", "..", "docs", "proposals", "gse155446_runs.tsv"),
              colClasses = list(character = "published_label"))
stopifnot(label %in% runs$published_label)

n_genes <- 33538L; n_pub_genes <- 26841L
n_barcodes <- as.integer(opt("barcodes", 6794880))   # smaller only to try the script out
t0 <- Sys.time(); say <- function(...) message(sprintf("[%5.0f s] ", as.numeric(Sys.time() - t0, units = "secs")), ...)

# Genes: Ensembl-like IDs, symbols with 13 mitochondrial genes and 20
# duplicated symbols (as in GRCh38-3.0.0); heavy-tailed expression levels,
# a fifth of genes never expressed.
gene_id <- sprintf("ENSG%011d", seq_len(n_genes))
symbol <- sprintf("GENE%05d", seq_len(n_genes))
symbol[1:13] <- sprintf("MT-%02d", 1:13)
symbol[1001:1020] <- symbol[2001:2020]
p <- rlnorm(n_genes, 0, 2.4); p[sample.int(n_genes, n_genes %/% 5)] <- 0; p[1:13] <- p[1:13] * 40
p <- p / sum(p)

# n UMIs per cell (lognormal, median --umis) spread over genes by p.
umis <- function(totals, prob) {
  dt <- data.table(j = rep.int(seq_along(totals), totals),
                   i = sample.int(length(prob), sum(totals), replace = TRUE, prob = prob))
  dt[, .(x = .N), keyby = .(j, i)]
}
cell_totals <- function(n) pmax(500L, as.integer(round(rlnorm(n, log(median_umis), 0.5))))

# Barcodes: every 16-mer index -> ACGT, as two 8-mers from a lookup table.
kmers8 <- do.call(paste0, rev(expand.grid(rep(list(c("A", "C", "G", "T")), 8), stringsAsFactors = FALSE)))
barcode <- function(idx0) paste0(kmers8[idx0 %/% 65536L + 1L], kmers8[idx0 %% 65536L + 1L])
wl <- barcode(seq_len(n_barcodes) - 1L)
pick <- sample.int(n_barcodes, n_cells + n_ambient)
cell_bc <- pick[seq_len(n_cells)]; amb_bc <- pick[-seq_len(n_cells)]
say("barcodes ready")

cells <- umis(cell_totals(n_cells), p)
say(sprintf("called cells: %d entries, median %d genes per cell", nrow(cells),
            as.integer(median(cells[, .N, by = j]$N))))
amb <- umis(1L + rpois(n_ambient, 2), p)
raw <- rbind(cells[, .(i, j = cell_bc[j], x)], amb[, .(i, j = amb_bc[j], x)])
setkey(raw, j, i)

write_mex <- function(dir, trip, barcodes) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  mtx <- file.path(dir, "matrix.mtx")
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "%",
               sprintf("%d %d %d", n_genes, length(barcodes), nrow(trip))), mtx)
  fwrite(trip[, .(i, j, x)], mtx, append = TRUE, sep = " ", col.names = FALSE)
  fwrite(data.table(barcodes), file.path(dir, "barcodes.tsv"), col.names = FALSE)
  fwrite(data.table(gene_id, symbol, "Gene Expression"), file.path(dir, "features.tsv"),
         sep = "\t", col.names = FALSE)
}
write_mex(file.path(out, "raw"), raw, wl)
write_mex(file.path(out, "filtered"), cells, wl[cell_bc])
say(sprintf("raw matrix: %d entries; filtered: %d cells", nrow(raw), n_cells))
rm(raw, amb); invisible(gc())

# Published matrix: 26,841 genes with unique symbols (mitochondrial ones
# included), 39,946 cells. Label `label` = a subset of our called cells
# with our counts; the other labels are drawn from the same model.
keep_g <- sort(c(1:13, sample(setdiff(which(!duplicated(symbol) & !(symbol %in% symbol[duplicated(symbol)])), 1:13),
                              n_pub_genes - 13L)))
row_of <- match(seq_len(n_genes), keep_g)           # our gene -> published row (NA if not kept)
n_lab <- runs[published_label == label, published_cells]
mine <- sort(sample.int(n_cells, n_lab))
pub <- list(cells[j %in% mine, .(i = row_of[i], j = match(j, mine), x)][!is.na(i)])
pub_ids <- paste0(label, "_", wl[cell_bc[mine]])
offset <- n_lab
pp <- p[keep_g] / sum(p[keep_g])
for (k in which(runs$published_label != label)) {
  n <- runs$published_cells[k]
  d <- umis(cell_totals(n), pp)
  pub[[length(pub) + 1]] <- d[, .(i, j = j + offset, x)]
  pub_ids <- c(pub_ids, paste0(runs$published_label[k], "_", barcode(sample.int(4^16 - 1, n))))
  offset <- offset + n
}
pub <- rbindlist(pub); setkey(pub, i, j)
stopifnot(offset == 39946L, length(pub_ids) == offset)
say(sprintf("published: %d genes x %d cells, %d non-zero (%.1f%%)", n_pub_genes, offset, nrow(pub),
            100 * nrow(pub) / (as.numeric(n_pub_genes) * offset)))

csv <- file.path(out, "published_counts.csv")
writeLines(paste(c("", pub_ids), collapse = ","), csv)
for (a in seq(1L, n_pub_genes, by = 1000L)) {
  b <- min(a + 999L, n_pub_genes)
  dense <- matrix(0L, nrow = b - a + 1L, ncol = offset)
  s <- pub[i >= a & i <= b]
  dense[cbind(s$i - a + 1L, s$j)] <- s$x
  fwrite(data.table(gene = symbol[keep_g[a:b]], dense), csv, append = TRUE, col.names = FALSE)
}
say("dense CSV written; compressing")
if (system2("gzip", c("-f", shQuote(csv))) != 0) stop("gzip failed")
say(sprintf("done: %s", out))
