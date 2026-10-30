#!/usr/bin/env Rscript
# Combine per-sample concordance results into one table and one verdict.
#
# Usage: summarise_concordance.R --lib DIR --out-prefix PREFIX JSON [JSON ...]
# Writes PREFIX.tsv (one row per sample x criterion) and PREFIX.verdict.json.
# The overall verdict passes only if every sample passes every criterion.
argv <- commandArgs(trailingOnly = TRUE)
lib <- argv[which(argv == "--lib") + 1]
for (f in list.files(lib, pattern = "\\.R$", full.names = TRUE)) source(f)
opt_end <- 0
while (opt_end + 2 <= length(argv) && startsWith(argv[opt_end + 1], "--")) opt_end <- opt_end + 2
a <- parse_args(argv[seq_len(opt_end)], required = c("lib", "out_prefix"))
files <- argv[-seq_len(opt_end)]
if (!length(files)) stop("no concordance JSON files given")

rows <- lapply(files, function(p) {
  j <- jsonlite::read_json(p, simplifyVector = TRUE)
  cbind(sample = j$sample, label = j$label, j$criteria, stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$sample, tab$criterion), ]
utils::write.table(tab, paste0(a$out_prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

per_sample <- tapply(tab$pass, tab$sample, all)
write_json(list(
  passed = all(per_sample),
  n_samples = length(per_sample),
  n_samples_passed = sum(per_sample),
  failing_samples = names(per_sample)[!per_sample],
  failing_criteria = unique(tab$criterion[!tab$pass])
), paste0(a$out_prefix, ".verdict.json"))
message(sprintf("concordance: %d/%d samples passed all criteria", sum(per_sample), length(per_sample)))
