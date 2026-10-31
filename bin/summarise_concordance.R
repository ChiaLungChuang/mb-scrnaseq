#!/usr/bin/env Rscript
# Combine per-sample concordance results into one table and one verdict.
#
# Usage: summarise_concordance.R --lib DIR --out-prefix PREFIX
#            --preregistered true|false --criteria-sha256 HEX --code-sha256 HEX
#            JSON [JSON ...]
# Writes PREFIX.tsv (one row per sample x criterion) and PREFIX.verdict.json.
# The overall verdict passes only if every sample passes every criterion.
# The verdict also says whether the criteria were pre-registered, so a test
# run's verdict can never be mistaken for a reportable one.
argv <- commandArgs(trailingOnly = TRUE)
lib <- argv[which(argv == "--lib") + 1]
for (f in list.files(lib, pattern = "\\.R$", full.names = TRUE)) source(f)
opt_end <- 0
while (opt_end + 2 <= length(argv) && startsWith(argv[opt_end + 1], "--")) opt_end <- opt_end + 2
a <- parse_args(argv[seq_len(opt_end)],
                required = c("lib", "out_prefix", "preregistered", "criteria_sha256", "code_sha256"))
if (!a$preregistered %in% c("true", "false")) stop("--preregistered must be true or false")
files <- argv[-seq_len(opt_end)]
if (!length(files)) stop("no concordance JSON files given")

rows <- lapply(files, function(p) {
  j <- jsonlite::read_json(p, simplifyVector = TRUE)
  crit <- j$criteria
  if (is.null(crit$observed)) crit$observed <- NA_real_   # all-NA metrics are omitted from JSON rows
  crit <- crit[, c("criterion", "metric", "observed", "rule", "pass")]
  cbind(sample = j$sample, label = j$label, crit, stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$sample, tab$criterion), ]
utils::write.table(tab, paste0(a$out_prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

per_sample <- tapply(tab$pass, tab$sample, all)
write_json(list(
  preregistered = a$preregistered == "true",
  criteria_sha256 = a$criteria_sha256,
  code_sha256 = a$code_sha256,
  passed = all(per_sample),
  n_samples = length(per_sample),
  n_samples_passed = sum(per_sample),
  failing_samples = names(per_sample)[!per_sample],
  failing_criteria = unique(tab$criterion[!tab$pass])
), paste0(a$out_prefix, ".verdict.json"))
message(sprintf("concordance: %d/%d samples passed all criteria", sum(per_sample), length(per_sample)))
