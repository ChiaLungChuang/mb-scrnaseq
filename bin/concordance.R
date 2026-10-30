#!/usr/bin/env Rscript
# Compare one sample's STARsolo counts with the published counts.
#
# Usage: concordance.R --lib DIR --sample ID --label LABEL --raw-dir DIR
#                      --filtered-dir DIR --published RDS --criteria YAML
#                      --out-prefix PREFIX
# Writes PREFIX.concordance.json (metrics + per-criterion verdicts) and
# PREFIX.per_cell.tsv.gz. A failed criterion is a result, not a crash: the
# step exits 0 and the verdict is recorded.
argv <- commandArgs(trailingOnly = TRUE)
lib <- argv[which(argv == "--lib") + 1]
for (f in list.files(lib, pattern = "\\.R$", full.names = TRUE)) source(f)
a <- parse_args(argv, required = c("lib", "sample", "label", "raw_dir", "filtered_dir",
                                   "published", "criteria", "out_prefix"))

crit_all <- yaml::read_yaml(a$criteria)
crit <- crit_all$phase1_concordance
if (is.null(crit)) stop("criteria file has no phase1_concordance block: ", a$criteria)

ours <- read_mex(a$raw_dir)
called <- read_mex(a$filtered_dir)$barcodes
published <- readRDS(a$published)

res <- concordance_metrics(ours, published, a$label, called_barcodes = called)
ev <- evaluate_criteria(res$metrics, crit)

write_json(list(
  sample = a$sample,
  label = a$label,
  criteria_status = crit_all$status,
  passed = all(ev$pass),
  metrics = res$metrics,
  criteria = ev
), paste0(a$out_prefix, ".concordance.json"))
con <- gzfile(paste0(a$out_prefix, ".per_cell.tsv.gz"), "w")
utils::write.table(res$per_cell, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

message(sprintf("%s: %s (%d/%d criteria passed)", a$sample,
                if (all(ev$pass)) "PASS" else "FAIL", sum(ev$pass), nrow(ev)))
