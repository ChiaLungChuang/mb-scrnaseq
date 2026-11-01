#!/usr/bin/env Rscript
# Run the phase 1 report (lib/R/report.R) on finished synthetic runs and
# check each report against that run's own outputs. CI runs it on the
# nf-test output directories; locally:
#   Rscript tests/check_phase1_report.R .nf-test/tests/*/output
# Directories without a concordance verdict (runs refused before
# concordance) are skipped. nf-test writes the trace to
# <test>/meta/trace.csv (-with-trace) instead of pipeline_info/trace.tsv;
# that file is used when present. Fails unless every report agrees with its run
# and the runs include at least one passing and one failing verdict (the
# swapped-labels test), so both branches of the report are exercised.
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
for (f in list.files(file.path(here, "..", "lib", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

dirs <- commandArgs(trailingOnly = TRUE)
dirs <- dirs[file.exists(file.path(dirs, "concordance", "concordance_summary.verdict.json"))]
if (!length(dirs)) stop("no finished runs (with a concordance verdict) among the directories given")

check <- function(ok, what) if (!isTRUE(ok)) stop(what, call. = FALSE)
seen <- character()
for (d in dirs) {
  nft_trace <- file.path(dirname(d), "meta", "trace.csv")
  trace <- if (!file.exists(file.path(d, "pipeline_info", "trace.tsv")) && file.exists(nft_trace)) nft_trace else NULL
  md <- phase1_report(d, trace = trace)
  v <- jsonlite::read_json(file.path(d, "concordance", "concordance_summary.verdict.json"))
  tab <- utils::read.delim(file.path(d, "concordance", "concordance_summary.tsv"), stringsAsFactors = FALSE)
  tr <- utils::read.delim(if (is.null(trace)) file.path(d, "pipeline_info", "trace.tsv") else trace, stringsAsFactors = FALSE)

  verdict <- grep("^\\*\\*(PASSED|FAILED)\\*\\*", md, value = TRUE)
  check(length(verdict) == 1 && startsWith(verdict, if (isTRUE(v$passed)) "**PASSED**" else "**FAILED**"),
        sprintf("%s: verdict line '%s' does not match the verdict file", d, paste(verdict, collapse = " / ")))
  check(any(grepl("Not pre-registered", md, fixed = TRUE)) == !isTRUE(v$preregistered),
        sprintf("%s: pre-registration statement does not match the verdict file", d))

  head_row <- grep("^\\| sample \\| label \\| C1_", md)
  check(length(head_row) == 1, sprintf("%s: no verdict table", d))
  for (s in unique(tab$sample)) {
    row <- grep(paste0("^\\| ", s, " \\|"), md, value = TRUE)[1]
    cells <- trimws(strsplit(row, "|", fixed = TRUE)[[1]])
    n_fail <- sum(grepl("^\\*\\*FAIL\\*\\*", cells))
    want <- sum(!tab$pass[tab$sample == s]) + if (all(tab$pass[tab$sample == s])) 0 else 1
    check(n_fail == want, sprintf("%s: %s has %d FAIL cells in the report, the table says %d", d, s, n_fail, want))
  }
  n_ran <- sum(tr$status != "CACHED")
  check(any(grepl(sprintf("This run: %d tasks", n_ran), md, fixed = TRUE)),
        sprintf("%s: the report does not count the trace's %d tasks", d, n_ran))
  check(any(grepl(v$code_sha256, md, fixed = TRUE)), sprintf("%s: code SHA-256 missing", d))

  seen <- c(seen, if (isTRUE(v$passed)) "passing" else "failing")
  cat(sprintf("ok  %s: %s, %d samples, %d tasks\n", d, if (isTRUE(v$passed)) "PASSED" else "FAILED",
              length(unique(tab$sample)), n_ran))
}
check(all(c("passing", "failing") %in% seen),
      sprintf("need at least one passing and one failing run; got: %s", paste(seen, collapse = ", ")))
cat(sprintf("phase 1 report checked on %d runs (%d passing, %d failing)\n",
            length(seen), sum(seen == "passing"), sum(seen == "failing")))
