#!/usr/bin/env Rscript
# Run the R unit tests: Rscript tests/run_tests.R
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
res <- testthat::test_dir(file.path(here, "testthat"), reporter = "summary", stop_on_failure = TRUE)

# A skipped test checks nothing, and testthat skips (rather than fails) a test
# whose package cannot be loaded: if scuttle failed to install, the QC tests
# would be skipped and the run would still pass. No test here is meant to
# skip, so any skip fails the run.
skipped <- as.data.frame(res)
skipped <- skipped[skipped$skipped, c("file", "test"), drop = FALSE]
if (nrow(skipped) > 0) {
  message(paste(sprintf("SKIPPED: %s: %s", skipped$file, skipped$test), collapse = "\n"))
  stop(nrow(skipped), " test(s) skipped; a skip is a failure here (missing package?)", call. = FALSE)
}
invisible(res)
