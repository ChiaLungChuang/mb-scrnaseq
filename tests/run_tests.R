#!/usr/bin/env Rscript
# Run the R unit tests: Rscript tests/run_tests.R
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
res <- testthat::test_dir(file.path(here, "testthat"), reporter = "summary", stop_on_failure = TRUE)
invisible(res)
