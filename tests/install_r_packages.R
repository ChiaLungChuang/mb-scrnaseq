#!/usr/bin/env Rscript
# Install the R packages the unit tests need: Rscript tests/install_r_packages.R
#
# CI's r-unit-tests job runs this on a cache miss. The cache key includes
# this file's hash, so editing the package list here refreshes the cache.
# BiocManager::install() only warns when a package fails to build, so every
# package is loaded at the end: a failed install stops here with an error.
cran <- c("data.table", "jsonlite", "yaml", "testthat", "withr", "BiocManager")
bioc <- c("scuttle")
install.packages(cran)
BiocManager::install(bioc, ask = FALSE, update = FALSE)
for (p in c("Matrix", cran, bioc)) library(p, character.only = TRUE)
cat("R packages installed and loadable:", paste(c("Matrix", cran, bioc), collapse = ", "), "\n")
