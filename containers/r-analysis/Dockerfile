# R environment for the analysis steps (QC, published-matrix import, concordance).
#
# rocker/r-ver pins CRAN to a dated Posit Package Manager snapshot for its R
# version, and Bioconductor is pinned to the release that matches R 4.4, so
# rebuilding this image later installs the same package versions.
#
# Build:  docker build -t mb-scrnaseq-r:dev containers/r-analysis
# CI builds it on every push and publishes it to GHCR on version tags.
FROM rocker/r-ver:4.4.3

LABEL org.opencontainers.image.source="https://github.com/ChiaLungChuang/mb-scrnaseq" \
      org.opencontainers.image.description="R environment for mb-scrnaseq" \
      org.opencontainers.image.licenses="MIT"

# procps provides `ps`, which Nextflow uses to record task resource usage.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled data.table jsonlite yaml BiocManager testthat withr \
 && Rscript -e 'BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)' \
 && Rscript -e 'BiocManager::install("scuttle", ask = FALSE, update = FALSE)' \
 && Rscript -e 'for (p in c("Matrix", "data.table", "jsonlite", "yaml", "scuttle")) library(p, character.only = TRUE)' \
 && rm -rf /tmp/downloaded_packages
