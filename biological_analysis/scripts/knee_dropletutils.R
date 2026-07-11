#!/usr/bin/env Rscript
# Gold-standard barcode-rank knee/inflection via DropletUtils::barcodeRanks (Lun et al. 2019).
# barcodeRanks fits a smoothing spline to the log-log rank/total-count curve, then reports the knee
# (point of maximum curvature) and the inflection (point of minimum gradient, the cell/empty cliff).
# We print both as a UMI threshold and as the corresponding rank (number of barcodes at or above it).
#
#   Rscript knee_dropletutils.R name1=/path/to/Gene/raw name2=/path/to/Gene/raw ...
#
# Each path is a STARsolo Gene/raw dir containing matrix.mtx(.gz), barcodes.tsv(.gz), features.tsv(.gz).
# Requires R with DropletUtils + Matrix (Bioconductor). If absent:
#   Rscript -e 'install.packages("BiocManager"); BiocManager::install("DropletUtils")'

suppressMessages({library(DropletUtils); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: knee_dropletutils.R name=Gene/raw [name=Gene/raw ...]")

for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  name <- kv[1]; d <- kv[2]
  mtx <- list.files(d, pattern = "^matrix\\.mtx", full.names = TRUE)[1]
  con <- if (grepl("\\.gz$", mtx)) gzfile(mtx) else mtx
  m   <- readMM(con)                                  # genes x barcodes
  br  <- barcodeRanks(m)                              # default lower=100 (ambient cutoff)
  kc  <- metadata(br)$knee
  ic  <- metadata(br)$inflection
  tot <- br$total
  cat(sprintf("%-10s  knee: UMI>=%.0f at rank %d   |   inflection: UMI>=%.0f at rank %d   |   n_barcodes %d\n",
              name, kc, sum(tot >= kc), ic, sum(tot >= ic), length(tot)))
}
