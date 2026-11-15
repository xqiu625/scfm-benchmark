#!/usr/bin/env Rscript
# Converts Wilk/GSE174072's double-gzip-compressed Seq-Well RDS matrices
# (list of exon/intron/spanning dgCMatrix, genes x cells, gene SYMBOL rownames
# only -- no Ensembl ID available, unlike every other cohort's raw data) into
# plain MatrixMarket + barcodes.tsv + features.tsv triples, one per GSM, so
# the existing Python build-script pattern (symbol-based Ensembl-ID
# intersection logic doesn't apply here -- falls back to symbol intersection
# directly) can consume them uniformly.
#
# Uses EXON counts only, not exon+intron+spanning summed: this is standard
# scRNA-seq (PBMC, not single-nucleus), where intronic/spanning reads are
# pre-mRNA/ambient background rather than mature transcript signal -- exon
# also has the broadest gene panel and highest total UMI count of the three
# layers in a spot check (33,109 genes / 13.69M UMI vs intron's 23,423 genes
# / 9.39M UMI and spanning's 20,301 genes / 0.37M UMI on GSM5285683).
#
# The double-gzip is real, not a mistake: `file` on the .rds.gz shows valid
# gzip, but readRDS() on it (even via gzfile()) fails "unknown input format"
# -- the DECOMPRESSED content is ITSELF gzip-compressed (R's own internal
# saveRDS(compress=TRUE) gzip, then the whole file gzip'd again for GEO
# distribution). Fixed by piping through zcat twice before readRDS.
#
# Usage: Rscript convert_wilk_rds_to_mtx.R <raw_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(Matrix))

rds_files <- list.files(raw_dir, pattern = "_cell\\.counts\\.matrices\\.rds\\.gz$", full.names = TRUE)
cat(sprintf("%d RDS files found\n", length(rds_files)))

for (rds_path in rds_files) {
  base <- basename(rds_path)
  base <- sub("_cell\\.counts\\.matrices\\.rds\\.gz$", "", base)
  ch_dir <- file.path(out_dir, base)
  if (file.exists(file.path(ch_dir, "matrix.mtx"))) {
    cat(sprintf("%s: already converted, skipping\n", base))
    next
  }
  dir.create(ch_dir, showWarnings = FALSE, recursive = TRUE)

  tmp <- tempfile(fileext = ".rds")
  system(sprintf("zcat %s | zcat > %s", shQuote(rds_path), shQuote(tmp)))
  x <- tryCatch(readRDS(tmp), error = function(e) NULL)
  file.remove(tmp)
  if (is.null(x) || is.null(x$exon)) {
    cat(sprintf("%s: FAILED to load or no exon matrix, skipping\n", base))
    next
  }
  m <- x$exon  # genes x cells, dgCMatrix
  writeMM(m, file.path(ch_dir, "matrix.mtx"))
  writeLines(colnames(m), file.path(ch_dir, "barcodes.tsv"))
  writeLines(rownames(m), file.path(ch_dir, "features.tsv"))
  cat(sprintf("%s: %d genes x %d cells written\n", base, nrow(m), ncol(m)))
}
cat("Done.\n")
