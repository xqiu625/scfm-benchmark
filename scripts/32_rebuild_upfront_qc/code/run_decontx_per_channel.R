#!/usr/bin/env Rscript
# Phase 1 DecontX, replacing scripts/31_donor_level_revision/code/run_decontx_multicelltype.R
# for the rebuild. Two changes from that retrospective-audit version:
#   1. z comes from harmonized_lineage_labels.py (identical marker panel/rule
#      applied to every channel in every cohort), never a cohort's own
#      published annotation.
#   2. batch= is now supplied (one channel/GEM-well id per cell). decontX
#      documents (celda 1.26.0, checked 2026-09-15) that supplying batch
#      fits contamination separately per batch -- "cells run in different
#      channels or assays should be considered different batches" -- so one
#      call per cohort's pooled matrix is still a per-channel model
#      underneath, not a joint fit across channels.
#
# Multiplexed-pool cohorts (COMBAT, Liu/GSE161918) set channel_id to the
# pool/lane id, not the donor id -- see project_302006_rebuild_architecture
# memory for why (cross-donor doublets are real events in those libraries).

suppressPackageStartupMessages({
  library(celda)
  library(SingleCellExperiment)
  library(Matrix)
})

# Two independent 32-bit ceilings in R/Matrix, both hit for real on
# Ren/GSE158055 (nnz = 2,319,108,599 > .Machine$integer.max = 2,147,483,647),
# not hypotheticals:
#   1. Matrix::readMM parses the MatrixMarket header (nrow ncol nnz) via
#      scan() with an integer type -- "scan() expected 'an integer', got
#      '2319108599'".
#   2. Matrix::sparseMatrix()'s DEFAULT output is CsparseMatrix (dgCMatrix,
#      compressed-sparse-column), whose "p" slot holds cumulative per-column
#      non-zero counts -- p's last entry equals nnz, which itself overflows
#      int32 here, independent of how the data was read in. A single
#      data.table::fread() call on the full triplet body ALSO breaks
#      separately (data.table/data.frame row count is capped by the same
#      2^31-1 limit via R's internal row-names machinery: "long vectors not
#      supported yet" / "NAs introduced by coercion to integer range").
#
# Fix for both: chunk fread() calls to stay under the row-count ceiling
# (each chunk's OWN data.table construction succeeds), accumulate into
# preallocated plain atomic vectors (R's long-vector support since R 3.0
# handles vectors longer than 2^31 elements fine -- it is specifically the
# data.frame/data.table row-indexing and the CSC "p" slot that don't), and
# build a TRIPLET-form sparse matrix (repr="T", dgTMatrix -- no cumulative
# pointer, i/j values are bounded only by n_row/n_col which are tiny here)
# instead of letting sparseMatrix() default to CSC.
read_large_mtx <- function(path, chunk_size = 1e9) {
  con <- file(path, "r")
  line <- readLines(con, n = 1)
  stopifnot(startsWith(line, "%%MatrixMarket"))
  n_header_lines <- 1
  repeat {
    line <- readLines(con, n = 1)
    n_header_lines <- n_header_lines + 1
    if (!startsWith(line, "%")) break
  }
  close(con)
  header <- as.numeric(strsplit(trimws(line), "\\s+")[[1]])
  n_row <- as.integer(header[1]); n_col <- as.integer(header[2]); nnz <- header[3]
  cat(sprintf("  read_large_mtx: %d rows x %d cols, %.0f non-zero entries -- chunked data.table::fread\n",
              n_row, n_col, nnz))
  i_vec <- integer(nnz); j_vec <- integer(nnz); x_vec <- double(nnz)
  pos <- 0
  n_chunks <- ceiling(nnz / chunk_size)
  for (c_idx in seq_len(n_chunks)) {
    this_n <- min(chunk_size, nnz - pos)
    dat <- data.table::fread(path, skip = n_header_lines + pos, header = FALSE, nrows = this_n,
                              col.names = c("i", "j", "x"),
                              colClasses = list(integer = 1:2, double = 3))
    stopifnot(nrow(dat) == this_n)
    rng <- (pos + 1):(pos + this_n)
    i_vec[rng] <- dat$i; j_vec[rng] <- dat$j; x_vec[rng] <- dat$x
    pos <- pos + this_n
    rm(dat); gc(verbose = FALSE)
    cat(sprintf("    chunk %d/%d: %.0f / %.0f rows read\n", c_idx, n_chunks, pos, nnz), flush = TRUE)
  }
  stopifnot(pos == nnz)
  sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(n_row, n_col), repr = "T")
}

# Symmetric problem on the way out: decontXcounts() output has the same
# sparsity scale as the input, so writeMM()/CSC casting would hit the same
# "p" slot overflow writing it back. Normalizes to TsparseMatrix (safe,
# doesn't need a cumulative pointer) and writes in chunks via base R
# connections -- no data.table/data.frame construction at all on the write
# path, sidestepping the row-count ceiling entirely rather than needing a
# second chunking scheme for it.
write_large_mtx <- function(m, path, chunk_size = 5e7) {
  mt <- as(m, "TsparseMatrix")
  n_row <- nrow(mt); n_col <- ncol(mt); nnz <- length(mt@x)
  cat(sprintf("  write_large_mtx: %d rows x %d cols, %.0f non-zero entries\n", n_row, n_col, nnz))
  con <- file(path, "w")
  writeLines("%%MatrixMarket matrix coordinate real general", con)
  writeLines("%", con)
  writeLines(sprintf("%d %d %.0f", n_row, n_col, nnz), con)
  n_chunks <- ceiling(nnz / chunk_size)
  for (c_idx in seq_len(n_chunks)) {
    lo <- (c_idx - 1) * chunk_size + 1
    hi <- min(c_idx * chunk_size, nnz)
    rng <- lo:hi
    writeLines(sprintf("%d %d %.6g", mt@i[rng] + 1L, mt@j[rng] + 1L, mt@x[rng]), con)
    cat(sprintf("    write chunk %d/%d\n", c_idx, n_chunks), flush = TRUE)
  }
  close(con)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: run_decontx_per_channel.R <cohort_io_dir> <cohort_tag>\n",
       "  Expects <cohort_io_dir>/counts.mtx, genes.csv, cells.csv\n",
       "  cells.csv must have columns: cell_id, channel_id, harmonized_lineage")
}
IO_DIR <- args[1]
COHORT_TAG <- args[2]

INT32_MAX <- .Machine$integer.max  # 2,147,483,647 -- CSC "p" slot ceiling

cat(sprintf("=== %s ===\n", COHORT_TAG))
cat("Loading filtered count matrix...\n")
counts_t <- read_large_mtx(file.path(IO_DIR, "counts.mtx"))
if (length(counts_t@x) < INT32_MAX) {
  counts <- as(counts_t, "CsparseMatrix")  # small enough: CSC is more efficient for downstream ops
} else {
  cat(sprintf("  nnz=%.0f exceeds CSC's int32 \"p\"-slot ceiling -- keeping TsparseMatrix\n",
              length(counts_t@x)))
  counts <- counts_t
}
rm(counts_t)
genes <- read.csv(file.path(IO_DIR, "genes.csv"))$gene
cells <- read.csv(file.path(IO_DIR, "cells.csv"))
stopifnot(all(c("cell_id", "channel_id", "harmonized_lineage") %in% colnames(cells)))
rownames(counts) <- make.unique(genes)
colnames(counts) <- cells$cell_id
cat(sprintf("  %d genes x %d cells, %d channels\n",
            nrow(counts), ncol(counts), length(unique(cells$channel_id))))
cat("Harmonized lineage composition:\n")
print(table(cells$harmonized_lineage))
cat("Cells per channel (summary):\n")
print(summary(as.numeric(table(cells$channel_id))))

sce <- SingleCellExperiment(assays = list(counts = counts))

cat("Running decontX per-channel (z = harmonized labels, batch = channel_id)...\n")
set.seed(20260914)
sce <- decontX(sce, z = cells$harmonized_lineage, batch = cells$channel_id)

est <- metadata(sce)$decontX$estimates
cat(sprintf("\nContamination fraction summary: mean=%.4f, median=%.4f, max=%.4f\n",
            mean(sce$decontX_contamination), median(sce$decontX_contamination),
            max(sce$decontX_contamination)))
cat("Contamination fraction BY channel (first 20):\n")
print(head(tapply(sce$decontX_contamination, cells$channel_id, mean), 20))
cat("Contamination fraction BY harmonized lineage:\n")
print(tapply(sce$decontX_contamination, cells$harmonized_lineage, mean))

decont_counts <- decontXcounts(sce)
decont_out_path <- file.path(IO_DIR, paste0(COHORT_TAG, "_decontx_counts.mtx"))
decont_nnz <- length(as(decont_counts, "TsparseMatrix")@x)
if (decont_nnz < INT32_MAX) {
  writeMM(as(decont_counts, "CsparseMatrix"), decont_out_path)
} else {
  write_large_mtx(decont_counts, decont_out_path)
}
write.csv(data.frame(cell_id = colnames(sce),
                      channel_id = cells$channel_id,
                      harmonized_lineage = cells$harmonized_lineage,
                      decontX_contamination = sce$decontX_contamination),
          file.path(IO_DIR, paste0(COHORT_TAG, "_decontx_contamination_scores.csv")), row.names = FALSE)

cat(sprintf("SAVED %s_decontx_counts.mtx, %s_decontx_contamination_scores.csv\n", COHORT_TAG, COHORT_TAG))
