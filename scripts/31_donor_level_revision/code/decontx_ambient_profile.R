#!/usr/bin/env Rscript
# External review point B/G-diagnostic: report the top genes of DecontX's
# estimated ambient/contamination profile (metadata(sce)$decontX$estimates
# $...$estRmat / eta), which is independently informative and was not saved
# by the original run_decontx.R (only decontX_contamination and
# decontX_clusters were extracted). Reruns decontX with the same inputs and
# seed -- ~2 min -- so this captures the full estimates object this time.
#
# Also directly addresses review point B's factual premise: this cohort
# (10,737 cells) is platelets-only (no heterogeneous background population
# in the input matrix), so DecontX's ambient-profile estimate is identified
# only from within-platelet compositional variation, not from a genuine
# platelet-vs-other-celltype contrast. Report that limitation alongside the
# profile itself rather than only alongside the correction's correctness.

suppressPackageStartupMessages({
  library(celda)
  library(SingleCellExperiment)
  library(Matrix)
})

IO_DIR <- "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_io"

cat("Loading exported matrix...\n")
counts <- as(readMM(file.path(IO_DIR, "counts.mtx")), "CsparseMatrix")
genes <- read.csv(file.path(IO_DIR, "genes.csv"))$gene
cells <- read.csv(file.path(IO_DIR, "cells.csv"))
rownames(counts) <- make.unique(genes)
colnames(counts) <- cells$cell_id
cat(sprintf("  %d genes x %d cells (platelets only -- no heterogeneous background cell type in this matrix)\n",
            nrow(counts), ncol(counts)))

sce <- SingleCellExperiment(assays = list(counts = counts))

cat("Running decontX (same seed as the original run)...\n")
set.seed(20260913)
sce <- decontX(sce)

est <- metadata(sce)$decontX$estimates
cat("Names under metadata(sce)$decontX$estimates:\n")
print(names(est))
run_name <- names(est)[1]
cat(sprintf("Using estimate set: %s\n", run_name))
sub <- est[[run_name]]
cat("Fields available:\n")
print(names(sub))

# eta = estimated contamination (ambient) distribution over genes, per cluster
# phi = estimated native distribution over genes, per cluster
eta <- sub$estNativeCounts  # placeholder, corrected below if absent
if (!is.null(sub$eta)) {
  eta_mat <- sub$eta
} else if (!is.null(sub$estimates$eta)) {
  eta_mat <- sub$estimates$eta
} else {
  eta_mat <- NULL
}

if (is.null(eta_mat)) {
  cat("eta not found at expected path -- dumping full structure for manual inspection\n")
  str(sub, max.level = 2)
} else {
  # eta_mat: genes x clusters (or a vector if 1 cluster/background)
  if (is.null(dim(eta_mat))) {
    eta_df <- data.frame(gene = rownames(counts), eta = as.numeric(eta_mat))
  } else {
    eta_mean <- rowMeans(eta_mat)
    eta_df <- data.frame(gene = rownames(counts), eta = eta_mean)
  }
  eta_df <- eta_df[order(-eta_df$eta), ]
  write.csv(eta_df, file.path(IO_DIR, "decontx_ambient_profile_top_genes.csv"), row.names = FALSE)
  cat("Top 30 genes in the estimated ambient/contamination profile:\n")
  print(head(eta_df, 30))
}

cat(sprintf("\nContamination fraction summary: mean=%.4f, median=%.4f, max=%.4f\n",
            mean(sce$decontX_contamination), median(sce$decontX_contamination),
            max(sce$decontX_contamination)))
cat("DONE\n")
