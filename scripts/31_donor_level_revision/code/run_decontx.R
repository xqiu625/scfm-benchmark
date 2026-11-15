#!/usr/bin/env Rscript
# Run DecontX (celda) on the canonical within-COVID cohort exported by
# export_for_decontx.py. DecontX estimates per-cell ambient-RNA
# contamination directly from a FILTERED (cell-called) count matrix --
# unlike SoupX/CellBender, it does not require the raw unfiltered/
# empty-droplet matrix, which is not available for these already-cell-
# called public source datasets (the feasibility gap this whole check was
# built around).
#
# No cluster labels (z) are supplied -- DecontX's documented default
# behavior is to auto-cluster the cells internally (via celda) when z is
# omitted; this is the standard, most common usage pattern and avoids an
# extra, potentially outcome-correlated clustering step of our own.

suppressPackageStartupMessages({
  library(celda)
  library(SingleCellExperiment)
  library(Matrix)
})

IO_DIR <- "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_io"

cat("Loading exported matrix...\n")
counts <- as(readMM(file.path(IO_DIR, "counts.mtx")), "CsparseMatrix")  # genes x cells
genes <- read.csv(file.path(IO_DIR, "genes.csv"))$gene
cells <- read.csv(file.path(IO_DIR, "cells.csv"))
rownames(counts) <- make.unique(genes)
colnames(counts) <- cells$cell_id
cat(sprintf("  %d genes x %d cells\n", nrow(counts), ncol(counts)))

sce <- SingleCellExperiment(assays = list(counts = counts))

cat("Running decontX (no z supplied -- internal auto-clustering)...\n")
set.seed(20260913)
sce <- decontX(sce)

cat(sprintf("Contamination fraction summary: mean=%.4f, median=%.4f, max=%.4f\n",
            mean(sce$decontX_contamination), median(sce$decontX_contamination),
            max(sce$decontX_contamination)))

decont_counts <- decontXcounts(sce)
writeMM(as(decont_counts, "CsparseMatrix"), file.path(IO_DIR, "decontx_counts.mtx"))
write.csv(data.frame(cell_id = colnames(sce),
                      decontX_contamination = sce$decontX_contamination,
                      decontX_clusters = sce$decontX_clusters),
          file.path(IO_DIR, "decontx_contamination_scores.csv"), row.names = FALSE)

cat("SAVED decontx_counts.mtx and decontx_contamination_scores.csv\n")
