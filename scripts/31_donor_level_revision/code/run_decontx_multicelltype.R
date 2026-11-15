#!/usr/bin/env Rscript
# Multi-cell-type DecontX rerun, step 3: run decontX on the combined
# platelet + genuine non-platelet-lineage object, THIS TIME supplying real
# cell-type labels (z) instead of relying on auto-clustering -- since we now
# have actual curated/marker-based lineage identities (Platelet/T/B/NK/Mono/
# DC), not just a platelet-only object with no real heterogeneity to
# cluster on.

suppressPackageStartupMessages({
  library(celda)
  library(SingleCellExperiment)
  library(Matrix)
})

IO_DIR <- "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_multicelltype"

cat("Loading combined multi-cell-type matrix...\n")
counts <- as(readMM(file.path(IO_DIR, "multicelltype_counts.mtx")), "CsparseMatrix")
genes <- read.csv(file.path(IO_DIR, "multicelltype_genes.csv"))$gene
cells <- read.csv(file.path(IO_DIR, "multicelltype_cells.csv"))
rownames(counts) <- make.unique(genes)
colnames(counts) <- cells$cell_id
cat(sprintf("  %d genes x %d cells\n", nrow(counts), ncol(counts)))
cat("Lineage composition:\n")
print(table(cells$lineage))

sce <- SingleCellExperiment(assays = list(counts = counts))

cat("Running decontX WITH real lineage labels as z...\n")
set.seed(20260914)
sce <- decontX(sce, z = cells$lineage)

est <- metadata(sce)$decontX$estimates
run_name <- names(est)[1]
sub <- est[[run_name]]
eta_mat <- sub$eta
if (!is.null(eta_mat)) {
  if (is.null(dim(eta_mat))) {
    eta_df <- data.frame(gene = rownames(counts), eta = as.numeric(eta_mat))
  } else {
    eta_df <- data.frame(gene = rownames(counts), eta = rowMeans(eta_mat))
  }
  eta_df <- eta_df[order(-eta_df$eta), ]
  write.csv(eta_df, file.path(IO_DIR, "multicelltype_ambient_profile_top_genes.csv"), row.names = FALSE)
  cat("Top 30 genes in the estimated ambient/contamination profile (WITH real lineage labels):\n")
  print(head(eta_df, 30))
} else {
  cat("eta not found at expected path -- dumping structure\n")
  str(sub, max.level = 2)
}

cat(sprintf("\nContamination fraction summary: mean=%.4f, median=%.4f, max=%.4f\n",
            mean(sce$decontX_contamination), median(sce$decontX_contamination),
            max(sce$decontX_contamination)))
cat("Contamination fraction BY lineage:\n")
print(tapply(sce$decontX_contamination, cells$lineage, mean))

decont_counts <- decontXcounts(sce)
writeMM(as(decont_counts, "CsparseMatrix"), file.path(IO_DIR, "multicelltype_decontx_counts.mtx"))
write.csv(data.frame(cell_id = colnames(sce),
                      decontX_contamination = sce$decontX_contamination,
                      lineage = cells$lineage),
          file.path(IO_DIR, "multicelltype_decontx_contamination_scores.csv"), row.names = FALSE)

cat("SAVED multicelltype_decontx_counts.mtx, multicelltype_decontx_contamination_scores.csv, "
    , "multicelltype_ambient_profile_top_genes.csv\n")
