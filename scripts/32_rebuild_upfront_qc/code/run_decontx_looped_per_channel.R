#!/usr/bin/env Rscript
# DecontX for cohorts too large for run_decontx_per_channel.R's single-call
# batch= approach -- see split_matrix_per_channel.py and
# project_302006_rebuild_architecture memory for why: celda's own internal
# per-batch subsetting on one unified matrix fails once total non-zeros
# exceed R/Matrix's 32-bit ceiling (2,147,483,647), hit for real on
# Ren/GSE158055 (2,319,108,599 non-zeros). This script gets the SAME
# statistical result -- decontX(z=harmonized_labels, batch=channel_id) is
# already documented to fit each batch separately -- by doing that batch
# split ourselves outside R, one small per-channel matrix at a time, instead
# of asking celda to split one matrix too large for R to hold as a single
# sparse object.
#
# Checkpointed and parallel-safe: each channel writes its own small output
# file immediately (not accumulated in memory for one big write at the end),
# and any channel whose output already exists is skipped on start -- this is
# what makes it safe to (a) re-submit after a partition walltime kill and
# pick up where it left off, and (b) run many array tasks concurrently, each
# claiming a subset of channels via SLURM_ARRAY_TASK_ID with no coordination
# needed, since no two tasks are ever assigned the same channel and a
# half-written file is never left in place (see the atomic-rename note
# below). Per-channel DecontX is small and fast (max ~28M non-zeros/channel
# for Ren, typically much less) -- once the per-channel split exists, this no
# longer needs highmem or a long walltime, so it can run on the `short`
# partition (2h cap) across many parallel array tasks instead of one long
# sequential highmem job.
#
# Usage: run_decontx_looped_per_channel.R <split_dir> <genes_csv> <out_dir> <cohort_tag> [task_id] [n_tasks]
#   <split_dir>: output of split_matrix_per_channel.py (one <channel_id>/
#                subdirectory per channel, each with counts.mtx + cells.csv)
#   <genes_csv>: the cohort's genes.csv (shared across all channels)
#   <out_dir>:   where to write per-channel <cohort_tag>_<channel_id>_decontx_counts.mtx
#                and _decontx_contamination_scores.csv (deliberately NOT one
#                combined counts file -- that would reintroduce the same
#                ceiling on the way out; combine_decontx_scores.R combines the
#                per-channel contamination-score CSVs afterward, which is safe
#                since that table has one row per CELL, not per non-zero)
#   [task_id], [n_tasks]: optional 0-indexed task id and total task count for
#                a parallel array job -- this task processes channels where
#                (index - 1) %% n_tasks == task_id. Omit both to process every
#                channel in one sequential run (the original behavior).

suppressPackageStartupMessages({
  library(celda)
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: run_decontx_looped_per_channel.R <split_dir> <genes_csv> <out_dir> <cohort_tag> [task_id] [n_tasks]")
}
SPLIT_DIR <- args[1]
GENES_CSV <- args[2]
OUT_DIR <- args[3]
COHORT_TAG <- args[4]
TASK_ID <- if (length(args) >= 6) as.integer(args[5]) else NULL
N_TASKS <- if (length(args) >= 6) as.integer(args[6]) else NULL
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

genes <- make.unique(read.csv(GENES_CSV)$gene)
all_channels <- sort(list.dirs(SPLIT_DIR, full.names = FALSE, recursive = FALSE))
all_channels <- all_channels[all_channels != ""]

if (!is.null(TASK_ID)) {
  keep <- ((seq_along(all_channels) - 1) %% N_TASKS) == TASK_ID
  channels <- all_channels[keep]
  cat(sprintf("=== %s: task %d/%d, %d of %d total channels assigned ===\n",
              COHORT_TAG, TASK_ID, N_TASKS, length(channels), length(all_channels)))
} else {
  channels <- all_channels
  cat(sprintf("=== %s: %d channels (single sequential run) ===\n", COHORT_TAG, length(channels)))
}

n_done <- 0
n_skipped <- 0
for (i in seq_along(channels)) {
  ch <- channels[i]
  out_mtx <- file.path(OUT_DIR, paste0(COHORT_TAG, "_", ch, "_decontx_counts.mtx"))
  out_scores <- file.path(OUT_DIR, paste0(COHORT_TAG, "_", ch, "_decontx_scores.csv"))
  if (file.exists(out_mtx) && file.exists(out_scores)) {
    n_skipped <- n_skipped + 1
    next
  }

  ch_dir <- file.path(SPLIT_DIR, ch)
  counts <- as(readMM(file.path(ch_dir, "counts.mtx")), "CsparseMatrix")
  cells <- read.csv(file.path(ch_dir, "cells.csv"))
  stopifnot(nrow(cells) == ncol(counts))
  rownames(counts) <- genes
  colnames(counts) <- cells$cell_id
  stopifnot(length(unique(cells$channel_id)) == 1, cells$channel_id[1] == ch)

  sce <- SingleCellExperiment(assays = list(counts = counts))
  set.seed(20260914)  # matches run_decontx_per_channel.R's seed for consistency
  sce <- decontX(sce, z = cells$harmonized_lineage, verbose = FALSE)

  # Write to a temp path then rename -- atomic on the same filesystem, so a
  # killed/timed-out job never leaves a half-written file that a later
  # checkpointed re-run would mistake for "done" via the file.exists() check
  # above.
  tmp_mtx <- paste0(out_mtx, ".tmp")
  writeMM(as(decontXcounts(sce), "CsparseMatrix"), tmp_mtx)
  file.rename(tmp_mtx, out_mtx)

  scores <- data.frame(cell_id = colnames(sce), channel_id = ch,
                        harmonized_lineage = cells$harmonized_lineage,
                        decontX_contamination = sce$decontX_contamination)
  tmp_scores <- paste0(out_scores, ".tmp")
  write.csv(scores, tmp_scores, row.names = FALSE)
  file.rename(tmp_scores, out_scores)

  n_done <- n_done + 1
  cat(sprintf("%s  [%d/%d assigned] %s: %d cells, mean contam=%.4f (skipped %d already-done)\n",
              format(Sys.time()), i, length(channels), ch, ncol(counts),
              mean(sce$decontX_contamination), n_skipped), flush = TRUE)
}

cat(sprintf("\nDONE this task: %d channels processed, %d already-done channels skipped.\n", n_done, n_skipped))
