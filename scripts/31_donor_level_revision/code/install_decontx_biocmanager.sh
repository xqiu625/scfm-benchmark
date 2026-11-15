#!/bin/bash
#SBATCH --job-name=install_decontx_bioc
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/install_decontx_bioc_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/install_decontx_bioc_%j.err

module load R/4.5.2
module load gcc/12.2.0

R_LIB=/bigdata/godziklab/shared/Xinru/302006/r_libs_decontx
mkdir -p "$R_LIB"
rm -rf "$R_LIB/00LOCK-bluster" "$R_LIB/00LOCK-scran" "$R_LIB/00LOCK-celda"

Rscript -e '
r_lib <- "/bigdata/godziklab/shared/Xinru/302006/r_libs_decontx"
.libPaths(c(r_lib, .libPaths()))
library(BiocManager)
BiocManager::install(c("bluster", "scran", "celda", "SingleCellExperiment", "Matrix"),
                      lib = r_lib, update = FALSE, ask = FALSE)
ok <- requireNamespace("celda", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE)
if (ok) cat("DECONTX_ENV_INSTALL_SUCCESS\n") else cat("DECONTX_ENV_INSTALL_FAILED\n")
'
