#!/bin/bash
#SBATCH --job-name=convert_wilk_rds
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/convert_wilk_rds_%j.log
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

module load gcc/16.2.0
module load R/4.5.2
export R_LIBS_USER=/bigdata/godziklab/shared/Xinru/R

Rscript /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/convert_wilk_rds_to_mtx.R \
  /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE174072_wilk \
  /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE174072_wilk/converted
