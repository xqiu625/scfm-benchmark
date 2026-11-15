#!/bin/bash
#SBATCH --job-name=qiu_decontx
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/qiu_decontx_%j.log
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

module load gcc/16.2.0
module load R/4.5.2
export R_LIBS_USER=/bigdata/godziklab/shared/Xinru/R

Rscript /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/run_decontx_per_channel.R \
  /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/qiu_decontx_input \
  qiu
