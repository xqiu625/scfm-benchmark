#!/bin/bash
#SBATCH --job-name=pooled_secondary_cv
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/pooled_secondary_cv_%j.log
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/pca200_pooled_secondary_cv.py \
  --which "$1" --lineage "$2"
