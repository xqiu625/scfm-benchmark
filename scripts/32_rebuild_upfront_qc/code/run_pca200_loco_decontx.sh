#!/bin/bash
#SBATCH --job-name=pca200_loco_decontx
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/pca200_loco_decontx_%j.log
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=30-00:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/pca200_loco_before_after_decontx.py --which decontx
