#!/bin/bash
#SBATCH --job-name=subsample_plt_mono
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/subsample_plt_mono_%j.log
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/subsample_platelet_monocyte_cells.py \
  --which "$1"
