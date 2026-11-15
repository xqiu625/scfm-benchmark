#!/bin/bash
#SBATCH --job-name=check_zero
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/check_zero_%j.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/adhoc_check_zero.py
