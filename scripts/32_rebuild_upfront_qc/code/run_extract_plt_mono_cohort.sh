#!/bin/bash
#SBATCH --job-name=extract_plt_mono
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/extract_plt_mono_%j.log
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/assemble_platelet_monocyte_benchmark_input.py \
  --ckpt-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_checkpoints \
  --cohort "$1"
