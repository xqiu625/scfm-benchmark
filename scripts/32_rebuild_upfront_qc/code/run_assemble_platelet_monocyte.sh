#!/bin/bash
#SBATCH --job-name=assemble_plt_mono
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/assemble_plt_mono_%j.log
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=220G
#SBATCH --time=08:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/assemble_platelet_monocyte_benchmark_input.py \
  --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input
