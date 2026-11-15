#!/bin/bash
#SBATCH --job-name=build_stephenson_input
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/build_stephenson_input_%j.log
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=220G
#SBATCH --time=08:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/build_stephenson_decontx_input.py \
  --h5ad /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/Stephenson_EMTAB10026/covid_portal_210320_with_raw.h5ad \
  --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/stephenson_decontx_input
