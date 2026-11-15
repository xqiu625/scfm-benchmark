#!/bin/bash
#SBATCH --job-name=build_liu_input
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/build_liu_input_%j.log
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=220G
#SBATCH --time=04:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/build_liu_decontx_input.py   --singlets-dir /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE161918_liu/singlets   --h5-dir /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE161918_liu/rna_h5   --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/liu_decontx_input
