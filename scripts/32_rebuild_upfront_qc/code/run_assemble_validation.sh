#!/bin/bash
#SBATCH --job-name=assemble_validation
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/assemble_validation_%j.log
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --time=04:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/assemble_within_covid_validation.py   --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/within_covid_validation
