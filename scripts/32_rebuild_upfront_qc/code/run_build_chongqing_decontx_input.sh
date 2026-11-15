#!/bin/bash
#SBATCH --job-name=build_chongqing_input
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/build_chongqing_input_%j.log
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=08:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/build_chongqing_decontx_input.py \
  --raw-dir /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/chongqing_sepsis_GSE279452/extracted_full \
  --labels /bigdata/godziklab/shared/Xinru/302006/302006_updategit/reports/manuscript/revision/literature_audit/canonical_donor_labels.csv \
  --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/chongqing_decontx_input
