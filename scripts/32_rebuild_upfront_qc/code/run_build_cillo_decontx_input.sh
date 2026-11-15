#!/bin/bash
#SBATCH --job-name=build_cillo_input
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/build_cillo_input_%j.log
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=08:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/build_cillo_decontx_input.py \
  --raw-dir /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE180578_cillo \
  --h5ad /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE180578_cillo/GSE180578_cillo_covid19_study_aggregrated_annotated_data.h5ad \
  --manifest /bigdata/godziklab/shared/Xinru/302006/302006_updategit/reports/manuscript/revision/literature_audit/T4_outcome_severity_verification/cillo_pbmc_manifest.csv \
  --inventory /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/channel_inventories/cillo_gse180578_channel_inventory.csv \
  --work-dir /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE180578_cillo/_extract_scratch \
  --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/cillo_decontx_input
