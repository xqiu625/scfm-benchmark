#!/bin/bash
#SBATCH --job-name=build_qiu_input
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/build_qiu_input_%j.log
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/build_qiu_decontx_input.py \
  --raw-dir /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted \
  --inventory /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/channel_inventories/qiu_gse167363_channel_inventory.csv \
  --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/qiu_decontx_input
