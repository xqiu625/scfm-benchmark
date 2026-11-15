#!/bin/bash
#SBATCH --job-name=build_shards
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/build_shards_%j.log
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00

source activate scgpt_finetune
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/build_extraction_shards.py \
  --which "$1" --n-shards "$2"
