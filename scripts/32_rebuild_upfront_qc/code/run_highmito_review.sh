#!/bin/bash
#SBATCH --job-name=highmito_review
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/highmito_review_%j.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/adhoc_highmito_platelet_identity.py
