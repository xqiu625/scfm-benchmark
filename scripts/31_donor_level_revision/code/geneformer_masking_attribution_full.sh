#!/bin/bash
#SBATCH --job-name=geneformer_mask_full
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/geneformer_mask_full_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/geneformer_mask_full_%j.err

set -euo pipefail
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scpeft_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python3 -u scripts/31_donor_level_revision/code/geneformer_masking_attribution.py
