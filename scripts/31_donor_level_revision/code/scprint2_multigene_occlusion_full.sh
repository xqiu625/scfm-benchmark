#!/bin/bash
#SBATCH --job-name=scprint2_occl_full
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/scprint2_occl_full_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/scprint2_occl_full_%j.err

set -euo pipefail
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface

python3 -u scripts/31_donor_level_revision/code/scprint2_multigene_occlusion.py
