#!/bin/bash
#SBATCH --job-name=test_uce_mask
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/test_uce_mask_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/test_uce_mask_%j.err

set -euo pipefail
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
export MASKING_TEST_MAX_CELLS=40
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
PATH="$ENV_DIR/bin:$PATH" LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}" PYTHONNOUSERSITE=1 \
    python3 -u scripts/31_donor_level_revision/code/uce_masking_attribution.py
