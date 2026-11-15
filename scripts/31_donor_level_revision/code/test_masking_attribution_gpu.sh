#!/bin/bash
#SBATCH --job-name=test_masking_attr_g
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/test_masking_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/test_masking_%j.err

set -euo pipefail
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
export MASKING_TEST_MAX_CELLS=40

echo "=== Geneformer masking attribution (small test) ==="
ENV_DIR_GF=/bigdata/godziklab/shared/Xinru/302006/scpeft_env
PATH="$ENV_DIR_GF/bin:$PATH" LD_LIBRARY_PATH="$ENV_DIR_GF/lib:${LD_LIBRARY_PATH:-}" PYTHONNOUSERSITE=1 \
    python3 -u scripts/31_donor_level_revision/code/geneformer_masking_attribution.py

echo ""
echo "=== STATE masking attribution (small test) ==="
ENV_DIR_ST=/bigdata/godziklab/shared/Xinru/302006/state_env
PATH="$ENV_DIR_ST/bin:$PATH" LD_LIBRARY_PATH="$ENV_DIR_ST/lib:${LD_LIBRARY_PATH:-}" PYTHONNOUSERSITE=1 \
    python3 -u scripts/31_donor_level_revision/code/state_masking_attribution.py
