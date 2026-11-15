#!/bin/bash
#SBATCH --job-name=jobS2_lobo_pool
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS2_lobo_pool_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS2_lobo_pool_%j.err

set -euo pipefail
echo "=== Round-4 primary-test correction, item 5: concatenated-pooled LOBO ==="
date

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python scripts/31_donor_level_revision/code/jobS2_pooled_lobo.py

echo "DONE_JOBS2_LOBO_POOL"
date
