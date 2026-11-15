#!/bin/bash
#SBATCH --job-name=0006_seed_stability
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/seed_stability_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/seed_stability_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== Ranking seed-stability (8 CV seeds x 5 FMs, within-disease COVID+Sepsis LR) ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
echo "Python: $(which python)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/30_review_response_stats/code/seed_stability.py

echo "Completed: $(date)"
