#!/bin/bash
#SBATCH --job-name=0006_subspace
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/subspace_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/subspace_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL
set -euo pipefail
echo "=== cross-disease subspace alignment (roadmap H + N) ==="; echo "Start: $(date)"
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"; export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"; export PYTHONNOUSERSITE=1
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python scripts/30_review_response_stats/code/subspace_alignment.py
echo "Completed: $(date)"
