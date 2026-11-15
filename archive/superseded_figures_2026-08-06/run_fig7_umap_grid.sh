#!/bin/bash
#SBATCH --job-name=fig7_umap
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/fig7_umap_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/fig7_umap_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== Figure 7 — per-FM mouse UMAP grid ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/26_figures/code/generate_fig7_umap_grid.py \
    --output-dir reports/manuscript/figures/fig7

echo "Completed: $(date)"
