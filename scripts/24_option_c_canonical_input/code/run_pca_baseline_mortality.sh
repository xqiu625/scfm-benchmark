#!/bin/bash
#SBATCH --job-name=pca_baseline
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/pca_baseline_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/pca_baseline_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== PCA baseline on raw gene expression (mortality binary) ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/24_option_c_canonical_input/code/pca_baseline_mortality.py \
    --input  /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad \
    --output /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_pca_baseline/pca_baseline_summary.csv

echo "Completed: $(date)"
