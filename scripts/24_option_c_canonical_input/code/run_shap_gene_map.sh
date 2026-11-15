#!/bin/bash
#SBATCH --job-name=shap_gene_map
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/shap_gene_map_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/shap_gene_map_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== SHAP gene mapping (re-run after gene-map bug fix) ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/24_option_c_canonical_input/code/per_fm_shap_gene_map.py \
    --shap-dir      /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_shap \
    --embedding-dir /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c \
    --raw-input    /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad \
    --output-dir   /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_shap

echo "Completed: $(date)"
