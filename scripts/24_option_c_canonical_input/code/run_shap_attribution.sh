#!/bin/bash
#SBATCH --job-name=shap_per_fm
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/shap_per_fm_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/shap_per_fm_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== Per-FM SHAP attribution on mortality LR (Option C 6-FM) ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1

# shap dependency (idempotent)
python -c "import shap" 2>/dev/null || \
    uv pip install --python "$ENV_DIR/bin/python" shap

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/24_option_c_canonical_input/code/per_fm_shap_attribution.py \
    --embedding-dir /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c \
    --raw-input    /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad \
    --output-dir   /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_shap

echo "Completed: $(date)"
