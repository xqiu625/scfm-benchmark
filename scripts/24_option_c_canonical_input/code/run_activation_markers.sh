#!/bin/bash
#SBATCH --job-name=activation_markers
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/activation_markers_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/activation_markers_%j.err

set -euo pipefail

echo "=== Platelet activation marker analysis (Option C) ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
echo "Python: $(which python)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/24_option_c_canonical_input/code/score_activation_markers.py \
    --embedding-dir /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c \
    --raw-input    /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad \
    --label-col    covid_severity \
    --fatal-tokens fatal \
    --output-dir   /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_activation_markers

echo "Completed: $(date)"
