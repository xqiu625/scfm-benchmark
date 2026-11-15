#!/bin/bash
#SBATCH --job-name=cosine_fatal
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/cosine_fatal_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/cosine_fatal_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Cosine similarity validation: cross-disease fatal directions ==="
echo "Start: $(date)"
echo "Host:  $(hostname)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "ERROR: $ENV_DIR/bin/python not found"
    exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
echo "Python: $(which python)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/24_option_c_canonical_input/code/cosine_fatal_direction_validation.py \
    --embedding-dir /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c \
    --output-dir    /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_cosine_validation \
    --label-col     covid_severity \
    --disease-col   Disease \
    --covid-survived  "mild,moderate,recovered,severe" \
    --sepsis-survived "severe"

echo "Completed: $(date)"
