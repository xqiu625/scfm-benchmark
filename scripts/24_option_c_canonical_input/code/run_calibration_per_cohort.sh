#!/bin/bash
#SBATCH --job-name=0006_calibration_cohort
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/calibration_cohort_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/calibration_cohort_%j.err

set -euo pipefail
echo "=== Calibration + per-cohort AUC (CZI workshop alignment) ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/24_option_c_canonical_input/code/calibration_per_cohort.py \
    --embedding-dir /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c \
    --label-col covid_severity \
    --group-col Sample_ID \
    --cohort-col Data_NO \
    --output-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_calibration

echo "Completed: $(date)"
