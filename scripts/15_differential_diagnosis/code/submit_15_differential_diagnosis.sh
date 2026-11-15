#!/bin/bash
#SBATCH --job-name=15_diff_diag
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/15_differential_diagnosis_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/15_differential_diagnosis_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC-15: Differential Diagnosis of Inflammatory Conditions (5-class)
# ============================================================================
# Task 2 of Two-Task Clinical Framework
# 5 Classes: Healthy, COVID-19, Sepsis, SLE, Pneumonia
# Excluded: Influenza (n=33), Non_COVID (catch-all)
# ============================================================================

echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "=================================================="

# Setup environment
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Use 302006_env_clean (recommended for analysis)
export ENV_PATH="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"

# Disable user site-packages to avoid conflicts with ~/.local
export PYTHONNOUSERSITE=1

# Verify environment
echo "Python: $(which python)"
python -c "import scanpy; print('Scanpy version:', scanpy.__version__)"
python -c "import sklearn; print('Sklearn version:', sklearn.__version__)"

# Base directory (git repo with code)
BASE_DIR="/bigdata/godziklab/shared/Xinru/302006/302006_updategit"
# Data directory (embeddings, results)
DATA_DIR="/bigdata/godziklab/shared/Xinru/302006"

# Run differential diagnosis benchmark
echo ""
echo "Running 5-class Differential Diagnosis Benchmark..."
echo "Classes: Healthy, COVID-19, Sepsis, SLE, Pneumonia"
echo ""

python scripts/15_differential_diagnosis/code/run_15_differential_diagnosis.py \
    --base-dir "$DATA_DIR" \
    --cv-folds 5 \
    --n-jobs 16

echo ""
echo "=================================================="
echo "End Time: $(date)"
echo "=================================================="

# Copy results to analyses folder structure
echo ""
echo "Syncing results to analyses folder..."
RESULTS_DIR="$DATA_DIR/05_RESULTS/15_differential_diagnosis"
ANALYSES_DIR="$BASE_DIR/scripts/15_differential_diagnosis"

if [ -d "$RESULTS_DIR" ]; then
    cp -r "$RESULTS_DIR"/* "$ANALYSES_DIR/results/" 2>/dev/null || true
    echo "Results copied to $ANALYSES_DIR/results/"
fi

echo "Job complete!"
