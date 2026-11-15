#!/bin/bash
#SBATCH --job-name=21_cvae
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/21_cvae_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/21_cvae_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC 21: Conditional VAE for Severity Transition Generation
# ============================================================================

echo "=================================================="
echo "ACC 21: Conditional VAE"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/21_conditional_generation/figures

echo "Python: $(which python)"

# Step 1: Train cVAE and run transition experiment
python scripts/21_conditional_generation/code/conditional_vae.py \
    --epochs 200 \
    --beta 0.5 \
    --latent-dim 128 \
    --batch-size 256

# Step 2: Evaluate generated embeddings
echo ""
echo "Evaluating generated embeddings..."
python scripts/21_conditional_generation/code/evaluate_generation.py

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

exit $EXIT_CODE
