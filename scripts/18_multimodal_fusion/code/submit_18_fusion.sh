#!/bin/bash
#SBATCH --job-name=18_fusion
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/18_fusion_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/18_fusion_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC 18: Multimodal Fusion of Foundation Model Embeddings
# ============================================================================
# 4 strategies: early, late, attention-based, expression+embedding
# Attention fusion requires GPU; others are CPU-only
# ============================================================================

echo "=================================================="
echo "ACC 18: Multimodal Fusion"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/18_multimodal_fusion/figures

echo "Python: $(which python)"
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

python scripts/18_multimodal_fusion/code/multimodal_fusion.py \
    --strategy all \
    --attention-epochs 100 \
    --common-dim 256

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

exit $EXIT_CODE
