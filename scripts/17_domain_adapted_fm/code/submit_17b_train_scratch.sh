#!/bin/bash
#SBATCH --job-name=17b_platelet_tf
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17b_platelet_tf_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17b_platelet_tf_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC 17B: Train Platelet Transformer from Scratch
# ============================================================================
# 6-layer, 384-dim, 6-head transformer (~20M params)
# Pre-train on platelet (48K) + X-Atlas (500K) cells with MGEM
# ============================================================================

echo "=================================================="
echo "ACC 17B: Platelet Transformer Training from Scratch"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Activate state_env (PyTorch 2.5.1)
export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/17_domain_adapted_fm/platelet_transformer
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/checkpoints/17_platelet_transformer

echo "Python: $(which python)"
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

# Run training
python scripts/17_domain_adapted_fm/code/platelet_transformer.py \
    --epochs 100 \
    --batch-size 64 \
    --lr 3e-4 \
    --d-model 384 \
    --n-layers 6 \
    --n-heads 6 \
    --no-xatlas

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

exit $EXIT_CODE
