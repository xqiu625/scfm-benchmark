#!/bin/bash
#SBATCH --job-name=17a_geneformer_pretrain
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17a_geneformer_pretrain_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17a_geneformer_pretrain_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC 17A: Continue Pre-Training Geneformer on Platelet Data
# ============================================================================
# Domain-adapt Geneformer-12L-30M checkpoint using MLM objective
# on 47,808 platelet cells (50 epochs, ~3-5 hours on A100)
# ============================================================================

echo "=================================================="
echo "ACC 17A: Geneformer Continued Pre-Training"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Start Time: $(date)"
echo "=================================================="

# Setup environment
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Activate scpeft_env (has transformers + torch)
export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/scpeft_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH=${ENV_PATH}/lib:$LD_LIBRARY_PATH

# Create output directories
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/17_domain_adapted_fm/geneformer_adapted
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/checkpoints/17_geneformer_adapted
mkdir -p /bigdata/godziklab/shared/Xinru/302006/06_LOGS

echo "Python: $(which python)"
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')"
python -c "from transformers import BertForMaskedLM; print('transformers OK')"

# Run continued pre-training
python scripts/17_domain_adapted_fm/code/geneformer_continue_pretrain.py \
    --epochs 50 \
    --batch-size 8 \
    --lr 1e-5 \
    --mask-ratio 0.15

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

exit $EXIT_CODE
