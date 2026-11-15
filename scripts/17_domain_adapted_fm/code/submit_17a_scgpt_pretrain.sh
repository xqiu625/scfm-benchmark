#!/bin/bash
#SBATCH --job-name=17a_scgpt_pretrain
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17a_scgpt_pretrain_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17a_scgpt_pretrain_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC 17A: Continue Pre-Training scGPT on Platelet Data
# ============================================================================
# Domain-adapt scGPT whole-human checkpoint using MGEM objective
# on 47,808 platelet cells (50 epochs, ~3-5 hours on A100)
# ============================================================================

echo "=================================================="
echo "ACC 17A: scGPT Continued Pre-Training"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Start Time: $(date)"
echo "=================================================="

# Setup environment
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Activate scpeft_env for scGPT
source ~/.bashrc
conda activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH

# Disable user site-packages
export PYTHONNOUSERSITE=1

# Fix CXXABI_1.3.15 error: ensure conda's libstdc++ is found before system's
export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/302006/scpeft_env/lib:$LD_LIBRARY_PATH

# Create output directories
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/17_domain_adapted_fm/scgpt_adapted
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/checkpoints/17_scgpt_adapted
mkdir -p /bigdata/godziklab/shared/Xinru/302006/06_LOGS

echo "Python: $(which python)"
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

# Run continued pre-training
python scripts/17_domain_adapted_fm/code/scgpt_continue_pretrain.py \
    --epochs 50 \
    --batch-size 32 \
    --lr 5e-5 \
    --mask-ratio 0.15

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

exit $EXIT_CODE
