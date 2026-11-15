#!/bin/bash
#SBATCH --job-name=19b_fsdp
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:a100:2
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/19b_fsdp_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/19b_fsdp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

echo "=================================================="
echo "ACC 19B: FSDP Distributed Fine-Tuning"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/19_model_optimization/fsdp

# Single GPU baseline
echo "Running single-GPU baseline..."
python scripts/19_model_optimization/code/fsdp_finetune.py --n-gpus 1 --epochs 15

# Multi-GPU FSDP (2 GPUs)
echo ""
echo "Running 2-GPU FSDP..."
torchrun --nproc_per_node=2 scripts/19_model_optimization/code/fsdp_finetune.py --n-gpus 2 --epochs 15

echo "Exit Code: $?"
echo "End Time: $(date)"
