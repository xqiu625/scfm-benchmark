#!/bin/bash
#SBATCH --job-name=finetune_v2
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/finetune_v2_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/finetune_v2_%j.log
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=48:00:00

echo "============================================================"
echo "ENHANCED FINE-TUNING V2 (Focal Loss + Optuna + Ensemble)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# Install optuna if not available
$PYTHON -c "import optuna" 2>/dev/null || pip install optuna

# Run all tasks with optimization and ensemble
for TASK in binary 3class 6class; do
    echo ""
    echo ">>> Running ${TASK} task with optimization..."
    $PYTHON unified_finetune_v2.py --model all --task $TASK --strategy all --optimize --ensemble
done

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
