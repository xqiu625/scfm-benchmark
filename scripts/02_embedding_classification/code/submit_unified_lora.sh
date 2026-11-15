#!/bin/bash
#SBATCH --job-name=unified_lora
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/unified_lora_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/unified_lora_%j.log
#SBATCH --partition=short_gpu
#SBATCH --account=godziklab
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00

echo "============================================================"
echo "UNIFIED LoRA FINE-TUNING - All Models"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# Run LoRA fine-tuning for all models, all tasks
for TASK in binary 3class 6class; do
    echo ""
    echo ">>> LoRA Fine-tuning: All models - ${TASK}..."
    $PYTHON unified_lora_finetune.py --model all --task $TASK --lora-r 8 --lora-alpha 16
done

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
