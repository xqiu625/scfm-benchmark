#!/bin/bash
#SBATCH --job-name=lora_finetune
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/lora_finetune_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/lora_finetune_%j.log
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00

echo "============================================================"
echo "LoRA FINE-TUNING - Geneformer"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# Run LoRA fine-tuning for all tasks
for TASK in binary 3class 6class; do
    echo ""
    echo ">>> LoRA Fine-tuning: ${TASK}..."
    $PYTHON geneformer_lora_finetune.py --task $TASK --lora-r 8 --lora-alpha 16
done

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
