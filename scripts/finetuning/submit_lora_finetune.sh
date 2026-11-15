#!/bin/bash
#SBATCH --job-name=state_lora
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/state_lora_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/state_lora_%j.log
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --account=godziklab
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=6:00:00

echo "============================================================"
echo "STATE SE-600M LoRA FINE-TUNING (PEFT)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "============================================================"

# Activate environment
source activate /bigdata/godziklab/shared/Xinru/302006/state_env
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

echo ""
echo "GPU Status:"
nvidia-smi

# LoRA configuration
# r=16: ~2% trainable parameters
# max-genes=512: More context than full fine-tuning (was 256)

echo ""
echo "============================================================"
echo "LoRA Binary Classification (r=16, max_genes=512)"
echo "============================================================"
python state_lora_finetune.py --task binary --epochs 20 --max-genes 512 --lora-r 16

echo ""
echo "============================================================"
echo "LoRA 3-Class Classification"
echo "============================================================"
python state_lora_finetune.py --task 3class --epochs 20 --max-genes 512 --lora-r 16

echo ""
echo "============================================================"
echo "LoRA 6-Class Classification"
echo "============================================================"
python state_lora_finetune.py --task 6class --epochs 20 --max-genes 512 --lora-r 16

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
