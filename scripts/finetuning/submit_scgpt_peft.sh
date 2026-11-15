#!/bin/bash
#SBATCH --job-name=scgpt_peft
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scgpt_peft_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scgpt_peft_%j.log
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --account=godziklab
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256gb
#SBATCH --time=72:00:00

echo "============================================================"
echo "scGPT PEFT FINE-TUNING (scPEFT Methodology)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "============================================================"

# Activate environment
source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
export PYTHONNOUSERSITE=1
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH
export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/302006/scpeft_env/lib:$LD_LIBRARY_PATH

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

echo ""
echo "GPU Status:"
nvidia-smi

# Run with different PEFT types
echo ""
echo "============================================================"
echo "scGPT LoRA Binary Classification"
echo "============================================================"
python scgpt_peft_finetune.py --task binary --peft LORA --epochs 20 --resume

echo ""
echo "============================================================"
echo "scGPT LoRA 3-Class Classification"
echo "============================================================"
python scgpt_peft_finetune.py --task 3class --peft LORA --epochs 20 --resume

echo ""
echo "============================================================"
echo "scGPT LoRA 6-Class Classification"
echo "============================================================"
python scgpt_peft_finetune.py --task 6class --peft LORA --epochs 20 --resume

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
