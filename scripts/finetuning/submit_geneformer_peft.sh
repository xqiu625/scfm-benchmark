#!/bin/bash
#SBATCH --job-name=geneformer_peft
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/geneformer_peft_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/geneformer_peft_%j.log
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
echo "GENEFORMER PEFT FINE-TUNING (scPEFT Methodology)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "============================================================"

# Activate Geneformer environment
source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
export PYTHONNOUSERSITE=1
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH
export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/302006/scpeft_env/lib:$LD_LIBRARY_PATH

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

echo ""
echo "GPU Status:"
nvidia-smi

# Run with LoRA
echo ""
echo "============================================================"
echo "Geneformer LoRA Binary Classification"
echo "============================================================"
python geneformer_peft_finetune.py --task binary --epochs 25 --lr 0.0008 --resume

echo ""
echo "============================================================"
echo "Geneformer LoRA 3-Class Classification"
echo "============================================================"
python geneformer_peft_finetune.py --task 3class --epochs 25 --lr 0.0008 --resume

echo ""
echo "============================================================"
echo "Geneformer LoRA 6-Class Classification"
echo "============================================================"
python geneformer_peft_finetune.py --task 6class --epochs 25 --lr 0.0008 --resume

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
