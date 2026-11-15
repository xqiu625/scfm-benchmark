#!/bin/bash
#SBATCH --job-name=state_true_ft
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/state_true_finetune_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/state_true_finetune_%j.log
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --account=godziklab
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=12:00:00

echo "============================================================"
echo "STATE SE-600M TRUE END-TO-END FINE-TUNING"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "============================================================"

# Activate STATE environment
source activate /bigdata/godziklab/shared/Xinru/302006/state_env
export PYTHONNOUSERSITE=1

# Set working directory
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

# Check GPU
echo ""
echo "GPU Status:"
nvidia-smi

# Check that required files exist
echo ""
echo "Checking required files..."
SE_MODEL="/bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M/se600m_epoch16.ckpt"
PROTEIN_EMB="/bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M/protein_embeddings.pt"

if [ -f "$SE_MODEL" ]; then
    echo "SE-600M checkpoint: OK"
else
    echo "ERROR: SE-600M checkpoint not found at $SE_MODEL"
    exit 1
fi

if [ -f "$PROTEIN_EMB" ]; then
    echo "Protein embeddings: OK ($(du -h $PROTEIN_EMB | cut -f1))"
else
    echo "ERROR: Protein embeddings not found at $PROTEIN_EMB"
    exit 1
fi

# Run TRUE fine-tuning (gradients through SE-600M)
echo ""
echo "============================================================"
echo "Running TRUE Fine-Tuning (Binary Classification)"
echo "Backpropagating through 603M parameter SE encoder"
echo "============================================================"
python state_true_finetune_v3.py --task binary --epochs 15 --batch-size 16 --lr 1e-5

echo ""
echo "============================================================"
echo "Running TRUE Fine-Tuning (3-Class Classification)"
echo "============================================================"
python state_true_finetune_v3.py --task 3class --epochs 15 --batch-size 16 --lr 1e-5

echo ""
echo "============================================================"
echo "Running TRUE Fine-Tuning (6-Class Classification)"
echo "============================================================"
python state_true_finetune_v3.py --task 6class --epochs 15 --batch-size 16 --lr 1e-5

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
