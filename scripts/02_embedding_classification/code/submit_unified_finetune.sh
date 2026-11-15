#!/bin/bash
#SBATCH --job-name=unified_finetune
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/unified_finetune_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/unified_finetune_%j.log
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00

echo "============================================================"
echo "UNIFIED FINE-TUNING BENCHMARK"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# Run binary task (all models, all strategies)
echo ""
echo ">>> Running BINARY task..."
$PYTHON unified_finetune.py --model all --task binary --strategy all

# Run 3-class task (all models, all strategies)
echo ""
echo ">>> Running 3-CLASS task..."
$PYTHON unified_finetune.py --model all --task 3class --strategy all

# Run 6-class task (all models, all strategies)
echo ""
echo ">>> Running 6-CLASS task..."
$PYTHON unified_finetune.py --model all --task 6class --strategy all

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
