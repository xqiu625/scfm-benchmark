#!/bin/bash
#SBATCH --job-name=ablation_scgpt_3class
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_scgpt_3class_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_scgpt_3class_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=48:00:00

echo "============================================================"
echo "EMBEDDING ABLATION - scGPT (3-class only, extended time)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time limit: 48 hours"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# Only run 3-class task (binary already completed in job 22072392)
$PYTHON embedding_ablations.py --models scGPT --tasks 3-class

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
