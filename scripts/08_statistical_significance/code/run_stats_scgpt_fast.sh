#!/bin/bash
#SBATCH --job-name=stats_scgpt_fast
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stats_scgpt_fast_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stats_scgpt_fast_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=4:00:00

echo "============================================================"
echo "OPTIMIZED STATISTICAL SIGNIFICANCE - scGPT"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time limit: 4 hours"
echo "Memory: 64GB"
echo "CPUs: 16 (parallel bootstrap)"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

$PYTHON statistical_significance_fast.py \
    --models scGPT \
    --tasks binary 3-class \
    --n-bootstrap 100 \
    --pca-dims 256

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
