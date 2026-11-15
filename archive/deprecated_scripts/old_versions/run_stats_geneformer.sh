#!/bin/bash
#SBATCH --job-name=stats_gf
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stats_geneformer_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stats_geneformer_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256gb
#SBATCH --time=72:00:00

echo "============================================================"
echo "STATISTICAL SIGNIFICANCE - Geneformer (72h, 256GB)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time limit: 72 hours"
echo "Memory: 256GB"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

$PYTHON statistical_significance.py --models Geneformer --n-bootstrap 500 --tasks binary 3-class

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
