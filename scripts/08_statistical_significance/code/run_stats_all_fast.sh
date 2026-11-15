#!/bin/bash
#SBATCH --job-name=stats_all_fast
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stats_all_fast_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stats_all_fast_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=8:00:00

echo "============================================================"
echo "OPTIMIZED STATISTICAL SIGNIFICANCE - ALL MODELS"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time limit: 8 hours"
echo "Memory: 128GB"
echo "CPUs: 16 (parallel bootstrap)"
echo "============================================================"
echo ""
echo "Optimizations:"
echo "  - Parallel bootstrap (16 cores)"
echo "  - Fast lbfgs solver (vs saga)"
echo "  - PCA reduction to 256 dims"
echo "  - 100 bootstrap iterations"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# Run ALL models at once (most efficient - loads data once)
$PYTHON statistical_significance_fast.py \
    --models STATE UCE scGPT Geneformer \
    --tasks binary 3-class \
    --n-bootstrap 100 \
    --pca-dims 256

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
