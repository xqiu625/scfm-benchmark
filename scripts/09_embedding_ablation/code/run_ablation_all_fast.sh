#!/bin/bash
#SBATCH --job-name=ablation_all_fast
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_all_fast_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_all_fast_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=4:00:00

echo "============================================================"
echo "OPTIMIZED EMBEDDING ABLATION - ALL MODELS"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time limit: 4 hours"
echo "Memory: 128GB"
echo "CPUs: 8 (parallel CV folds)"
echo "============================================================"
echo ""
echo "Optimizations:"
echo "  - Parallel CV folds (8 cores)"
echo "  - Fast lbfgs solver (vs saga)"
echo "  - All 4 ablation types"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# Run ALL models with ALL ablations
$PYTHON embedding_ablations_fast.py \
    --models STATE UCE scGPT Geneformer \
    --tasks binary 3-class \
    --ablations dim_reduction remove_top_pca normalization random_dims

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
