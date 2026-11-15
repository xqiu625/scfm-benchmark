#!/bin/bash
#SBATCH --job-name=abl_uce_fast
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_uce_fast_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_uce_fast_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00

echo "============================================================"
echo "OPTIMIZED EMBEDDING ABLATION - UCE"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time limit: 2 hours"
echo "Memory: 64GB"
echo "CPUs: 8 (parallel CV folds)"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

$PYTHON embedding_ablations_fast.py \
    --models UCE \
    --tasks binary 3-class \
    --ablations dim_reduction remove_top_pca normalization random_dims

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
