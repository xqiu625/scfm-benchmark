#!/bin/bash
#SBATCH --job-name=batch_effect
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/batch_effect_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/batch_effect_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch

# Batch Effect Analysis: 6 Models Comparison
# Metrics: ARI, NMI, ASW, kBET, iLISI

echo "=========================================="
echo "BATCH EFFECT ANALYSIS"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /bigdata/godziklab/shared/Xinru/302006/state_env

# Navigate to code directory
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/12_umap_visualization/code/selected_code

# Run the batch effect analysis script
python batch_effect_analysis.py

echo ""
echo "=========================================="
echo "COMPLETE"
echo "End time: $(date)"
echo "=========================================="
