#!/bin/bash
#SBATCH --job-name=manuscript_umap
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/manuscript_umap_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/manuscript_umap_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=short

# Manuscript UMAP Figure Generation
# Generates 6-model x 3-scheme UMAP visualizations

echo "=========================================="
echo "MANUSCRIPT UMAP FIGURES"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /bigdata/godziklab/shared/Xinru/302006/state_env

# Navigate to code directory
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/12_umap_visualization/code

# Run the manuscript figure generation script
python manuscript_umap_figures.py

echo ""
echo "=========================================="
echo "COMPLETE"
echo "End time: $(date)"
echo "=========================================="
