#!/bin/bash
#SBATCH --job-name=debug_umap
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/debug_umap_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/debug_umap_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=short

echo "=========================================="
echo "DEBUG UMAP - Individual Models"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

eval "$(conda shell.bash hook)"
conda activate /bigdata/godziklab/shared/Xinru/302006/state_env

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/12_umap_visualization/code/selected_code
python debug_umap_individual.py

echo "End time: $(date)"
