#!/bin/bash
#SBATCH --job-name=cluster_sep
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/cluster_separation_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/cluster_separation_%j.log
#SBATCH --partition=short_gpu
#SBATCH --account=godziklab
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=1:00:00

echo "============================================================"
echo "CLUSTER SEPARATION ANALYSIS"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/12_umap_visualization/code/selected_code

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

$PYTHON cluster_separation_analysis.py

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
