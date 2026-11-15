#!/bin/bash
#SBATCH --job-name=abl_gf_rand
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_gf_random_dims_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/ablation_gf_random_dims_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=24:00:00

echo "============================================================"
echo "EMBEDDING ABLATION - Geneformer (random_dims only)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time limit: 24 hours"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

$PYTHON embedding_ablations.py --models Geneformer --tasks 3-class --ablations random_dims

echo ""
echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
