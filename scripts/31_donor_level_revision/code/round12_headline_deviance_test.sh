#!/bin/bash
#SBATCH --job-name=r12_headdev
#SBATCH --partition=epyc
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-3
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/r12_headdev_%A_%a.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/r12_headdev_%A_%a.err

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1
REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO/scripts/31_donor_level_revision/code"
mkdir -p "$REPO/scripts/31_donor_level_revision/logs"
python3 -u round12_headline_deviance_test.py "$SLURM_ARRAY_TASK_ID"
