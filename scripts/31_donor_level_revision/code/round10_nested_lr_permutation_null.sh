#!/bin/bash
#SBATCH --job-name=r10_nestedLR
#SBATCH --partition=epyc
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-5
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/r10_nestedLR_%A_%a.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/r10_nestedLR_%A_%a.err

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1
REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO/scripts/31_donor_level_revision/code"
mkdir -p "$REPO/scripts/31_donor_level_revision/logs"
python3 -u round10_nested_lr_permutation_null.py "$SLURM_ARRAY_TASK_ID"
