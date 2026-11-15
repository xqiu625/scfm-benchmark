#!/bin/bash
#SBATCH --job-name=jobS_integrate
#SBATCH --partition=highmem
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_integrate_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_integrate_%j.err

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO"
mkdir -p "$REPO/scripts/31_donor_level_revision/logs"

python3 -u "$REPO/scripts/31_donor_level_revision/code/integrate_chongqing_sepsis.py"
