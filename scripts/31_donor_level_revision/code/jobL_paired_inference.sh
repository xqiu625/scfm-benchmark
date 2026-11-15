#!/bin/bash
#SBATCH --job-name=jobL_paired
#SBATCH --partition=epyc
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobL_paired_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobL_paired_%j.err
#SBATCH --requeue

set -euo pipefail
export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO/scripts/31_donor_level_revision/code"
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/paired_inference

python3 -u jobL_paired_inference.py
