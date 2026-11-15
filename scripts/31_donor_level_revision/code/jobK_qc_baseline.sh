#!/bin/bash
#SBATCH --job-name=jobK_qc
#SBATCH --partition=epyc
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobK_qc_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobK_qc_%j.err
#SBATCH --requeue

set -euo pipefail

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO/scripts/31_donor_level_revision/code"

mkdir -p "$REPO/scripts/31_donor_level_revision/logs"
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/qc_baseline

python3 -u jobK_qc_baseline.py
