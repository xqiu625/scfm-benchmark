#!/bin/bash
#SBATCH --job-name=jobS_adaptcurve
#SBATCH --partition=epyc
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_adaptcurve_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_adaptcurve_%j.err

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1
REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO/scripts/31_donor_level_revision/code"
mkdir -p "$REPO/scripts/31_donor_level_revision/logs"
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/adaptation_curve_bootstrap_chongqing
python3 -u run_adaptation_curve_chongqing.py
