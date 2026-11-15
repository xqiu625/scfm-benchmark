#!/bin/bash
#SBATCH --job-name=jobH_adaptcurve
#SBATCH --partition=epyc
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobH_step4_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobH_step4_%j.err
#SBATCH --requeue

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO"
mkdir -p "$REPO/scripts/31_donor_level_revision/logs"
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level/adaptation_curve_bootstrap_patient

python3 -u "$REPO/scripts/31_donor_level_revision/code/run_adaptation_curve_bootstrap_patient.py"
