#!/bin/bash
#SBATCH --job-name=hvg_marker_full
#SBATCH --partition=epyc
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/hvg_marker_full_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/hvg_marker_full_%j.err
export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python3 -u scripts/31_donor_level_revision/code/hvg_marker_baselines_patient_oof_v2.py --n-seeds 10
