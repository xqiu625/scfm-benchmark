#!/bin/bash
#SBATCH --job-name=pseudobulk_en_c1
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/pseudobulk_en_c1_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/pseudobulk_en_c1_%j.err
export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python3 -u scripts/31_donor_level_revision/code/pseudobulk_elasticnet_patient_oof.py --seeds 42,43,44,45,46
