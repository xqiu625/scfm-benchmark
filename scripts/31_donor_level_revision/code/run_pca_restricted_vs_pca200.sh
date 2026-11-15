#!/bin/bash
#SBATCH --job-name=pca_restricted_cmp
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/pca_restricted_cmp_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/pca_restricted_cmp_%j.err

source ~/.bashrc
conda activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code
python pca_restricted_vs_pca200.py
