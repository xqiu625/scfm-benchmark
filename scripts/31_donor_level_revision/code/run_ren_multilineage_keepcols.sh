#!/bin/bash
#SBATCH --job-name=ren_ml_keepcols
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/ren_ml_keepcols_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/ren_ml_keepcols_%j.err

source ~/.bashrc
conda activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
python /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code/ren_multilineage_keepcols.py
