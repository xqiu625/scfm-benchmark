#!/bin/bash
#SBATCH --job-name=decontx_reaudit
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/decontx_reaudit_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/decontx_reaudit_%j.err

source ~/.bashrc
conda activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
python /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code/decontx_reaudit.py
