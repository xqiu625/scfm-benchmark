#!/bin/bash
#SBATCH --job-name=r6_purity_audit
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/r6_purity_audit_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/r6_purity_audit_%j.err

source ~/.bashrc
conda activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
python /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code/r6_purity_and_gene_audit.py
