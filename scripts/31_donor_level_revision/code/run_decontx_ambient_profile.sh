#!/bin/bash
#SBATCH --job-name=decontx_ambient
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/decontx_ambient_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/decontx_ambient_%j.err

module load R/4.5.2
export R_LIBS=/bigdata/godziklab/shared/Xinru/302006/r_libs_decontx:$R_LIBS
Rscript /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code/decontx_ambient_profile.R
