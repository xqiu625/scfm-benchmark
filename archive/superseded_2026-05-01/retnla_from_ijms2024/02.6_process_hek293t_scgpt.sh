#!/bin/bash
#SBATCH --job-name=hek293t_scgpt
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=2
#SBATCH --mem=400G
#SBATCH --time=23:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/logs/hek293t_scgpt_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/logs/hek293t_scgpt_%j.err

cd /bigdata/godziklab/shared/Xinru/302006
mkdir -p logs retnla

module load singularity

singularity exec \
    --bind /bigdata/godziklab/shared/Xinru:/bigdata/godziklab/shared/Xinru \
    scgpt_0.1.7.sif \
    python process_hek293t_scgpt.py