#!/bin/bash
#SBATCH --job-name=mouse_scgpt
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=400G  
#SBATCH --time=12:00:00  
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/logs/mouse_scgpt_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/logs/mouse_scgpt_%j.err

cd /bigdata/godziklab/shared/Xinru/302006
mkdir -p logs retnla

# Set cache directories
export NUMBA_CACHE_DIR=/bigdata/godziklab/shared/Xinru/temp_cache
export PYTHONUSERBASE=/bigdata/godziklab/shared/Xinru/temp_python
mkdir -p /bigdata/godziklab/shared/Xinru/temp_cache
mkdir -p /bigdata/godziklab/shared/Xinru/temp_python

module load singularity

echo "Processing Mouse RETNLA data with scGPT..."

# Run in Singularity container
singularity exec \
    --bind /bigdata/godziklab/shared/Xinru:/bigdata/godziklab/shared/Xinru \
    scgpt_0.1.7.sif \
    python process_mouse_scgpt.py

echo "Mouse processing completed!"
ls -la /bigdata/godziklab/shared/Xinru/302006/retnla/mouse_*