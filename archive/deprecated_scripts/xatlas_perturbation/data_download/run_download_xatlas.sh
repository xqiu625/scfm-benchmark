#!/bin/bash
#SBATCH --job-name=dl_xatlas
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/download_xatlas_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/download_xatlas_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256gb
#SBATCH --time=24:00:00

echo "============================================================"
echo "Download X-Atlas/Orion Dataset from Hugging Face"
echo "Started at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "============================================================"

# IMPORTANT: Set HuggingFace cache to bigdata (avoid home directory quota)
export HF_HOME="/bigdata/godziklab/shared/Xinru/302006/temp/hf_cache"
export HF_DATASETS_CACHE="/bigdata/godziklab/shared/Xinru/302006/temp/hf_cache/datasets"
export TRANSFORMERS_CACHE="/bigdata/godziklab/shared/Xinru/302006/temp/hf_cache/transformers"
mkdir -p $HF_HOME $HF_DATASETS_CACHE $TRANSFORMERS_CACHE

echo "HF_HOME: $HF_HOME"
echo "HF_DATASETS_CACHE: $HF_DATASETS_CACHE"

# Activate environment
source ~/.bashrc
conda activate /bigdata/godziklab/shared/Xinru/302006/state_env

# Install datasets library if needed
pip install datasets --quiet

# Verify environment
echo "Python: $(which python)"
python -c "from datasets import load_dataset; print('datasets library OK')"

# Navigate to project
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Output directory
OUTPUT_DIR="/bigdata/godziklab/shared/Xinru/302006/01_DATA/xatlas"
mkdir -p $OUTPUT_DIR

echo ""
echo "============================================================"
echo "Downloading HCT116 (~4M cells)"
echo "============================================================"

python scripts/data/download_xatlas_orion.py \
    --cell-line HCT116 \
    --output-dir $OUTPUT_DIR \
    --format h5ad \
    --chunk-size 500000 \
    --batch-size 50000

echo ""
echo "============================================================"
echo "Downloading HEK293T (~4M cells)"
echo "============================================================"

python scripts/data/download_xatlas_orion.py \
    --cell-line HEK293T \
    --output-dir $OUTPUT_DIR \
    --format h5ad \
    --chunk-size 500000 \
    --batch-size 50000

echo ""
echo "============================================================"
echo "Download Complete"
echo "============================================================"
echo "Files saved to: $OUTPUT_DIR"
ls -lh $OUTPUT_DIR

echo "============================================================"
echo "Completed at $(date)"
echo "============================================================"
