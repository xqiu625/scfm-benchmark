#!/bin/bash
#SBATCH --job-name=umap_disease_batch
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/umap_disease_batch_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/umap_disease_batch_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# UMAP by Disease and Batch Analysis
# ============================================================================
# Investigates whether scGPT crescents correspond to batches
# ============================================================================

echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "=================================================="

# Setup environment
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Use 302006_env_clean
export ENV_PATH="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

# Verify environment
echo "Python: $(which python)"
python -c "import scanpy; print('Scanpy:', scanpy.__version__)"
python -c "import umap; print('UMAP:', umap.__version__)"

# Run analysis
echo ""
echo "Generating UMAP visualizations by disease and batch..."
echo ""

python scripts/12_umap_visualization/code/umap_by_disease_batch.py

echo ""
echo "=================================================="
echo "End Time: $(date)"
echo "=================================================="

# Copy to main figures folder
cp -v scripts/12_umap_visualization/figures/umap_5models_by_disease.png figures/ 2>/dev/null || true
cp -v scripts/12_umap_visualization/figures/umap_5models_by_batch.png figures/ 2>/dev/null || true
cp -v scripts/12_umap_visualization/figures/scgpt_detailed_analysis.png figures/ 2>/dev/null || true

echo "Job complete!"
