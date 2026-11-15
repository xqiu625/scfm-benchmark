#!/bin/bash
#SBATCH --job-name=12_umap_viz
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/12_umap_visualization_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/12_umap_visualization_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC-12: UMAP Visualization for All 5 Foundation Models
# ============================================================================
# Models: STATE, UCE, scGPT (cleaned), TranscriptFormer, Geneformer
# Outputs: UMAP visualizations (6-class severity and binary)
# scGPT now uses cleaned embeddings (Orion data removed)
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

# Verify cleaned scGPT file exists
echo ""
echo "Checking cleaned scGPT file..."
ls -lh /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/scgpt/covid_scgpt_clean.h5ad

# Base directories
BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
GIT_DIR="/bigdata/godziklab/shared/Xinru/302006/302006_updategit"
RESULTS_DIR="$BASE_DIR/05_RESULTS/covid_severity_benchmark/figures"

# Run main UMAP visualization (6-class and binary)
echo ""
echo "=========================================="
echo "Generating UMAP visualizations (6-class, binary)"
echo "Models: STATE, UCE, scGPT (cleaned), TranscriptFormer, Geneformer"
echo "=========================================="

python scripts/12_umap_visualization/code/create_embedding_visualizations.py

# Run disease/batch UMAP visualization
echo ""
echo "=========================================="
echo "Generating UMAP by disease and batch"
echo "=========================================="

python scripts/12_umap_visualization/code/umap_by_disease_batch.py

echo ""
echo "=================================================="
echo "End Time: $(date)"
echo "=================================================="

# Copy figures to git repo
echo ""
echo "Copying figures..."

# From results dir
cp -v "$RESULTS_DIR/umap_5models"*.png "$GIT_DIR/figures/" 2>/dev/null || true
cp -v "$RESULTS_DIR/umap_5models"*.pdf "$GIT_DIR/figures/" 2>/dev/null || true
cp -v "$RESULTS_DIR/"*.png "$GIT_DIR/scripts/12_umap_visualization/figures/" 2>/dev/null || true

# From analyses dir (disease/batch UMAPs)
cp -v "$GIT_DIR/scripts/12_umap_visualization/figures/umap_5models"*.png "$GIT_DIR/figures/" 2>/dev/null || true

echo ""
echo "Figures saved to:"
echo "  - $RESULTS_DIR"
echo "  - $GIT_DIR/scripts/12_umap_visualization/figures/"
echo "  - $GIT_DIR/figures/"

echo ""
echo "Job complete!"
