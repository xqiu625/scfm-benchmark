#!/bin/bash
#SBATCH --job-name=scgpt_bp
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scgpt_best_practices_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scgpt_best_practices_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# scGPT Embedding Generation with Best Practices
# ============================================================================
# Following official scGPT documentation:
# - Proper preprocessing with normalization and log1p
# - HVG selection (1200 genes)
# - Batch-aware processing (batch_key=Data_NO)
# - Value binning (51 bins)
# ============================================================================

echo "============================================================"
echo "scGPT EMBEDDING GENERATION (Best Practices)"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo ""

# Paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
GIT_DIR="$BASE_DIR/302006_updategit"
INPUT_FILE="$BASE_DIR/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
OUTPUT_DIR="$BASE_DIR/02_EMBEDDINGS/scgpt_best_practices"
OUTPUT_FILE="$OUTPUT_DIR/covid_scgpt_best_practices.h5ad"
MODEL_DIR="$BASE_DIR/models/scGPT_human/scGPT_human"
CONTAINER="$BASE_DIR/07_CONTAINERS/scgpt_0.1.7.sif"

# Setup environment
cd "$GIT_DIR"

# Load singularity module
module load singularity 2>/dev/null || true

# Verify container exists
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Singularity container not found: $CONTAINER"
    exit 1
fi
echo "Container: $CONTAINER"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Configuration:"
echo "  Input: $INPUT_FILE"
echo "  Output: $OUTPUT_FILE"
echo "  Model: $MODEL_DIR"
echo "  HVG: 1200 genes"
echo "  Bins: 51"
echo "  Batch key: Data_NO"
echo ""

# Check input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# Check model files exist
if [ ! -d "$MODEL_DIR" ]; then
    echo "ERROR: Model directory not found: $MODEL_DIR"
    exit 1
fi

echo "============================================================"
echo "Running scGPT Best Practices Pipeline"
echo "============================================================"
echo ""

# Run with Singularity container
# Use --cleanenv to prevent host environment from interfering
# Set PYTHONNOUSERSITE=1 to prevent loading packages from ~/.local
echo "Using Singularity container: $CONTAINER"
singularity exec --nv --cleanenv \
    --bind "$BASE_DIR:$BASE_DIR" \
    --bind "$GIT_DIR:$GIT_DIR" \
    --env PYTHONNOUSERSITE=1 \
    "$CONTAINER" \
    python "$GIT_DIR/analyses/16_embedding_generation_5models/code/run_scgpt_best_practices.py" \
        --input "$INPUT_FILE" \
        --output "$OUTPUT_FILE" \
        --model-dir "$MODEL_DIR" \
        --n-hvg 1200 \
        --n-bins 51 \
        --batch-size 64 \
        --batch-key "Data_NO"

# Check if output was created
if [ -f "$OUTPUT_FILE" ]; then
    echo ""
    echo "============================================================"
    echo "SUCCESS: Output file created"
    echo "============================================================"
    echo "Output: $OUTPUT_FILE"
    ls -lh "$OUTPUT_FILE"

    # Copy to analyses directory
    cp -v "$OUTPUT_FILE" "$GIT_DIR/analyses/16_embedding_generation_5models/scGPT_results/"
else
    echo ""
    echo "============================================================"
    echo "ERROR: Output file not created"
    echo "============================================================"
    exit 1
fi

echo ""
echo "============================================================"
echo "End Time: $(date)"
echo "============================================================"
