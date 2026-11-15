#!/bin/bash
# ===================================================
# Submit UCE Embedding Generation to SLURM
#
# This is a thin wrapper that submits the Python script to SLURM.
# The scientific logic is in generate_uce_embeddings.py.
#
# Usage:
#   ./submit_uce_embeddings.sh [dataset]
#   ./submit_uce_embeddings.sh covid
#   ./submit_uce_embeddings.sh all
#
# Author: Project 302006
# Date: January 2026
# ===================================================

set -e

# Configuration
DEFAULT_BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
BASE_DIR="${BASE_DIR:-$DEFAULT_BASE_DIR}"
DATASET="${1:-covid}"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/generate_uce_embeddings.py"

# Directories
LOGS_DIR="$BASE_DIR/logs"
JOBS_DIR="$BASE_DIR/jobs/uce"
CONDA_ENV="$BASE_DIR/../302006_env"

# SLURM configuration
PARTITION="gpu"
QOS="gpu"
MEMORY="250GB"
TIME="24:00:00"
GPUS=1

# Validate
if [ ! -d "$BASE_DIR" ]; then
    echo "ERROR: Base directory not found: $BASE_DIR"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Create directories
mkdir -p "$LOGS_DIR"
mkdir -p "$JOBS_DIR"

echo "======================================================="
echo "UCE EMBEDDING SUBMISSION"
echo "======================================================="
echo "Base directory: $BASE_DIR"
echo "Dataset:        $DATASET"
echo "Python script:  $PYTHON_SCRIPT"
echo "Conda env:      $CONDA_ENV"
echo "======================================================="

# Create job script
JOB_SCRIPT="$JOBS_DIR/uce_${DATASET}.sh"

cat > "$JOB_SCRIPT" << EOF
#!/bin/bash
#SBATCH --job-name=uce_${DATASET}
#SBATCH --partition=$PARTITION
#SBATCH --qos=$QOS
#SBATCH --time=$TIME
#SBATCH --mem=$MEMORY
#SBATCH --gres=gpu:$GPUS
#SBATCH --output=$LOGS_DIR/uce_${DATASET}_%j.out
#SBATCH --error=$LOGS_DIR/uce_${DATASET}_%j.err

set -e

echo "======================================================="
echo "UCE Embedding Generation: $DATASET"
echo "======================================================="
echo "Job ID: \$SLURM_JOB_ID"
echo "Node:   \$SLURM_NODELIST"
echo "GPU:    \$CUDA_VISIBLE_DEVICES"
echo "Start:  \$(date)"
echo ""

# Activate conda environment
source activate "$CONDA_ENV" || {
    echo "ERROR: Failed to activate conda environment"
    exit 1
}

echo "Python: \$(which python)"
echo "Conda:  \$CONDA_DEFAULT_ENV"
echo ""

# Run embedding generation
python "$PYTHON_SCRIPT" \\
    --base-dir "$BASE_DIR" \\
    --dataset "$DATASET" \\
    --batch-size 200

echo ""
echo "Completed: \$(date)"
echo "======================================================="
EOF

chmod +x "$JOB_SCRIPT"

# Submit job
JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")
echo ""
echo "Submitted job: $JOB_ID"
echo "Monitor: squeue -j $JOB_ID"
echo "Logs:    tail -f $LOGS_DIR/uce_${DATASET}_${JOB_ID}.out"
