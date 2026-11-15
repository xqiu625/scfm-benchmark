#!/bin/bash
# ===================================================
# Generate scGPT Embeddings - Batch Processing Script
#
# Processes data chunks in parallel using SLURM jobs.
# Uses batched inference for improved performance (~10x faster).
#
# Usage:
#   ./generate_scgpt_embeddings.sh /path/to/base/dir
#   ./generate_scgpt_embeddings.sh  # Uses default path
#
# Author: Project 302006
# Date: January 2026
# ===================================================

set -e  # Exit on any error

# Configuration - can be overridden by command line argument
DEFAULT_BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
BASE_DIR="${1:-$DEFAULT_BASE_DIR}"

# Validate base directory
if [ ! -d "$BASE_DIR" ]; then
    echo "ERROR: Base directory not found: $BASE_DIR"
    echo "Usage: $0 /path/to/base/dir"
    exit 1
fi

# Directory structure
CHUNKS_DIR="$BASE_DIR/data/gpu_chunks"
MODEL_DIR="$BASE_DIR/models/scGPT_human/scGPT_human"
JOBS_DIR="$BASE_DIR/jobs/scgpt_chunks"
LOGS_DIR="$BASE_DIR/logs"
CONTAINER="$BASE_DIR/07_CONTAINERS/scgpt_0.1.7.sif"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/run_scgpt_chunk.py"

# SLURM configuration
PARTITION="highmem"
MEMORY="400GB"
TIME="24:00:00"  # Reduced from 48h due to batched processing
CPUS=8
BATCH_SIZE=64  # Cells per batch for inference

echo "======================================================="
echo "scGPT EMBEDDING GENERATION - BATCH PROCESSING"
echo "======================================================="
echo "Base directory: $BASE_DIR"
echo "Chunks directory: $CHUNKS_DIR"
echo "Model directory: $MODEL_DIR"
echo "Container: $CONTAINER"
echo "Python script: $PYTHON_SCRIPT"
echo "Batch size: $BATCH_SIZE"
echo "======================================================="

# Validate required files
if [ ! -d "$MODEL_DIR" ]; then
    echo "ERROR: Model directory not found: $MODEL_DIR"
    exit 1
fi

if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Singularity container not found: $CONTAINER"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Create directories
mkdir -p "$JOBS_DIR"
mkdir -p "$LOGS_DIR"

# Count existing processed chunks
count_processed() {
    local prefix=$1
    local total=$2
    local count=0
    for i in $(seq 1 $total); do
        local output_file="${prefix}_chunk_$(printf '%03d' $i)_of_$(printf '%03d' $total)_scgpt.h5ad"
        if [ -f "$CHUNKS_DIR/$output_file" ]; then
            ((count++))
        fi
    done
    echo $count
}

covid_processed=$(count_processed "covid_platelets" 9)
hek293t_processed=$(count_processed "hek293t_platelets" 16)
hct116_processed=$(count_processed "hct116_platelets" 16)

echo ""
echo "Current processing status:"
echo "  COVID platelets:   $covid_processed/9"
echo "  HEK293T platelets: $hek293t_processed/16"
echo "  HCT116 platelets:  $hct116_processed/16"
echo "  Total: $((covid_processed + hek293t_processed + hct116_processed))/41"
echo ""

# Track submitted job IDs
job_ids=()

# Function to submit job for a single chunk
submit_chunk_job() {
    local dataset_prefix=$1
    local chunk_num=$2
    local total_chunks=$3
    local dataset_name=$4

    local chunk_file="${dataset_prefix}_chunk_$(printf '%03d' $chunk_num)_of_$(printf '%03d' $total_chunks).h5ad"
    local output_file="${dataset_prefix}_chunk_$(printf '%03d' $chunk_num)_of_$(printf '%03d' $total_chunks)_scgpt.h5ad"
    local input_path="$CHUNKS_DIR/$chunk_file"
    local output_path="$CHUNKS_DIR/$output_file"

    # Check if input exists
    if [ ! -f "$input_path" ]; then
        echo "  SKIP: Input not found: $chunk_file"
        return 1
    fi

    # Check if already processed
    if [ -f "$output_path" ]; then
        echo "  DONE: $output_file"
        return 0
    fi

    # Create job script
    local job_script="$JOBS_DIR/${dataset_prefix}_${chunk_num}.sh"

    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=scgpt_${dataset_prefix}_${chunk_num}
#SBATCH --partition=$PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEMORY
#SBATCH --time=$TIME
#SBATCH --output=$LOGS_DIR/scgpt_${dataset_prefix}_${chunk_num}_%j.out
#SBATCH --error=$LOGS_DIR/scgpt_${dataset_prefix}_${chunk_num}_%j.err

set -e

echo "======================================================="
echo "scGPT Processing: $dataset_name chunk $chunk_num/$total_chunks"
echo "======================================================="
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$SLURM_NODELIST"
echo "Start: \$(date)"
echo ""

# Input/output paths
INPUT_FILE="$input_path"
OUTPUT_FILE="$output_path"
MODEL_DIR="$MODEL_DIR"
PYTHON_SCRIPT="$PYTHON_SCRIPT"

echo "Input:  \$INPUT_FILE"
echo "Output: \$OUTPUT_FILE"
echo "Model:  \$MODEL_DIR"
echo ""

# Check if output already exists
if [ -f "\$OUTPUT_FILE" ]; then
    echo "Output already exists. Skipping."
    exit 0
fi

# Load environment
module load singularity 2>/dev/null || echo "Using system singularity"

# Run inference with batched processing
echo "Running scGPT inference with batch size $BATCH_SIZE..."
singularity exec --cleanenv \\
    --env PYTHONNOUSERSITE=1 \\
    --bind "$BASE_DIR:$BASE_DIR" \\
    "$CONTAINER" \\
    python "\$PYTHON_SCRIPT" \\
        --input "\$INPUT_FILE" \\
        --output "\$OUTPUT_FILE" \\
        --model-dir "\$MODEL_DIR" \\
        --batch-size $BATCH_SIZE

echo ""
echo "Completed: \$(date)"
echo "======================================================="
EOF

    chmod +x "$job_script"

    # Submit job
    local job_id=$(sbatch --parsable "$job_script")
    job_ids+=("$job_id")
    echo "  SUBMIT: $dataset_name chunk $chunk_num -> Job $job_id"

    # Brief pause to avoid overwhelming scheduler
    sleep 0.2
}

# Process all datasets
echo "Submitting jobs for unprocessed chunks..."
echo ""

echo "COVID Platelets (9 chunks):"
for i in $(seq 1 9); do
    submit_chunk_job "covid_platelets" $i 9 "COVID"
done

echo ""
echo "HEK293T Platelets (16 chunks):"
for i in $(seq 1 16); do
    submit_chunk_job "hek293t_platelets" $i 16 "HEK293T"
done

echo ""
echo "HCT116 Platelets (16 chunks):"
for i in $(seq 1 16); do
    submit_chunk_job "hct116_platelets" $i 16 "HCT116"
done

# Summary
echo ""
echo "======================================================="
echo "JOB SUBMISSION SUMMARY"
echo "======================================================="
echo "Total jobs submitted: ${#job_ids[@]}"
if [ ${#job_ids[@]} -gt 0 ]; then
    echo "Job IDs: ${job_ids[*]}"
fi
echo ""
echo "Configuration:"
echo "  Memory: $MEMORY"
echo "  Time:   $TIME"
echo "  CPUs:   $CPUS"
echo "  Batch:  $BATCH_SIZE cells/batch"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  watch -n 60 'ls -1 $CHUNKS_DIR/*_scgpt.h5ad 2>/dev/null | wc -l'"
echo ""
echo "Expected speedup: ~10x faster with batched processing"
echo "======================================================="
