#!/bin/bash
# ===================================================
# OPTIMIZED scGPT Batch Submission - Better Time Estimates
# Based on observed performance: 130 cells = 90 minutes
# ===================================================

echo "🚀 OPTIMIZED SCGPT BATCH SUBMISSION"
echo "==================================="

# Base directories
BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
DATA_DIR="$BASE_DIR/data/raw"
OUTPUT_DIR="$BASE_DIR/embeddings/scgpt_full"
MODEL_DIR="$BASE_DIR/models/scGPT_human/scGPT_human"
SCRIPT_PATH="$BASE_DIR/302006_updategit/scripts/scGPT/generate_scgpt_full_scale.py"
SINGULARITY_IMAGE="$BASE_DIR/scgpt_0.1.7.sif"

# Create output directory
mkdir -p $OUTPUT_DIR

# EXTENDED Dataset configurations with generous time buffers
# Performance: 130 cells = 90 minutes = ~0.7 min/cell
# Adding 50-100% time buffer for safety

datasets=(
    # Large epithelial datasets - EXTENDED time
    "hek293t_epithelial:hek293t_epithelial.h5ad:36:500"     # ~47K cells = ~330 min = 5.5h → 36h (6x buffer)
    "hct116_epithelial:hct116_epithelial.h5ad:32:500"      # ~45K cells = ~315 min = 5.25h → 32h (6x buffer)
    
    # Large COVID dataset - EXTENDED time
    "human_platelet_covid_severity:human_platelet_covid_severity.h5ad:28:450"  # ~48K cells → 28h
    
    # Smaller platelet datasets - EXTENDED time
    "hek293t_platelets:hek293t_platelets.h5ad:12:300"      # ~12K cells = ~84 min = 1.4h → 12h (8x buffer)
    "hct116_platelets:hct116_platelets.h5ad:12:300"        # ~11K cells = ~77 min = 1.3h → 12h (8x buffer)
)

echo "📋 EXTENDED Dataset estimates (very generous time buffers):"
for dataset in "${datasets[@]}"; do
    IFS=':' read -r name file time mem <<< "$dataset"
    echo "  • $name ($file) - Time: ${time}h, Mem: ${mem}GB"
done

echo "🐋 Using Singularity container: $SINGULARITY_IMAGE"

# Check current job status
echo -e "\n📊 Current job status:"
squeue -u $USER

echo -e "\n🔧 Creating optimized SLURM job scripts..."

# Generate individual job scripts for each dataset
for dataset in "${datasets[@]}"; do
    IFS=':' read -r name input_file est_hours mem_gb <<< "$dataset"
    
    job_script="${BASE_DIR}/jobs/scgpt_${name}_optimized.sh"
    mkdir -p "${BASE_DIR}/jobs"
    
    cat > "$job_script" << JOBEOF
#!/bin/bash
#SBATCH --job-name=scgpt_${name}_opt
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=${mem_gb}GB
#SBATCH --time=${est_hours}:00:00
#SBATCH --output=${BASE_DIR}/logs/scgpt_${name}_optimized_%j.out
#SBATCH --error=${BASE_DIR}/logs/scgpt_${name}_optimized_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xinru.qiu@ucr.edu

echo "🚀 Starting OPTIMIZED scGPT embedding generation for ${name}"
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$SLURM_NODELIST"
echo "Start time: \$(date)"
echo "Input: ${input_file}"
echo "Allocated: ${mem_gb}GB RAM, ${est_hours}h time"
echo "Container: ${SINGULARITY_IMAGE}"

# Load modules if needed
module load singularity 2>/dev/null || echo "Singularity module not needed"

# Change to working directory
cd ${BASE_DIR}

# Check if input file exists
if [ ! -f "${DATA_DIR}/${input_file}" ]; then
    echo "❌ Input file not found: ${DATA_DIR}/${input_file}"
    exit 1
fi

# Check if output already exists
if [ -f "${OUTPUT_DIR}/${name}_scgpt_full.h5ad" ]; then
    echo "⏭️  Output already exists: ${OUTPUT_DIR}/${name}_scgpt_full.h5ad"
    echo "File size: \$(du -h ${OUTPUT_DIR}/${name}_scgpt_full.h5ad | cut -f1)"
    echo "Skipping processing."
    exit 0
fi

# Pre-check dataset size
echo "📊 Checking dataset size..."
singularity exec ${SINGULARITY_IMAGE} python << 'SIZE_CHECK_EOF'
import scanpy as sc
from pathlib import Path

input_path = Path('${DATA_DIR}/${input_file}')
try:
    adata = sc.read_h5ad(input_path)
    print(f"📊 Dataset: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Estimate processing time
    cells = adata.shape[0]
    estimated_minutes = cells * 0.7  # Based on observed 0.7 min/cell
    estimated_hours = estimated_minutes / 60
    
    print(f"⏱️  Estimated processing time: {estimated_minutes:.1f} minutes ({estimated_hours:.1f} hours)")
    print(f"🕐 Time allocated: ${est_hours} hours")
    
    if estimated_hours > ${est_hours} * 0.9:
        print(f"⚠️  WARNING: Estimated time ({estimated_hours:.1f}h) close to limit (${est_hours}h)")
    else:
        print(f"✅ Time allocation looks sufficient")
        
except Exception as e:
    print(f"❌ Error checking dataset: {e}")
SIZE_CHECK_EOF

# Run scGPT embedding generation using Singularity
echo "🐋 Running scGPT with Singularity container..."

singularity exec --nv ${SINGULARITY_IMAGE} python << 'PYTHON_EOF'
import sys
import os
sys.path.append('${BASE_DIR}/302006_updategit/scripts/scGPT')

# Import required modules
import scgpt as scg
import scanpy as sc
import numpy as np
import torch
import json
import time
import gc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

print(f"🔧 Python version: {sys.version}")
print(f"🔧 PyTorch version: {torch.__version__}")
print(f"🔧 scGPT version: {scg.__version__}")

# Load the functions from your main script
exec(open('${SCRIPT_PATH}').read())

# Process this specific dataset
input_path = Path('${DATA_DIR}/${input_file}')
output_path = Path('${OUTPUT_DIR}/${name}_scgpt_full.h5ad')
model_dir = Path('${MODEL_DIR}')

print(f"📂 Processing: {input_path}")
print(f"💾 Output: {output_path}")
print(f"🤖 Model: {model_dir}")

# Check if input exists
if not input_path.exists():
    print(f"❌ Input file not found: {input_path}")
    sys.exit(1)

try:
    success = process_dataset_full_scale(input_path, output_path, model_dir, device='cpu')
    if success:
        print(f"✅ SUCCESS: {output_path}")
    else:
        print(f"❌ FAILED: {input_path}")
        sys.exit(1)
except Exception as e:
    print(f"❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

PYTHON_EOF

echo "Completed at: \$(date)"
echo "🎉 scGPT embedding generation finished for ${name}"
JOBEOF

    chmod +x "$job_script"
    echo "  ✅ Created: $job_script"
done

# Create logs directory
mkdir -p "${BASE_DIR}/logs"

# Cancel any running jobs with the old names
echo -e "\n🛑 Canceling any existing scGPT jobs..."
scancel --name=scgpt_hek293t_epithelial 2>/dev/null || true
scancel --name=scgpt_hct116_epithelial 2>/dev/null || true
scancel --name=scgpt_hek293t_platelets 2>/dev/null || true
scancel --name=scgpt_hct116_platelets 2>/dev/null || true
scancel --name=scgpt_human_platelet_covid_severity 2>/dev/null || true

echo -e "\n📤 Submitting OPTIMIZED jobs to SLURM..."

# Submit jobs - start with smaller datasets first
job_ids=()

# Submit in order of increasing time requirements
submit_order=(
    "hct116_platelets:hct116_platelets.h5ad:12:300"        # 12h
    "hek293t_platelets:hek293t_platelets.h5ad:12:300"      # 12h  
    "human_platelet_covid_severity:human_platelet_covid_severity.h5ad:28:450"  # 28h
    "hct116_epithelial:hct116_epithelial.h5ad:32:500"      # 32h
    "hek293t_epithelial:hek293t_epithelial.h5ad:36:500"    # 36h (longest)
)

for dataset in "${submit_order[@]}"; do
    IFS=':' read -r name input_file est_hours mem_gb <<< "$dataset"
    
    job_script="${BASE_DIR}/jobs/scgpt_${name}_optimized.sh"
    
    # Submit job
    job_id=$(sbatch --parsable "$job_script")
    job_ids+=($job_id)
    echo "  📤 Submitted ${name}: Job ID $job_id (${est_hours}h, ${mem_gb}GB)"
    
    # Small delay to avoid overwhelming scheduler
    sleep 2
done

echo -e "\n📊 Job Summary:"
echo "Job IDs: ${job_ids[*]}"
echo "Total jobs: ${#job_ids[@]}"

echo -e "\n🔍 Monitor jobs with:"
echo "squeue -u \$USER"
echo "watch -n 30 'squeue -u \$USER'"

echo -e "\n📁 Check logs at:"
echo "${BASE_DIR}/logs/scgpt_*_optimized_*.out"

echo -e "\n⏱️  REVISED completion time estimates:"
total_hours=0
for dataset in "${datasets[@]}"; do
    IFS=':' read -r name input_file est_hours mem_gb <<< "$dataset"
    echo "  • $name: ~${est_hours} hours"
    total_hours=$((total_hours + est_hours))
done

echo -e "\n🎯 All OPTIMIZED scGPT jobs submitted!"
echo "Total estimated time: ${total_hours} hours (running in parallel)"
echo "Expected completion: within 20-24 hours"