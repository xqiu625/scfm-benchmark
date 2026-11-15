#!/bin/bash
#SBATCH --job-name=esm_multimodal
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=esm_multimodal_%j.out
#SBATCH --error=esm_multimodal_%j.err

# ============================================================================
# ACC 22: ESM Protein Multimodal Fusion Pipeline
#
# Runs steps 1-5:
#   1. Identify DEGs (Wilcoxon rank-sum)
#   2. Retrieve UniProt protein sequences
#   3. Generate ESM-2 protein embeddings (GPU)
#   4. Build per-cell protein embeddings (expression-weighted)
#   5. 5-modality attention fusion
#
# Usage:
#   sbatch submit_22_esm_multimodal.sh
#   sbatch submit_22_esm_multimodal.sh --start-step 3   # resume from step 3
# ============================================================================

set -euo pipefail

# Parse optional start step
START_STEP=${1:-1}

# Environment setup
export PYTHONNOUSERSITE=1
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean

# Paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results"
mkdir -p "${RESULTS_DIR}"

echo "============================================"
echo "ACC 22: ESM Protein Multimodal Fusion"
echo "Start step: ${START_STEP}"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'none')"
echo "============================================"

# Step 1: Identify DEGs
if [ "${START_STEP}" -le 1 ]; then
    echo ""
    echo ">>> Step 1: Identify DEGs"
    python "${SCRIPT_DIR}/step1_identify_degs.py" \
        --base-dir "${BASE_DIR}" \
        --n-top 500 \
        --output-dir "${RESULTS_DIR}"
    echo "Step 1 complete: $(date)"
fi

# Step 2: Get protein sequences from UniProt
if [ "${START_STEP}" -le 2 ]; then
    echo ""
    echo ">>> Step 2: Get protein sequences"
    python "${SCRIPT_DIR}/step2_get_protein_sequences.py" \
        --input "${RESULTS_DIR}/deg_gene_list.txt" \
        --output-dir "${RESULTS_DIR}"
    echo "Step 2 complete: $(date)"
fi

# Step 3: Generate ESM-2 embeddings (GPU-intensive)
if [ "${START_STEP}" -le 3 ]; then
    echo ""
    echo ">>> Step 3: Generate ESM-2 embeddings"
    python "${SCRIPT_DIR}/step3_generate_esm_embeddings.py" \
        --model esm2_t33_650M_UR50D \
        --batch-size 4 \
        --input "${RESULTS_DIR}/deg_proteins.fasta" \
        --output-dir "${RESULTS_DIR}"
    echo "Step 3 complete: $(date)"
fi

# Step 4: Build per-cell protein embeddings
if [ "${START_STEP}" -le 4 ]; then
    echo ""
    echo ">>> Step 4: Build per-cell protein embeddings"
    python "${SCRIPT_DIR}/step4_build_cell_protein_embeddings.py" \
        --base-dir "${BASE_DIR}" \
        --output-dir "${RESULTS_DIR}"
    echo "Step 4 complete: $(date)"
fi

# Step 5: 5-modality attention fusion
if [ "${START_STEP}" -le 5 ]; then
    echo ""
    echo ">>> Step 5: 5-modality fusion"
    python "${SCRIPT_DIR}/step5_multimodal_fusion_5model.py" \
        --base-dir "${BASE_DIR}" \
        --esm-results-dir "${RESULTS_DIR}" \
        --output-dir "${RESULTS_DIR}" \
        --strategy all \
        --attention-epochs 100 \
        --common-dim 256
    echo "Step 5 complete: $(date)"
fi

echo ""
echo "============================================"
echo "Pipeline complete: $(date)"
echo "Results in: ${RESULTS_DIR}"
echo "============================================"
