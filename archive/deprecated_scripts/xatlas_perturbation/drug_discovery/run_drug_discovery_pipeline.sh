#!/bin/bash
# =============================================================================
# Drug Discovery Pipeline for COVID-19/Sepsis Therapeutic Targets
# =============================================================================
# This pipeline:
# 1. Retrieves protein sequences for top therapeutic candidates
# 2. Downloads AlphaFold structures
# 3. Generates ESM-2 protein embeddings
# 4. Searches drug databases for existing drugs
# 5. Prepares structures for molecular docking
#
# Prerequisites:
#   pip install requests pandas biopython fair-esm torch
#   Optional: arc-state, AutoDock Vina, DiffDock
#
# Usage:
#   bash run_drug_discovery_pipeline.sh [--candidates CSV_PATH] [--output-dir DIR]
# =============================================================================

set -e  # Exit on error

# Default parameters
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
OUTPUT_DIR="${PROJECT_DIR}/results/drug_discovery"
CANDIDATES_CSV=""
TOP_N=10

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --candidates)
            CANDIDATES_CSV="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --top-n)
            TOP_N="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "DRUG DISCOVERY PIPELINE"
echo "============================================================"
echo "Output directory: $OUTPUT_DIR"
echo "Top N candidates: $TOP_N"
if [ -n "$CANDIDATES_CSV" ]; then
    echo "Candidates file: $CANDIDATES_CSV"
fi
echo ""

# -----------------------------------------------------------------------------
# Step 1: Retrieve Protein Sequences
# -----------------------------------------------------------------------------
echo "============================================================"
echo "STEP 1: Retrieving protein sequences from UniProt"
echo "============================================================"

PROTEIN_DIR="$OUTPUT_DIR/protein_data"
mkdir -p "$PROTEIN_DIR"

if [ -n "$CANDIDATES_CSV" ]; then
    python "$SCRIPT_DIR/protein_structure/get_protein_sequences.py" \
        --candidates-csv "$CANDIDATES_CSV" \
        --output-dir "$PROTEIN_DIR" \
        --top-n "$TOP_N"
else
    python "$SCRIPT_DIR/protein_structure/get_protein_sequences.py" \
        --output-dir "$PROTEIN_DIR" \
        --top-n "$TOP_N"
fi

# Find the summary file
SUMMARY_FILE=$(ls -t "$PROTEIN_DIR"/target_summary_*.csv | head -1)
FASTA_FILE=$(ls -t "$PROTEIN_DIR"/target_sequences_*.fasta | head -1)

echo "Summary file: $SUMMARY_FILE"
echo "FASTA file: $FASTA_FILE"

# -----------------------------------------------------------------------------
# Step 2: Download AlphaFold Structures
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "STEP 2: Downloading AlphaFold structures"
echo "============================================================"

STRUCTURE_DIR="$OUTPUT_DIR/structures"
mkdir -p "$STRUCTURE_DIR"

python "$SCRIPT_DIR/protein_structure/download_alphafold_structures.py" \
    --summary "$SUMMARY_FILE" \
    --output "$STRUCTURE_DIR"

# -----------------------------------------------------------------------------
# Step 3: Generate ESM-2 Embeddings
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "STEP 3: Generating ESM-2 protein embeddings"
echo "============================================================"

ESM_DIR="$OUTPUT_DIR/esm_embeddings"
mkdir -p "$ESM_DIR"

# Check if ESM is installed
if python -c "import esm" 2>/dev/null; then
    python "$SCRIPT_DIR/protein_structure/generate_esm_embeddings.py" \
        --fasta "$FASTA_FILE" \
        --output-dir "$ESM_DIR"
else
    echo "WARNING: ESM not installed. Skipping embedding generation."
    echo "Install with: pip install fair-esm"
fi

# -----------------------------------------------------------------------------
# Step 4: Search Drug Databases
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "STEP 4: Searching drug databases"
echo "============================================================"

DRUG_DIR="$OUTPUT_DIR/drug_search"
mkdir -p "$DRUG_DIR"

if [ -n "$CANDIDATES_CSV" ]; then
    python "$SCRIPT_DIR/drug_repurposing/search_drug_databases.py" \
        --candidates-csv "$CANDIDATES_CSV" \
        --output-dir "$DRUG_DIR" \
        --top-n "$TOP_N"
else
    python "$SCRIPT_DIR/drug_repurposing/search_drug_databases.py" \
        --output-dir "$DRUG_DIR" \
        --top-n "$TOP_N"
fi

# -----------------------------------------------------------------------------
# Step 5: Prepare for Molecular Docking
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "STEP 5: Preparing structures for molecular docking"
echo "============================================================"

DOCKING_DIR="$OUTPUT_DIR/docking"
mkdir -p "$DOCKING_DIR"

python "$SCRIPT_DIR/drug_repurposing/prepare_docking.py" \
    --structures-dir "$STRUCTURE_DIR" \
    --output-dir "$DOCKING_DIR" \
    --method both

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "============================================================"
echo ""
echo "Output directories:"
echo "  Protein sequences: $PROTEIN_DIR"
echo "  AlphaFold structures: $STRUCTURE_DIR"
echo "  ESM embeddings: $ESM_DIR"
echo "  Drug search results: $DRUG_DIR"
echo "  Docking inputs: $DOCKING_DIR"
echo ""
echo "Key output files:"
echo "  - $(ls -t "$PROTEIN_DIR"/target_summary_*.csv 2>/dev/null | head -1)"
echo "  - $(ls -t "$DRUG_DIR"/drug_summary_*.csv 2>/dev/null | head -1)"
echo "  - $(ls -t "$STRUCTURE_DIR"/structure_quality_*.csv 2>/dev/null | head -1)"
echo ""
echo "Next steps:"
echo "  1. Review drug_summary.csv for repurposing candidates"
echo "  2. Run molecular docking with prepared structures"
echo "  3. Validate top hits with literature search"
echo ""
