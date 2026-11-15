#!/bin/bash
# Fix duplicate Ensembl IDs and regenerate embeddings for HEK293T and HCT116
# COVID dataset already completed successfully

echo "=========================================="
echo "FIX DUPLICATE GENES - CELL LINE DATASETS"
echo "=========================================="

# Activate conda environment
echo ""
echo "Step 1: Activating transcriptformer environment..."
eval "$(conda shell.bash hook)"
conda activate transcriptformer

# Verify environment
echo ""
echo "Step 2: Verifying environment..."
python --version

# Re-prepare only HEK293T and HCT116 datasets
echo ""
echo "Step 3: Re-preparing HEK293T and HCT116 with duplicate removal..."
python - <<'EOF'
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from datetime import datetime
import sys

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'fix_cell_lines_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def fix_dataset(dataset_name, input_file):
    """Fix duplicate Ensembl IDs in a single dataset"""
    logger.info(f"\n{'='*80}")
    logger.info(f"Fixing: {dataset_name}")
    logger.info(f"{'='*80}")

    # Load data
    logger.info(f"Loading {input_file.name}...")
    adata = sc.read_h5ad(input_file)
    logger.info(f"  Shape before: {adata.shape}")

    # Check for duplicate Ensembl IDs
    logger.info("\nChecking for duplicate Ensembl IDs...")
    ensembl_ids = adata.var['ensembl_id'].values
    unique_ids, unique_indices = np.unique(ensembl_ids, return_index=True)

    n_duplicates = len(ensembl_ids) - len(unique_ids)
    if n_duplicates > 0:
        logger.warning(f"  Found {n_duplicates} duplicate Ensembl IDs")
        logger.warning(f"  Keeping first occurrence of each duplicate")

        # Keep only first occurrence of each Ensembl ID
        adata = adata[:, unique_indices].copy()
        logger.info(f"  Shape after removing duplicates: {adata.shape}")
    else:
        logger.info(f"  No duplicate Ensembl IDs found")

    # Save fixed data
    output_file = Path(f"/bigdata/godziklab/shared/Xinru/302006/01_DATA/transcriptformer_ready/{dataset_name}_transcriptformer_ready.h5ad")
    logger.info(f"\nSaving fixed data to: {output_file}")
    adata.write_h5ad(output_file)

    logger.info(f"✓ {dataset_name} fixed successfully")
    logger.info(f"  Final shape: {adata.shape}")

    return adata

# Process datasets
base_dir = Path("/bigdata/godziklab/shared/Xinru/302006/01_DATA/transcriptformer_ready")

datasets = {
    'hek293t': 'hek293t_transcriptformer_ready.h5ad',
    'hct116': 'hct116_transcriptformer_ready.h5ad'
}

for dataset_name, filename in datasets.items():
    input_file = base_dir / filename
    if input_file.exists():
        try:
            adata = fix_dataset(dataset_name, input_file)
        except Exception as e:
            logger.error(f"Failed to fix {dataset_name}: {str(e)}")
            continue
    else:
        logger.error(f"File not found: {input_file}")

logger.info("\n" + "="*80)
logger.info("DUPLICATE FIXING COMPLETE")
logger.info("="*80)
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Data fixing failed"
    exit 1
fi

# Generate embeddings for HEK293T
echo ""
echo "Step 4: Generating TranscriptFormer embeddings for HEK293T..."
python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --dataset hek293t --batch-size 16

if [ $? -ne 0 ]; then
    echo "ERROR: HEK293T embedding generation failed"
    exit 1
fi

# Generate embeddings for HCT116
echo ""
echo "Step 5: Generating TranscriptFormer embeddings for HCT116..."
python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --dataset hct116 --batch-size 16

if [ $? -ne 0 ]; then
    echo "ERROR: HCT116 embedding generation failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "ALL DATASETS COMPLETE!"
echo "=========================================="
echo ""
echo "Generated embeddings:"
echo "  COVID:   /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/covid_transcriptformer.h5ad"
echo "  HEK293T: /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/hek293t_transcriptformer.h5ad"
echo "  HCT116:  /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/hct116_transcriptformer.h5ad"
echo ""
echo "Next steps:"
echo "  1. Verify all embedding shapes and quality"
echo "  2. Merge with UCE and scGPT embeddings"
echo "  3. Run three-model comparison analysis"
echo "=========================================="
