#!/bin/bash
# FINAL FIX: Remove Entrez IDs (plain numbers) and keep only valid Ensembl IDs (ENSG*)
# Issue: mygene.info returned mix of Ensembl IDs and Entrez IDs
# Solution: Filter to keep only IDs starting with "ENSG"

echo "=========================================="
echo "FINAL FIX: REMOVE ENTREZ IDS"
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

# Re-prepare HEK293T and HCT116 - keep only valid Ensembl IDs
echo ""
echo "Step 3: Filtering for valid Ensembl IDs (ENSG* only)..."
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
        logging.FileHandler(f'fix_entrez_ids_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def fix_dataset(dataset_name, input_file):
    """Keep only valid Ensembl IDs (starting with ENSG)"""
    logger.info(f"\n{'='*80}")
    logger.info(f"Fixing: {dataset_name}")
    logger.info(f"{'='*80}")

    # Load data
    logger.info(f"Loading {input_file.name}...")
    adata = sc.read_h5ad(input_file)
    logger.info(f"  Shape before: {adata.shape}")

    # Filter for valid Ensembl IDs only
    logger.info("\nFiltering for valid Ensembl IDs (ENSG* only)...")
    ensembl_ids = adata.var['ensembl_id'].values

    def is_valid_ensembl(eid):
        if pd.isna(eid):
            return False
        eid_str = str(eid)
        return eid_str.startswith('ENSG')

    valid_mask = np.array([is_valid_ensembl(eid) for eid in ensembl_ids])

    n_valid = valid_mask.sum()
    n_invalid = (~valid_mask).sum()

    logger.info(f"  Valid Ensembl IDs (ENSG*): {n_valid}")
    logger.info(f"  Invalid IDs (Entrez/other): {n_invalid}")

    # Show examples of invalid IDs
    if n_invalid > 0:
        invalid_ids = ensembl_ids[~valid_mask]
        logger.info(f"  Example invalid IDs: {list(invalid_ids[:5])}")

    # Keep only valid genes
    adata = adata[:, valid_mask].copy()
    logger.info(f"  Shape after filtering: {adata.shape}")

    # Now strip version numbers and check for duplicates
    logger.info("\nStripping version numbers and checking duplicates...")
    ensembl_ids_clean = adata.var['ensembl_id'].values
    ensembl_ids_no_version = np.array([str(eid).split('.')[0] for eid in ensembl_ids_clean])

    # Check for duplicates
    unique_ids, unique_indices = np.unique(ensembl_ids_no_version, return_index=True)

    n_duplicates = len(ensembl_ids_no_version) - len(unique_ids)
    if n_duplicates > 0:
        logger.warning(f"  Found {n_duplicates} duplicate Ensembl IDs")
        logger.warning(f"  Keeping first occurrence")

        adata = adata[:, unique_indices].copy()
        adata.var['ensembl_id'] = ensembl_ids_no_version[unique_indices]
        logger.info(f"  Shape after removing duplicates: {adata.shape}")
    else:
        logger.info(f"  ✓ No duplicates found")
        adata.var['ensembl_id'] = ensembl_ids_no_version

    # Save fixed data
    output_file = Path(f"/bigdata/godziklab/shared/Xinru/302006/01_DATA/transcriptformer_ready/{dataset_name}_transcriptformer_ready.h5ad")
    logger.info(f"\nSaving fixed data to: {output_file}")
    adata.write_h5ad(output_file)

    logger.info(f"✓ {dataset_name} fixed successfully")
    logger.info(f"  Final shape: {adata.shape}")
    logger.info(f"  All IDs start with ENSG: {all(str(eid).startswith('ENSG') for eid in adata.var['ensembl_id'])}")

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
            import traceback
            traceback.print_exc()
            continue
    else:
        logger.error(f"File not found: {input_file}")

logger.info("\n" + "="*80)
logger.info("ENTREZ ID REMOVAL COMPLETE")
logger.info("="*80)
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Data fixing failed"
    exit 1
fi

# Verify the fix
echo ""
echo "Verifying fix..."
python scripts/foundation_models/TranscriptFormer/diagnose_duplicates.py

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
python - <<'EOF'
import scanpy as sc
from pathlib import Path

base_path = Path('/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer')
datasets = ['covid', 'hek293t', 'hct116']

print("\nFinal Embedding Summary:")
print("="*60)
for dataset in datasets:
    file_path = base_path / f"{dataset}_transcriptformer.h5ad"
    if file_path.exists():
        adata = sc.read_h5ad(file_path)
        print(f"{dataset.upper():8s}: {adata.shape[0]:6,} cells × {adata.obsm['X_transcriptformer'].shape[1]:4,} dims ✓")
    else:
        print(f"{dataset.upper():8s}: File not found")
print("="*60)
EOF

echo ""
echo "Next steps:"
echo "  1. Merge with UCE and scGPT embeddings"
echo "  2. Run three-model comparison analysis"
echo "=========================================="
