#!/bin/bash
# Final fix: Re-prepare HEK293T and HCT116 with version-stripped Ensembl IDs
# Then generate embeddings

echo "=========================================="
echo "FINAL FIX: VERSION-STRIPPED ENSEMBL IDS"
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

# Re-prepare HEK293T and HCT116 with version stripping
echo ""
echo "Step 3: Re-preparing HEK293T and HCT116 with version-stripped Ensembl IDs..."
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
        logging.FileHandler(f'final_fix_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def fix_dataset(dataset_name, input_file):
    """Fix Ensembl IDs by stripping version numbers and removing duplicates"""
    logger.info(f"\n{'='*80}")
    logger.info(f"Fixing: {dataset_name}")
    logger.info(f"{'='*80}")

    # Load data
    logger.info(f"Loading {input_file.name}...")
    adata = sc.read_h5ad(input_file)
    logger.info(f"  Shape before: {adata.shape}")

    # Strip version numbers from Ensembl IDs
    logger.info("\nStripping version numbers from Ensembl IDs...")
    ensembl_ids = adata.var['ensembl_id'].values
    ensembl_ids_no_version = np.array([str(eid).split('.')[0] if pd.notna(eid) else eid for eid in ensembl_ids])

    # Check for duplicates after version removal
    unique_ids, unique_indices = np.unique(ensembl_ids_no_version, return_index=True)

    n_duplicates = len(ensembl_ids_no_version) - len(unique_ids)
    if n_duplicates > 0:
        logger.warning(f"  Found {n_duplicates} duplicate Ensembl IDs (after removing version numbers)")
        logger.warning(f"  Keeping first occurrence of each duplicate")

        # Keep only first occurrence
        adata = adata[:, unique_indices].copy()
        adata.var['ensembl_id'] = ensembl_ids_no_version[unique_indices]
        logger.info(f"  Shape after removing duplicates: {adata.shape}")
    else:
        logger.info(f"  No duplicate Ensembl IDs found")
        adata.var['ensembl_id'] = ensembl_ids_no_version

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
            import traceback
            traceback.print_exc()
            continue
    else:
        logger.error(f"File not found: {input_file}")

logger.info("\n" + "="*80)
logger.info("VERSION STRIPPING AND DUPLICATE REMOVAL COMPLETE")
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
echo "Verification:"
python - <<'EOF'
import scanpy as sc
from pathlib import Path

base_path = Path('/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer')
datasets = ['covid', 'hek293t', 'hct116']

print("\nEmbedding Summary:")
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
