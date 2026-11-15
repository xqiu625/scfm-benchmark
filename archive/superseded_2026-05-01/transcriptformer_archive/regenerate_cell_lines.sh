#!/bin/bash
# Regenerate embeddings for HEK293T and HCT116 (duplicates already removed)
# Data has been fixed, just need to regenerate embeddings

echo "=========================================="
echo "REGENERATE CELL LINE EMBEDDINGS"
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

# Generate embeddings for HEK293T
echo ""
echo "Step 3: Generating TranscriptFormer embeddings for HEK293T..."
python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --dataset hek293t --batch-size 16

if [ $? -ne 0 ]; then
    echo "ERROR: HEK293T embedding generation failed"
    exit 1
fi

# Generate embeddings for HCT116
echo ""
echo "Step 4: Generating TranscriptFormer embeddings for HCT116..."
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
echo "  3. Evaluate TranscriptFormer performance"
echo "=========================================="
