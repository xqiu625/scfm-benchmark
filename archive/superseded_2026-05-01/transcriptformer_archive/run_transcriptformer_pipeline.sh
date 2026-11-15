#!/bin/bash
# Complete TranscriptFormer Pipeline
# 1. Re-prepare data with required metadata (assay, organism, tissue)
# 2. Generate embeddings using TF-Sapiens model

echo "=========================================="
echo "TRANSCRIPTFORMER COMPLETE PIPELINE"
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
python -c "import torch; print(f'PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

# Re-prepare data with required metadata
echo ""
echo "Step 3: Re-preparing data with required metadata..."
echo "This will add: assay, organism, tissue columns"
python scripts/foundation_models/TranscriptFormer/prepare_data_for_transcriptformer.py

if [ $? -ne 0 ]; then
    echo "ERROR: Data preparation failed"
    exit 1
fi

# Generate embeddings
echo ""
echo "Step 4: Generating TranscriptFormer embeddings..."
echo "This will process all three datasets using TF-Sapiens model"
python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --batch-size 16

if [ $? -ne 0 ]; then
    echo "ERROR: Embedding generation failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Data: /bigdata/godziklab/shared/Xinru/302006/01_DATA/transcriptformer_ready/"
echo "  Embeddings: /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/"
echo ""
echo "Next steps:"
echo "  1. Verify embedding shapes and quality"
echo "  2. Merge with UCE and scGPT embeddings"
echo "  3. Run three-model comparison analysis"
echo "=========================================="
