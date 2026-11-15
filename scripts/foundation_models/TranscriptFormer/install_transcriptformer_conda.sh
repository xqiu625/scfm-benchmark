#!/bin/bash
# Install TranscriptFormer using Conda (Python 3.11)
# For Project 302006

echo "=========================================="
echo "TRANSCRIPTFORMER INSTALLATION (CONDA)"
echo "=========================================="

# Create conda environment with Python 3.11
echo ""
echo "Step 1: Creating conda environment with Python 3.11..."
ENV_NAME="transcriptformer"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '${ENV_NAME}' already exists"
    echo "To remove and recreate: conda env remove -n ${ENV_NAME}"
    read -p "Remove and recreate? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n ${ENV_NAME}
        conda create -n ${ENV_NAME} python=3.11 -y
    fi
else
    conda create -n ${ENV_NAME} python=3.11 -y
fi

# Activate environment
echo ""
echo "Step 2: Activating environment..."
eval "$(conda shell.bash hook)"
conda activate ${ENV_NAME}

# Verify Python version
echo ""
echo "Step 3: Verifying Python version..."
python --version

# Install PyTorch with CUDA support
# CRITICAL: PyTorch must be ≤2.5.1 (≥2.6.0 causes pickle errors)
echo ""
echo "Step 4: Installing PyTorch 2.5.1 with CUDA 12.1..."
pip install torch==2.5.1 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install TranscriptFormer first to get its dependencies
echo ""
echo "Step 5: Installing TranscriptFormer..."
pip install transcriptformer

# Install exact versions required by TranscriptFormer 0.6.1
echo ""
echo "Step 6: Installing exact dependency versions..."
pip install numpy==2.2.6 pandas==2.2.2 scipy==1.15.3
pip install h5py==3.14 anndata==0.11.4 scanpy==1.11.2
pip install pytorch-lightning hydra-core
pip install mygene  # For Ensembl ID mapping

# Verify installation
echo ""
echo "Step 7: Verifying installation..."
echo "Checking PyTorch..."
python -c "import torch; print(f'PyTorch version: {torch.__version__}')"
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
if [ $? -eq 0 ]; then
    echo "✓ PyTorch installed successfully"
else
    echo "✗ PyTorch installation failed"
    exit 1
fi

echo ""
echo "Checking core packages..."
python -c "import scanpy; import anndata; import transcriptformer; print('All packages: OK')" 2>/dev/null
if [ $? -eq 0 ]; then
    echo "✓ All packages installed successfully"
else
    echo "✗ Package verification failed"
    exit 1
fi

# Download model weights
echo ""
echo "Step 8: Downloading TF-Sapiens model weights..."
echo "This will download ~5GB of model weights"
echo ""

# Check if transcriptformer CLI is available
if command -v transcriptformer &> /dev/null; then
    echo "Downloading model..."
    transcriptformer download tf-sapiens
else
    echo "Note: TranscriptFormer CLI not available"
    echo "You may need to download models manually or use Python API"
fi

echo ""
echo "=========================================="
echo "INSTALLATION COMPLETE!"
echo "=========================================="
echo ""
echo "Conda environment: ${ENV_NAME}"
echo ""
echo "To activate this environment:"
echo "  conda activate ${ENV_NAME}"
echo ""
echo "To deactivate:"
echo "  conda deactivate"
echo ""
echo "Next steps:"
echo "  1. Activate environment: conda activate ${ENV_NAME}"
echo "  2. Prepare data: python scripts/foundation_models/TranscriptFormer/prepare_data_for_transcriptformer.py"
echo "  3. Generate embeddings: python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py"
echo ""
echo "=========================================="
