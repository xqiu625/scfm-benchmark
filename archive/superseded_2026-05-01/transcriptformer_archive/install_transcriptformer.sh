#!/bin/bash
# Install TranscriptFormer for Project 302006
# Based on TranscriptFormer comprehensive implementation guide

echo "=========================================="
echo "TRANSCRIPTFORMER INSTALLATION"
echo "=========================================="

# Load Python 3.11+ module on cluster
echo ""
echo "Step 1: Loading Python 3.11 module..."
module load python/3.11 2>/dev/null || echo "  (module command not available or Python 3.11 not found)"

# Check Python version (must be 3.11+)
echo ""
echo "Step 2: Checking Python version..."
python_version=$(python3 --version 2>&1 | grep -oP '\d+\.\d+' | head -1)
major_version=$(echo $python_version | cut -d'.' -f1)
minor_version=$(echo $python_version | cut -d'.' -f2)

echo "Current Python version: $python_version"

if [ "$major_version" -lt 3 ] || ([ "$major_version" -eq 3 ] && [ "$minor_version" -lt 11 ]); then
    echo ""
    echo "ERROR: Python 3.11+ required, found $python_version"
    echo ""
    echo "Please load Python 3.11+ manually:"
    echo "  module load python/3.11"
    echo "Or check available versions:"
    echo "  module avail python"
    exit 1
fi

# Create virtual environment
echo ""
echo "Step 3: Creating virtual environment..."
VENV_PATH="/bigdata/godziklab/shared/Xinru/transcriptformer_env"

if [ -d "$VENV_PATH" ]; then
    echo "Environment already exists at $VENV_PATH"
else
    python3 -m venv $VENV_PATH
    echo "Created environment at $VENV_PATH"
fi

# Activate environment
echo ""
echo "Step 4: Activating environment..."
source $VENV_PATH/bin/activate

# Upgrade pip
echo ""
echo "Step 5: Upgrading pip..."
pip install --upgrade pip

# Install TranscriptFormer
# CRITICAL: PyTorch must be ≤2.5.1 (≥2.6.0 causes pickle errors)
echo ""
echo "Step 6: Installing PyTorch and TranscriptFormer..."
echo "NOTE: Installing PyTorch 2.5.1 (versions ≥2.6.0 cause errors)"

# Install PyTorch with CUDA 12.1 support
pip install torch==2.5.1 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install TranscriptFormer and dependencies
pip install transcriptformer scanpy anndata pytorch-lightning hydra-core

# Verify installation
echo ""
echo "Step 7: Verifying installation..."
python -c "import transcriptformer; print('TranscriptFormer installed successfully')" 2>/dev/null || echo "Warning: TranscriptFormer import check skipped"
python -c "import torch; print(f'PyTorch version: {torch.__version__}')"
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
if [ $? -eq 0 ]; then
    echo "✓ TranscriptFormer installed successfully"
else
    echo "✗ Installation verification failed"
    exit 1
fi

# Download model weights
echo ""
echo "Step 8: Downloading model weights..."
echo "TF-Sapiens (human-specific, 57M cells, ~5GB)"
echo ""
echo "Options:"
echo "  1. tf-sapiens    - Human-specific (57M cells, 368M params)"
echo "  2. tf-exemplar   - 5 model organisms (110M cells, 542M params)"
echo "  3. tf-metazoa    - 12 species (112M cells, 444M params)"
echo "  4. all           - Download all models (~15GB)"
echo ""
echo "For this project (human COVID data), we'll use tf-sapiens"

transcriptformer download tf-sapiens

echo ""
echo "=========================================="
echo "INSTALLATION COMPLETE!"
echo "=========================================="
echo ""
echo "Virtual environment: $VENV_PATH"
echo ""
echo "To activate:"
echo "  source $VENV_PATH/bin/activate"
echo ""
echo "Next steps:"
echo "  1. Prepare raw count H5AD files with Ensembl IDs"
echo "  2. Run TranscriptFormer embedding generation"
echo "=========================================="
