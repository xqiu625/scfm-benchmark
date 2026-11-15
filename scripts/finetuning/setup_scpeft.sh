#!/bin/bash
# =============================================================================
# Setup scPEFT for scGPT/Geneformer PEFT Fine-Tuning
#
# This script:
#   1. Clones the scPEFT repo (modified scGPT/Geneformer with adapter support)
#   2. Creates conda environment (Python 3.9, PyTorch 2.1+cu118)
#   3. Downloads scGPT whole-human pretrained checkpoint
#   4. Verifies installation
#
# Usage (on interactive GPU node):
#   bash setup_scpeft.sh
#
# After setup:
#   source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
#   cd /bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT
#   # Run tutorials or custom scripts
# =============================================================================

set -e

MODELS_DIR="/bigdata/godziklab/shared/Xinru/302006/03_MODELS"
SCPEFT_DIR="${MODELS_DIR}/scPEFT"
ENV_PREFIX="/bigdata/godziklab/shared/Xinru/302006/scpeft_env"

echo "============================================================"
echo "scPEFT Setup for COVID Severity Classification"
echo "============================================================"

# ---------------------------------------------------------
# Step 1: Clone scPEFT repo
# ---------------------------------------------------------
echo ""
echo "[Step 1/4] Cloning scPEFT repository..."

if [ -d "${SCPEFT_DIR}" ]; then
    echo "scPEFT directory already exists at ${SCPEFT_DIR}"
    echo "Pulling latest changes..."
    cd "${SCPEFT_DIR}" && git pull
else
    cd "${MODELS_DIR}"
    git clone https://github.com/coffee19850519/scPEFT.git
fi

echo "scPEFT repo ready at ${SCPEFT_DIR}"

# ---------------------------------------------------------
# Step 2: Create conda environment
# ---------------------------------------------------------
echo ""
echo "[Step 2/4] Setting up conda environment: ${ENV_PREFIX}"

if [ -d "${ENV_PREFIX}" ]; then
    echo "Environment already exists at ${ENV_PREFIX}. Skipping creation."
    echo "To recreate: conda env remove -p ${ENV_PREFIX} -y && bash setup_scpeft.sh"
else
    conda create --prefix ${ENV_PREFIX} python=3.9 -y
fi

source activate ${ENV_PREFIX}

# Install PyTorch 2.1 (matches scPEFT requirements.yaml)
echo "Installing PyTorch 2.1.0+cu118..."
pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118

# Install torchtext (required by scPEFT's bundled scGPT)
echo "Installing torchtext 0.16.0 (required by scGPT)..."
pip install torchtext==0.16.0 --index-url https://download.pytorch.org/whl/cu118

# Install key dependencies matching scPEFT requirements
echo "Installing scPEFT dependencies..."
pip install transformers==4.45.2
pip install peft==0.4.0
pip install accelerate==0.21.0
pip install "numpy<2" scipy scanpy scikit-learn loralib
pip install ipython flash-attn 2>/dev/null || echo "flash-attn install failed (optional, will use standard attention)"
pip install scvi-tools 2>/dev/null || echo "scvi-tools install failed (optional)"

# ---------------------------------------------------------
# Step 3: Download scGPT whole-human pretrained checkpoint
# ---------------------------------------------------------
echo ""
echo "[Step 3/4] Setting up scGPT pretrained model..."

SAVE_DIR="${SCPEFT_DIR}/save/scGPT_human"
mkdir -p "${SAVE_DIR}"

if [ -f "${SAVE_DIR}/best_model.pt" ]; then
    echo "scGPT checkpoint already exists at ${SAVE_DIR}"
else
    echo ""
    echo "============================================================"
    echo "MANUAL DOWNLOAD REQUIRED"
    echo "============================================================"
    echo ""
    echo "The scGPT whole-human checkpoint must be downloaded from Google Drive:"
    echo ""
    echo "  https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y"
    echo ""
    echo "Download these 3 files and place them in:"
    echo "  ${SAVE_DIR}/"
    echo ""
    echo "  - best_model.pt   (model weights)"
    echo "  - args.json        (model config)"
    echo "  - vocab.json       (gene vocabulary)"
    echo ""
    echo "Alternative: Download via AWS S3 (no authentication needed):"
    echo "  pip install awscli"
    echo "  aws s3 sync --no-sign-request s3://czi-scgenept-public/models/pretrained/scgpt ${SAVE_DIR}/"
    echo ""
    echo "Attempting AWS S3 download..."
    pip install awscli -q
    aws s3 sync --no-sign-request s3://czi-scgenept-public/models/pretrained/scgpt "${SAVE_DIR}/" 2>/dev/null || {
        echo "AWS S3 download failed. Trying gdown for Google Drive..."
        pip install gdown -q
        # Google Drive folder ID: 1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y
        gdown --folder "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y" -O "${SAVE_DIR}/" 2>/dev/null || {
            echo ""
            echo "WARNING: Automatic download failed."
            echo "Please manually download from the Google Drive link above."
        }
    }
fi

# Also symlink for easy access
ln -sf "${SAVE_DIR}" "${MODELS_DIR}/scgpt/scGPT_human" 2>/dev/null || true

# ---------------------------------------------------------
# Step 4: Verify installation
# ---------------------------------------------------------
echo ""
echo "[Step 4/4] Verifying installation..."
echo "============================================================"

# Add scPEFT to Python path and fix libstdc++ for conda env
export PYTHONPATH="${SCPEFT_DIR}:${PYTHONPATH}"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH}"

python -c "
import sys
sys.path.insert(0, '${SCPEFT_DIR}')

import torch
print(f'PyTorch: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'GPU: {torch.cuda.get_device_name(0)}')

import transformers
print(f'transformers: {transformers.__version__}')

import scanpy
print(f'scanpy: {scanpy.__version__}')

# Check scPEFT's bundled scGPT
try:
    from scgpt.model import TransformerModel
    print('scGPT TransformerModel: OK')
except Exception as e:
    print(f'scGPT TransformerModel: FAILED ({e})')

try:
    from scgpt.utils import PeftConfig
    print('PeftConfig: OK')
except Exception as e:
    print(f'PeftConfig: FAILED ({e})')

try:
    from scgpt.tokenizer.gene_tokenizer import GeneVocab
    print('GeneVocab: OK')
except Exception as e:
    print(f'GeneVocab: FAILED ({e})')

# Check checkpoint
import os
save_dir = '${SAVE_DIR}'
for f in ['best_model.pt', 'args.json', 'vocab.json']:
    path = os.path.join(save_dir, f)
    if os.path.exists(path):
        size_mb = os.path.getsize(path) / 1024 / 1024
        print(f'{f}: OK ({size_mb:.1f} MB)')
    else:
        print(f'{f}: MISSING')

print()
print('All checks complete!')
"

echo ""
echo "============================================================"
echo "Setup complete!"
echo "============================================================"
echo ""
echo "To use scPEFT:"
echo "  source activate ${ENV_PREFIX}"
echo "  export PYTHONPATH=${SCPEFT_DIR}:\$PYTHONPATH"
echo "  export LD_LIBRARY_PATH=${ENV_PREFIX}/lib:\$LD_LIBRARY_PATH"
echo "  cd ${SCPEFT_DIR}"
echo ""
echo "Run scGPT PEFT tutorial:"
echo "  cd tutorial_peft"
echo "  jupyter notebook Tutorial_Identification.ipynb"
echo ""
echo "Or run custom COVID severity script:"
echo "  cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning"
echo "  python scgpt_peft_finetune.py --peft LORA --epochs 20 --all-tasks"
echo ""
