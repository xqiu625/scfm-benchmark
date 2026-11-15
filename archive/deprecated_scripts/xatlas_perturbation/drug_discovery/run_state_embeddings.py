#!/usr/bin/env python3
"""
Generate STATE embeddings for COVID severity data.

This script patches PyTorch to work with PyTorch 2.6+ weights_only security changes,
then uses the STATE Python API directly.

Usage:
    python run_state_embeddings.py

Run on GPU node:
    srun -p short_gpu --gres=gpu:1 --mem=128gb -c 4 --time=2:00:00 --pty bash -l
    conda activate /bigdata/godziklab/shared/Xinru/302006/state_env
    cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
    python scripts/drug_discovery/state_model/run_state_embeddings.py

Alternative (requires PyTorch 2.5.x):
    pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121
    state emb transform --model-folder /path/to/SE-600M --checkpoint checkpoint.ckpt --input input.h5ad --output output.h5ad
"""

import os
import sys

# =============================================================================
# CRITICAL: Patch torch.load BEFORE any other imports
# This fixes PyTorch 2.6+ weights_only=True default issue
# =============================================================================
import torch

_original_torch_load = torch.load

def _patched_torch_load(f, *args, **kwargs):
    """Force weights_only=False for all torch.load calls."""
    kwargs['weights_only'] = False
    return _original_torch_load(f, *args, **kwargs)

# Apply patch globally
torch.load = _patched_torch_load
print("Patched torch.load to use weights_only=False")

# Also patch Lightning's internal loader
try:
    import lightning.fabric.utilities.cloud_io as cloud_io

    def _patched_pl_load(path, map_location=None, weights_only=None):
        return _original_torch_load(path, map_location=map_location, weights_only=False)

    cloud_io._load = _patched_pl_load
    print("Patched Lightning cloud_io._load")
except ImportError:
    print("Lightning not found, skipping patch")

# Patch lightning.pytorch.core.saving if it exists
try:
    import lightning.pytorch.core.saving as pl_saving
    original_pl_load = pl_saving.pl_load if hasattr(pl_saving, 'pl_load') else None

    if original_pl_load:
        def _patched_saving_load(path, map_location=None, weights_only=None):
            return _original_torch_load(path, map_location=map_location, weights_only=False)
        pl_saving.pl_load = _patched_saving_load
        print("Patched lightning.pytorch.core.saving.pl_load")
except ImportError:
    pass

# =============================================================================
# Now import everything else
# =============================================================================
import numpy as np
import scanpy as sc
from pathlib import Path

# Paths
MODEL_FOLDER = "/bigdata/godziklab/shared/Xinru/302006/state_models/SE-600M"
CHECKPOINT = os.path.join(MODEL_FOLDER, "se600m_epoch4.ckpt")
CONFIG = os.path.join(MODEL_FOLDER, "config.yaml")
PROTEIN_EMBEDDINGS = os.path.join(MODEL_FOLDER, "protein_embeddings.pt")

INPUT_FILE = "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/uce/human_platelet_covid_severity_uce_adata.h5ad"
OUTPUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/state"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "covid_state.h5ad")

BATCH_SIZE = 32
EMBED_KEY = "X_state"


def main():
    print("=" * 60)
    print("STATE EMBEDDING GENERATION")
    print("=" * 60)

    # Check CUDA
    print(f"\nPyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"CUDA device: {torch.cuda.get_device_name(0)}")

    # Verify paths
    for path, name in [(MODEL_FOLDER, "Model folder"), (CHECKPOINT, "Checkpoint"),
                       (CONFIG, "Config"), (PROTEIN_EMBEDDINGS, "Protein embeddings"),
                       (INPUT_FILE, "Input file")]:
        if not os.path.exists(path):
            print(f"ERROR: {name} not found: {path}")
            sys.exit(1)

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load input data
    print(f"\nLoading input data: {INPUT_FILE}")
    adata = sc.read_h5ad(INPUT_FILE)
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")

    # Import STATE modules
    print("\nImporting STATE...")
    try:
        from state.emb.inference import Inference
        from omegaconf import OmegaConf
    except ImportError as e:
        print(f"ERROR: Could not import STATE: {e}")
        print("Make sure arc-state is installed: pip install arc-state")
        sys.exit(1)

    # Load config
    print(f"\nLoading config: {CONFIG}")
    cfg = OmegaConf.load(CONFIG)

    # Load protein embeddings
    print(f"Loading protein embeddings: {PROTEIN_EMBEDDINGS}")
    protein_embeds = torch.load(PROTEIN_EMBEDDINGS, map_location='cpu', weights_only=False)

    # Create inference object
    print("\nCreating STATE inference object...")
    inferer = Inference(cfg=cfg, protein_embeds=protein_embeds)

    # Load model checkpoint
    print(f"Loading model checkpoint: {CHECKPOINT}")
    inferer.load_model(CHECKPOINT)
    print("  Model loaded successfully!")

    # Check available methods
    methods = [m for m in dir(inferer) if not m.startswith('_') and callable(getattr(inferer, m, None))]
    print(f"\nAvailable inference methods: {methods}")

    # Generate embeddings
    print(f"\nGenerating embeddings...")
    print("  This may take a while for 47,977 cells...")

    # encode_adata expects input path and output path
    # It writes embeddings directly to output file
    if hasattr(inferer, 'encode_adata'):
        print(f"  Using inferer.encode_adata('{INPUT_FILE}', '{OUTPUT_FILE}')...")
        inferer.encode_adata(INPUT_FILE, OUTPUT_FILE)
        print("  Encoding complete!")

        # Load the output to verify and get embedding info
        print(f"\nLoading output to verify...")
        adata_out = sc.read_h5ad(OUTPUT_FILE)

        # Find the embedding key (STATE might use different name)
        embed_keys = list(adata_out.obsm.keys())
        print(f"  Available obsm keys: {embed_keys}")

        if embed_keys:
            # Use the first embedding key found
            embed_key = embed_keys[0] if EMBED_KEY not in embed_keys else EMBED_KEY
            embeddings = adata_out.obsm[embed_key]
            print(f"  Embedding shape: {embeddings.shape}")
            print(f"  Embedding key: {embed_key}")
        else:
            print("  WARNING: No embeddings found in obsm!")
            embeddings = None
    else:
        print("  ERROR: encode_adata method not found!")
        sys.exit(1)

    # Summary
    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Input: {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    if embeddings is not None:
        print(f"Cells: {embeddings.shape[0]}")
        print(f"Embedding dim: {embeddings.shape[1]}")


if __name__ == "__main__":
    main()
