#!/usr/bin/env python3
"""
Generate STATE embeddings for X-Atlas/Orion perturbation data (HEK293T + HCT116).

This script extracts STATE embeddings for the genome-wide Perturb-seq data
to enable perturbation-to-disease mapping analysis.

Usage:
    python generate_perturbation_embeddings.py --cell-line HEK293T
    python generate_perturbation_embeddings.py --cell-line HCT116
    python generate_perturbation_embeddings.py --cell-line all

Requirements:
    - STATE model installed (pip install arc-state)
    - GPU with sufficient memory (recommended: 24GB+)
    - STATE model weights in state_models/SE-600M/

Author: Xinru Qiu
Date: December 2025
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime

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
import pandas as pd
import scanpy as sc
import anndata as ad

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def check_state_installation():
    """Check if STATE model is installed."""
    try:
        import state
        logger.info(f"STATE version: {state.__version__ if hasattr(state, '__version__') else 'unknown'}")
        return True
    except ImportError:
        logger.error("STATE not installed. Install with: pip install arc-state")
        return False


def load_xatlas_data(cell_line: str, base_dir: str) -> ad.AnnData:
    """
    Load X-Atlas perturbation data for a specific cell line.

    Args:
        cell_line: 'HEK293T' or 'HCT116'
        base_dir: Base directory containing data

    Returns:
        AnnData object with raw counts
    """
    data_dir = Path(base_dir) / "01_DATA" / "xatlas"
    embedding_dir = Path(base_dir) / "02_EMBEDDINGS"

    # Try different possible file locations
    possible_files = [
        # Original X-Atlas raw data locations
        data_dir / f"{cell_line.lower()}_raw.h5ad",
        data_dir / f"{cell_line.upper()}_raw.h5ad",
        data_dir / f"{cell_line}_perturbation.h5ad",
        data_dir / f"xatlas_{cell_line.lower()}.h5ad",
        # UCE embedding files (contain raw counts in X matrix)
        embedding_dir / "uce" / f"{cell_line.lower()}_platelets_uce_adata.h5ad",
        embedding_dir / "UCE" / f"{cell_line.lower()}_platelets_uce_adata.h5ad",
    ]

    for f in possible_files:
        if f.exists():
            logger.info(f"Loading {cell_line} data from {f}")
            adata = sc.read_h5ad(f)
            logger.info(f"  Shape: {adata.shape}")
            logger.info(f"  Columns: {list(adata.obs.columns[:10])}")

            # If loaded from UCE file, raw counts are in X matrix
            if 'uce' in str(f).lower():
                logger.info("  Source: UCE embedding file (raw counts in X matrix)")
                # Remove UCE embeddings to save memory, we'll add STATE embeddings
                if 'X_uce' in adata.obsm:
                    logger.info("  Removing existing UCE embeddings to save memory")
                    del adata.obsm['X_uce']

            return adata

    # If no file found, check if we need to extract from combined file
    combined_file = data_dir / "xatlas_combined.h5ad"
    if combined_file.exists():
        logger.info(f"Loading {cell_line} from combined file...")
        adata = sc.read_h5ad(combined_file)

        # Filter for cell line
        if 'cell_line' in adata.obs.columns:
            adata = adata[adata.obs['cell_line'] == cell_line.upper()].copy()
        elif 'dataset' in adata.obs.columns:
            adata = adata[adata.obs['dataset'] == cell_line.upper()].copy()

        logger.info(f"  Filtered shape: {adata.shape}")
        return adata

    raise FileNotFoundError(f"Could not find X-Atlas data for {cell_line}")


def generate_state_embeddings_inference(
    input_path: str,
    output_path: str,
    model_folder: str,
    checkpoint: str = None,
    config: str = None,
    protein_embeddings: str = None,
) -> np.ndarray:
    """
    Generate STATE embeddings using the Inference API (recommended method).

    Args:
        input_path: Path to input h5ad file
        output_path: Path to save output h5ad file
        model_folder: Path to STATE model folder (e.g., SE-600M)
        checkpoint: Path to checkpoint file (optional, auto-detected)
        config: Path to config.yaml (optional, auto-detected)
        protein_embeddings: Path to protein_embeddings.pt (optional, auto-detected)

    Returns:
        Embedding matrix (n_cells x embedding_dim)
    """
    from state.emb.inference import Inference
    from omegaconf import OmegaConf

    model_folder = Path(model_folder)

    # Auto-detect paths if not provided
    if checkpoint is None:
        checkpoint = model_folder / "se600m_epoch4.ckpt"
    if config is None:
        config = model_folder / "config.yaml"
    if protein_embeddings is None:
        protein_embeddings = model_folder / "protein_embeddings.pt"

    logger.info("Loading STATE model using Inference API...")
    logger.info(f"  Model folder: {model_folder}")
    logger.info(f"  Checkpoint: {checkpoint}")
    logger.info(f"  Config: {config}")

    # Verify paths
    for path, name in [(checkpoint, "Checkpoint"), (config, "Config"),
                       (protein_embeddings, "Protein embeddings")]:
        if not Path(path).exists():
            raise FileNotFoundError(f"{name} not found: {path}")

    # Load config
    cfg = OmegaConf.load(config)

    # Load protein embeddings
    logger.info(f"  Loading protein embeddings...")
    protein_embeds = torch.load(str(protein_embeddings), map_location='cpu', weights_only=False)

    # Create inference object
    logger.info("  Creating STATE inference object...")
    inferer = Inference(cfg=cfg, protein_embeds=protein_embeds)

    # Load model checkpoint
    logger.info(f"  Loading model checkpoint...")
    inferer.load_model(str(checkpoint))
    logger.info("  Model loaded successfully!")

    # Generate embeddings using encode_adata
    logger.info(f"  Generating embeddings (this may take a while)...")
    inferer.encode_adata(input_path, output_path)
    logger.info("  Encoding complete!")

    # Load output and extract embeddings
    adata_out = sc.read_h5ad(output_path)

    # Find STATE embedding key
    embed_keys = list(adata_out.obsm.keys())
    logger.info(f"  Available obsm keys: {embed_keys}")

    state_key = None
    for key in embed_keys:
        if 'state' in key.lower() or 'emb' in key.lower():
            state_key = key
            break

    if state_key is None and embed_keys:
        state_key = embed_keys[0]

    if state_key:
        embeddings = np.array(adata_out.obsm[state_key])
        logger.info(f"  STATE embedding shape: {embeddings.shape}")
        return embeddings
    else:
        raise ValueError("No embeddings found in STATE output")


def generate_state_embeddings_simple(
    adata: ad.AnnData,
    batch_size: int = 64
) -> np.ndarray:
    """
    Simplified STATE embedding generation using command-line interface.

    This is a fallback if the Python API doesn't work as expected.
    """
    import subprocess
    import tempfile

    # Save adata to temp file
    with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as f:
        temp_input = f.name
    with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as f:
        temp_output = f.name

    try:
        adata.write_h5ad(temp_input)

        # Run STATE CLI
        cmd = [
            'state', 'embed',
            '--input', temp_input,
            '--output', temp_output,
            '--batch-size', str(batch_size)
        ]

        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"STATE CLI failed: {result.stderr}")
            raise RuntimeError(result.stderr)

        # Load output
        adata_out = sc.read_h5ad(temp_output)

        # Find embedding key
        for key in ['X_state', 'X_STATE', 'X_emb']:
            if key in adata_out.obsm:
                return np.array(adata_out.obsm[key])

        raise ValueError("No embedding found in STATE output")

    finally:
        # Cleanup
        if os.path.exists(temp_input):
            os.remove(temp_input)
        if os.path.exists(temp_output):
            os.remove(temp_output)


def main():
    parser = argparse.ArgumentParser(
        description='Generate STATE embeddings for X-Atlas perturbation data'
    )
    parser.add_argument(
        '--cell-line', type=str, required=True,
        choices=['HEK293T', 'HCT116', 'all'],
        help='Cell line to process'
    )
    parser.add_argument(
        '--base-dir', type=str,
        default="/bigdata/godziklab/shared/Xinru/302006",
        help='Base directory'
    )
    parser.add_argument(
        '--model-folder', type=str,
        default="/bigdata/godziklab/shared/Xinru/302006/state_models/SE-600M",
        help='Path to STATE model folder (SE-600M)'
    )
    parser.add_argument(
        '--batch-size', type=int, default=64,
        help='Batch size for embedding generation'
    )
    parser.add_argument(
        '--output-dir', type=str, default=None,
        help='Output directory (default: base_dir/02_EMBEDDINGS/state/)'
    )
    args = parser.parse_args()

    # Check STATE installation
    if not check_state_installation():
        sys.exit(1)

    # Setup paths
    base_dir = Path(args.base_dir)
    output_dir = Path(args.output_dir) if args.output_dir else base_dir / "02_EMBEDDINGS" / "state"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine cell lines to process
    cell_lines = ['HEK293T', 'HCT116'] if args.cell_line == 'all' else [args.cell_line]

    for cell_line in cell_lines:
        logger.info("=" * 60)
        logger.info(f"Processing {cell_line}")
        logger.info("=" * 60)

        output_file = output_dir / f"{cell_line.lower()}_state.h5ad"

        # Skip if already exists
        if output_file.exists():
            logger.info(f"Output already exists: {output_file}")
            logger.info("Skipping... (delete file to regenerate)")
            continue

        try:
            # Load data
            adata = load_xatlas_data(cell_line, args.base_dir)

            # Save data to temp file for STATE processing
            import tempfile
            with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as f:
                temp_input = f.name

            logger.info(f"Saving input to temp file: {temp_input}")
            adata.write_h5ad(temp_input)

            # Generate embeddings using Inference API (recommended method)
            try:
                embeddings = generate_state_embeddings_inference(
                    input_path=temp_input,
                    output_path=str(output_file),
                    model_folder=args.model_folder,
                )
                logger.info(f"Saved: {output_file}")

            except Exception as e:
                logger.warning(f"Inference API failed: {e}")
                logger.info("Trying simple CLI fallback...")
                embeddings = generate_state_embeddings_simple(
                    adata,
                    batch_size=args.batch_size
                )
                # Add to AnnData
                adata.obsm['X_state'] = embeddings
                # Save
                adata.write_h5ad(output_file)
                logger.info(f"Saved: {output_file}")

            # Cleanup temp file
            if os.path.exists(temp_input):
                os.remove(temp_input)

        except Exception as e:
            logger.error(f"Failed to process {cell_line}: {e}")
            continue

    logger.info("")
    logger.info("=" * 60)
    logger.info("DONE")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
