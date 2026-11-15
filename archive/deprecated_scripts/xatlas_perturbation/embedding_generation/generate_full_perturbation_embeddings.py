#!/usr/bin/env python3
"""
Generate Foundation Model Embeddings for ALL Perturbation Cells

This script generates STATE or UCE embeddings for the complete HEK293T and HCT116
perturbation datasets from X-Atlas/Orion (8 million cells).

Supports loading chunked h5ad files from X-Atlas download.

Usage:
    python generate_full_perturbation_embeddings.py --model STATE --cell-lines HEK293T HCT116
    python generate_full_perturbation_embeddings.py --model UCE --cell-lines HEK293T
    python generate_full_perturbation_embeddings.py --model STATE --max-cells 500000

Author: Project 302006
Date: December 2025
"""

import os
import sys
import argparse
import logging
import gc
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
DATA_DIR = BASE_DIR / "01_DATA"
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"

# Raw data locations (try multiple possible paths, including X-Atlas chunks)
RAW_DATA_PATHS = {
    'HEK293T': [
        # X-Atlas chunked files (primary source - 4M cells)
        DATA_DIR / "xatlas" / "hek293t_xatlas_chunk*.h5ad",  # Glob pattern
        DATA_DIR / "xatlas" / "HEK293T_raw.h5ad",
        DATA_DIR / "xatlas" / "hek293t_raw.h5ad",
        DATA_DIR / "raw" / "HEK293T.h5ad",
        DATA_DIR / "HEK293T_subset_150k_quality.h5ad",
        BASE_DIR / "302006_git" / "data" / "HEK293T_subset_150k_quality.h5ad",
    ],
    'HCT116': [
        # X-Atlas chunked files (primary source - 4M cells)
        DATA_DIR / "xatlas" / "hct116_xatlas_chunk*.h5ad",  # Glob pattern
        DATA_DIR / "xatlas" / "HCT116_raw.h5ad",
        DATA_DIR / "xatlas" / "hct116_raw.h5ad",
        DATA_DIR / "raw" / "HCT116.h5ad",
        DATA_DIR / "HCT116_subset_150k_quality.h5ad",
        BASE_DIR / "302006_git" / "data" / "HCT116_subset_150k_quality.h5ad",
    ]
}


def find_raw_data(cell_line: str):
    """
    Find the raw data file(s) for a cell line.

    Returns:
        tuple: (path_or_paths, is_chunked)
               - For single file: (Path, False)
               - For chunked files: (list[Path], True)
    """
    for path in RAW_DATA_PATHS.get(cell_line, []):
        path_str = str(path)

        # Handle glob patterns for chunked files
        if '*' in path_str:
            parent = path.parent
            pattern = path.name
            if parent.exists():
                chunk_files = sorted(parent.glob(pattern))
                if chunk_files:
                    logger.info(f"Found {len(chunk_files)} X-Atlas chunks for {cell_line}")
                    for f in chunk_files:
                        logger.info(f"  - {f.name}")
                    return chunk_files, True
        else:
            # Regular single file
            if path.exists():
                logger.info(f"Found {cell_line} data at: {path}")
                return path, False

    # List available files if not found
    logger.error(f"Could not find raw data for {cell_line}")
    logger.info("Searched paths:")
    for path in RAW_DATA_PATHS.get(cell_line, []):
        logger.info(f"  {path}")

    # Try to find any h5ad files
    for search_dir in [DATA_DIR, DATA_DIR / "xatlas", DATA_DIR / "raw"]:
        if search_dir.exists():
            h5ad_files = list(search_dir.glob("*.h5ad"))
            if h5ad_files:
                logger.info(f"Available h5ad files in {search_dir}:")
                for f in h5ad_files[:10]:
                    logger.info(f"  {f.name}")

    return None, False


def load_chunked_data(chunk_files: list, max_cells: int = None, sample_strategy: str = 'balanced') -> ad.AnnData:
    """
    Load and concatenate chunked h5ad files.

    Args:
        chunk_files: List of chunk file paths
        max_cells: Maximum number of cells to load (None = all)
        sample_strategy: 'balanced' (equal from each chunk) or 'sequential' (first N cells)

    Returns:
        Concatenated AnnData object
    """
    logger.info(f"Loading {len(chunk_files)} chunk files...")

    if max_cells is not None:
        logger.info(f"  Max cells: {max_cells:,}")
        logger.info(f"  Strategy: {sample_strategy}")

        if sample_strategy == 'balanced':
            cells_per_chunk = max_cells // len(chunk_files)
            logger.info(f"  Cells per chunk: ~{cells_per_chunk:,}")

    adatas = []
    total_cells = 0

    for i, chunk_path in enumerate(chunk_files):
        logger.info(f"  Loading chunk {i}: {chunk_path.name}")
        adata = sc.read_h5ad(chunk_path)
        logger.info(f"    Cells: {adata.n_obs:,}")

        # Apply sampling if needed
        if max_cells is not None and sample_strategy == 'balanced':
            if adata.n_obs > cells_per_chunk:
                # Random sample
                np.random.seed(42 + i)
                indices = np.random.choice(adata.n_obs, cells_per_chunk, replace=False)
                adata = adata[indices].copy()
                logger.info(f"    Sampled: {adata.n_obs:,}")

        adatas.append(adata)
        total_cells += adata.n_obs

        # Check total cell limit
        if max_cells is not None and sample_strategy == 'sequential':
            if total_cells >= max_cells:
                logger.info(f"  Reached max cells ({max_cells:,}), stopping...")
                break

        # Clear memory
        gc.collect()

    # Make var_names unique in each chunk before concatenating
    logger.info("Making gene names unique in each chunk...")
    for i, adata in enumerate(adatas):
        if not adata.var_names.is_unique:
            logger.info(f"  Chunk {i}: Making {adata.n_vars} var_names unique")
            adata.var_names_make_unique()

    # Concatenate
    logger.info("Concatenating chunks...")
    if len(adatas) == 1:
        combined = adatas[0]
    else:
        # Use inner join to only keep common genes
        combined = ad.concat(adatas, join='inner')

    # Final sampling for sequential strategy
    if max_cells is not None and combined.n_obs > max_cells:
        np.random.seed(42)
        indices = np.random.choice(combined.n_obs, max_cells, replace=False)
        combined = combined[indices].copy()

    logger.info(f"  Total cells: {combined.n_obs:,}")
    logger.info(f"  Total genes: {combined.n_vars:,}")

    # Clean up
    del adatas
    gc.collect()

    return combined


def load_perturbation_data(cell_line: str, max_cells: int = None) -> ad.AnnData:
    """
    Load complete perturbation dataset.

    Args:
        cell_line: 'HEK293T' or 'HCT116'
        max_cells: Maximum cells to load (for memory management)

    Returns:
        AnnData object with perturbation data
    """
    path_or_paths, is_chunked = find_raw_data(cell_line)

    if path_or_paths is None:
        raise FileNotFoundError(f"No data found for {cell_line}")

    if is_chunked:
        # Load chunked X-Atlas files
        adata = load_chunked_data(path_or_paths, max_cells=max_cells)
    else:
        # Load single file
        logger.info(f"Loading {cell_line} data from {path_or_paths}")
        adata = sc.read_h5ad(path_or_paths)

        # Apply max_cells sampling if needed
        if max_cells is not None and adata.n_obs > max_cells:
            logger.info(f"Sampling {max_cells:,} cells from {adata.n_obs:,}")
            np.random.seed(42)
            indices = np.random.choice(adata.n_obs, max_cells, replace=False)
            adata = adata[indices].copy()

    logger.info(f"  Total cells: {adata.n_obs:,}")
    logger.info(f"  Total genes: {adata.n_vars:,}")
    logger.info(f"  Columns: {list(adata.obs.columns)[:10]}")

    # Find perturbation column
    pert_col = None
    for col in ['guide_target', 'gene_target', 'perturbation', 'target_gene']:
        if col in adata.obs.columns:
            pert_col = col
            break

    if pert_col:
        n_perts = adata.obs[pert_col].nunique()
        logger.info(f"  Perturbation column: {pert_col}")
        logger.info(f"  Unique perturbations: {n_perts:,}")

    return adata


def generate_state_embeddings(adata: ad.AnnData, output_path: Path):
    """Generate STATE embeddings for the data."""
    import torch

    # Patch torch.load for STATE compatibility
    _original_torch_load = torch.load
    def _patched_torch_load(f, *args, **kwargs):
        kwargs['weights_only'] = False
        return _original_torch_load(f, *args, **kwargs)
    torch.load = _patched_torch_load

    try:
        from state.emb.inference import Inference
        from omegaconf import OmegaConf
    except ImportError:
        logger.error("STATE not installed. Install with: pip install arc-state")
        return None

    # STATE model paths
    model_folder = BASE_DIR / "state_models" / "SE-600M"
    checkpoint = model_folder / "se600m_epoch4.ckpt"
    config = model_folder / "config.yaml"
    protein_embeddings = model_folder / "protein_embeddings.pt"

    if not all(p.exists() for p in [checkpoint, config, protein_embeddings]):
        logger.error("STATE model files not found")
        return None

    # Save to temp file
    import tempfile
    with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as f:
        temp_input = f.name

    logger.info("Saving data for STATE processing...")
    adata.write_h5ad(temp_input)

    # Load model
    logger.info("Loading STATE model...")
    cfg = OmegaConf.load(config)
    protein_embeds = torch.load(str(protein_embeddings), map_location='cpu', weights_only=False)
    inferer = Inference(cfg=cfg, protein_embeds=protein_embeds)
    inferer.load_model(str(checkpoint))

    # Generate embeddings
    logger.info("Generating STATE embeddings (this may take a while)...")
    inferer.encode_adata(temp_input, str(output_path))

    # Cleanup
    if os.path.exists(temp_input):
        os.remove(temp_input)

    logger.info(f"Saved STATE embeddings to: {output_path}")
    return output_path


def generate_uce_embeddings(adata: ad.AnnData, output_path: Path, cell_line: str):
    """Generate UCE embeddings for the data."""
    import subprocess
    import tempfile

    # UCE model path
    uce_dir = Path("/bigdata/godziklab/shared/Xinru/UCE")

    if not uce_dir.exists():
        logger.error(f"UCE directory not found: {uce_dir}")
        return None

    # Save to temp file
    with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as f:
        temp_input = f.name

    logger.info("Saving data for UCE processing...")
    adata.write_h5ad(temp_input)

    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Run UCE
    logger.info("Generating UCE embeddings...")
    cmd = [
        'python', str(uce_dir / 'eval_single_anndata.py'),
        '--adata_path', temp_input,
        '--dir', str(output_path.parent),
        '--species', 'human',
        '--batch_size', '200',
    ]

    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(uce_dir))

    if result.returncode != 0:
        logger.error(f"UCE failed: {result.stderr}")
        # Cleanup
        if os.path.exists(temp_input):
            os.remove(temp_input)
        return None

    # Find and rename output file
    # UCE typically outputs to the same directory with modified name
    possible_outputs = list(output_path.parent.glob("*.h5ad"))
    for out_file in possible_outputs:
        if 'uce' in out_file.name.lower() and out_file != output_path:
            out_file.rename(output_path)
            break

    # Cleanup
    if os.path.exists(temp_input):
        os.remove(temp_input)

    logger.info(f"Saved UCE embeddings to: {output_path}")
    return output_path


def main():
    parser = argparse.ArgumentParser(description='Generate embeddings for all perturbation cells')
    parser.add_argument('--model', type=str, required=True, choices=['STATE', 'UCE'],
                        help='Foundation model to use')
    parser.add_argument('--cell-lines', type=str, nargs='+', default=['HEK293T', 'HCT116'],
                        help='Cell lines to process')
    parser.add_argument('--max-cells', type=int, default=None,
                        help='Maximum cells per cell line (for memory/speed, default: all)')
    parser.add_argument('--base-dir', type=str, default="/bigdata/godziklab/shared/Xinru/302006",
                        help='Base directory')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing output files')
    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info(f"GENERATING {args.model} EMBEDDINGS FOR PERTURBATION CELLS")
    logger.info("=" * 80)
    logger.info(f"Cell lines: {args.cell_lines}")
    logger.info(f"Max cells: {args.max_cells if args.max_cells else 'ALL'}")

    # Output directory
    if args.model == 'STATE':
        output_dir = EMB_DIR / "state"
    else:
        output_dir = EMB_DIR / "uce"
    output_dir.mkdir(parents=True, exist_ok=True)

    for cell_line in args.cell_lines:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {cell_line}")
        logger.info(f"{'='*60}")

        # Output filename includes max_cells if specified
        if args.max_cells:
            output_file = output_dir / f"{cell_line.lower()}_xatlas_{args.max_cells//1000}k_{args.model.lower()}.h5ad"
        else:
            output_file = output_dir / f"{cell_line.lower()}_xatlas_full_{args.model.lower()}.h5ad"

        if output_file.exists() and not args.force:
            logger.info(f"Output already exists: {output_file}")
            logger.info("Skipping... (use --force to regenerate)")
            continue

        try:
            # Load data
            adata = load_perturbation_data(cell_line, max_cells=args.max_cells)

            # Generate embeddings
            if args.model == 'STATE':
                generate_state_embeddings(adata, output_file)
            else:
                generate_uce_embeddings(adata, output_file, cell_line)

            # Clear memory
            del adata
            gc.collect()

        except Exception as e:
            logger.error(f"Failed to process {cell_line}: {e}")
            import traceback
            traceback.print_exc()
            continue

    logger.info("\n" + "=" * 80)
    logger.info("DONE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
