#!/usr/bin/env python3
"""
Download X-Atlas/Orion Dataset from Hugging Face

Downloads the genome-wide Perturb-seq dataset containing 8 million cells
from HCT116 and HEK293T cell lines.

Dataset: https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion

Usage:
    # Download both cell lines
    python download_xatlas_orion.py --output-dir /path/to/output

    # Download only HCT116
    python download_xatlas_orion.py --cell-line HCT116 --output-dir /path/to/output

    # Download only HEK293T
    python download_xatlas_orion.py --cell-line HEK293T --output-dir /path/to/output

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
from scipy.sparse import csr_matrix, vstack
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def download_dataset(cell_line: str = None, streaming: bool = True):
    """
    Download X-Atlas/Orion dataset from Hugging Face.

    Args:
        cell_line: 'HCT116', 'HEK293T', or None for both
        streaming: If True, use streaming mode (memory efficient)

    Returns:
        Dataset object
    """
    from datasets import load_dataset

    logger.info("Loading X-Atlas/Orion dataset from Hugging Face...")
    logger.info("  Repository: Xaira-Therapeutics/X-Atlas-Orion")

    if cell_line:
        logger.info(f"  Cell line: {cell_line}")
        ds = load_dataset(
            "Xaira-Therapeutics/X-Atlas-Orion",
            split=cell_line,
            streaming=streaming
        )
    else:
        logger.info("  Cell lines: HCT116 + HEK293T (both)")
        ds = load_dataset(
            "Xaira-Therapeutics/X-Atlas-Orion",
            streaming=streaming
        )

    return ds


def download_gene_metadata():
    """Download gene metadata mapping."""
    from datasets import load_dataset

    logger.info("Loading gene metadata...")
    gene_metadata = load_dataset(
        "Xaira-Therapeutics/X-Atlas-Orion",
        "gene_metadata"
    )
    return gene_metadata


def convert_to_anndata_chunked(ds, gene_metadata, output_dir: Path, cell_line: str,
                                chunk_size: int = 500000, batch_size: int = 50000):
    """
    Convert Hugging Face dataset to AnnData format, saving in chunks.

    Saves every chunk_size cells as a separate h5ad file to avoid OOM.

    Args:
        ds: Hugging Face dataset (streaming)
        gene_metadata: Gene metadata with token to gene name mapping
        output_dir: Directory to save h5ad files
        cell_line: Cell line name for file naming
        chunk_size: Number of cells per output file
        batch_size: Number of cells to process at a time
    """
    import anndata as ad

    logger.info("Converting to AnnData format (chunked)...")
    logger.info(f"  Chunk size: {chunk_size:,} cells per file")

    # Get gene metadata as DataFrame
    if hasattr(gene_metadata, 'to_pandas'):
        gene_df = gene_metadata.to_pandas()
    else:
        gene_df = pd.DataFrame(gene_metadata['train'])

    n_genes = len(gene_df)
    logger.info(f"  Total genes: {n_genes:,}")

    # Create var DataFrame (same for all chunks)
    var_df = gene_df.set_index('gene_name')

    # Convert dataset to iterable
    if hasattr(ds, '__iter__'):
        data_iter = iter(ds)
    else:
        data_iter = iter(ds['train'] if 'train' in ds else ds)

    # Process in chunks
    chunk_idx = 0
    total_cells = 0
    chunk_files = []

    while True:
        # Collect cells for this chunk
        all_obs = []
        all_X_data = []
        all_X_indices = []
        all_X_indptr = [0]
        chunk_cells = 0

        logger.info(f"\n  Processing chunk {chunk_idx}...")

        try:
            while chunk_cells < chunk_size:
                # Get batch of cells
                batch_cells = []
                for _ in range(batch_size):
                    try:
                        row = next(data_iter)
                        batch_cells.append(row)
                    except StopIteration:
                        break

                if not batch_cells:
                    break

                # Process batch
                for cell in batch_cells:
                    obs_row = {
                        'cell_barcode': cell.get('cell_barcode', ''),
                        'sample': cell.get('sample', ''),
                        'guide_target': cell.get('guide_target', ''),
                        'gene_target': cell.get('gene_target', ''),
                        'n_genes_by_counts': cell.get('n_genes_by_counts', 0),
                        'total_counts': cell.get('total_counts', 0),
                        'total_counts_mt': cell.get('total_counts_mt', 0),
                        'pct_counts_mt': cell.get('pct_counts_mt', 0),
                        'pass_guide_filter': cell.get('pass_guide_filter', False),
                    }
                    all_obs.append(obs_row)

                    gene_tokens = cell.get('gene_token_id', [])
                    expressions = cell.get('gene_expression', [])

                    for token, expr in zip(gene_tokens, expressions):
                        all_X_indices.append(token)
                        all_X_data.append(expr)

                    all_X_indptr.append(len(all_X_data))

                chunk_cells += len(batch_cells)
                total_cells += len(batch_cells)

                if chunk_cells % 100000 == 0:
                    logger.info(f"    Processed {chunk_cells:,} cells in chunk, {total_cells:,} total...")

        except StopIteration:
            pass

        if chunk_cells == 0:
            break

        # Save this chunk
        logger.info(f"  Saving chunk {chunk_idx} ({chunk_cells:,} cells)...")

        X = csr_matrix(
            (all_X_data, all_X_indices, all_X_indptr),
            shape=(chunk_cells, n_genes)
        )

        obs_df = pd.DataFrame(all_obs)
        obs_df.index = obs_df['cell_barcode'].values

        adata = ad.AnnData(X=X, obs=obs_df, var=var_df.copy())

        chunk_path = output_dir / f"{cell_line.lower()}_xatlas_chunk{chunk_idx:02d}.h5ad"
        adata.write_h5ad(chunk_path)
        chunk_files.append(chunk_path)

        logger.info(f"  Saved: {chunk_path}")

        # Clear memory
        del adata, X, obs_df, all_obs, all_X_data, all_X_indices, all_X_indptr
        gc.collect()

        chunk_idx += 1

    logger.info(f"\n  Total cells processed: {total_cells:,}")
    logger.info(f"  Total chunks saved: {len(chunk_files)}")

    return chunk_files


def download_processed_h5ad(output_dir: Path):
    """
    Download pre-processed h5ad files from Figshare.

    These are the processed versions from the authors.
    URL: https://doi.org/10.25452/figshare.plus.29190726
    """
    logger.info("Downloading processed h5ad files from Figshare...")
    logger.info("  URL: https://doi.org/10.25452/figshare.plus.29190726")
    logger.info("")
    logger.info("  NOTE: Manual download may be required.")
    logger.info("  Visit the URL above and download the h5ad files directly.")
    logger.info(f"  Save to: {output_dir}")

    return None


def main():
    parser = argparse.ArgumentParser(
        description='Download X-Atlas/Orion dataset from Hugging Face'
    )
    parser.add_argument(
        '--cell-line', type=str, default=None,
        choices=['HCT116', 'HEK293T'],
        help='Cell line to download (default: both)'
    )
    parser.add_argument(
        '--output-dir', type=str,
        default="/bigdata/godziklab/shared/Xinru/302006/01_DATA/xatlas",
        help='Output directory'
    )
    parser.add_argument(
        '--format', type=str, default='h5ad',
        choices=['h5ad', 'parquet', 'raw'],
        help='Output format (default: h5ad)'
    )
    parser.add_argument(
        '--processed', action='store_true',
        help='Download processed h5ad from Figshare instead'
    )
    parser.add_argument(
        '--chunk-size', type=int, default=500000,
        help='Number of cells per output file (default: 500000)'
    )
    parser.add_argument(
        '--batch-size', type=int, default=50000,
        help='Batch size for processing (default: 50000)'
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 80)
    logger.info("X-ATLAS/ORION DATASET DOWNLOAD")
    logger.info("=" * 80)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cell line: {args.cell_line or 'both'}")
    logger.info(f"Format: {args.format}")
    logger.info(f"Chunk size: {args.chunk_size:,} cells")
    logger.info("")

    # Option 1: Download processed h5ad from Figshare
    if args.processed:
        download_processed_h5ad(output_dir)
        return

    # Option 2: Download from Hugging Face
    try:
        from datasets import load_dataset
    except ImportError:
        logger.error("Please install datasets: pip install datasets")
        sys.exit(1)

    # Download gene metadata first
    gene_metadata = download_gene_metadata()

    # Download dataset
    cell_lines = [args.cell_line] if args.cell_line else ['HCT116', 'HEK293T']

    for cell_line in cell_lines:
        logger.info(f"\n{'='*60}")
        logger.info(f"Downloading {cell_line}...")
        logger.info(f"{'='*60}")

        # Download with streaming to save memory
        ds = download_dataset(cell_line=cell_line, streaming=True)

        if args.format == 'h5ad':
            chunk_files = convert_to_anndata_chunked(
                ds, gene_metadata, output_dir, cell_line,
                chunk_size=args.chunk_size,
                batch_size=args.batch_size
            )
            logger.info(f"\nSaved {len(chunk_files)} chunk files for {cell_line}")

        elif args.format == 'parquet':
            output_path = output_dir / f"{cell_line.lower()}_xatlas.parquet"
            logger.info(f"Saving to parquet: {output_path}")
            # For parquet, need non-streaming
            ds_full = download_dataset(cell_line=cell_line, streaming=False)
            ds_full.to_parquet(str(output_path))

        else:  # raw
            logger.info("Keeping dataset in raw Hugging Face format")
            ds_full = download_dataset(cell_line=cell_line, streaming=False)
            output_path = output_dir / f"{cell_line.lower()}_xatlas"
            ds_full.save_to_disk(str(output_path))

    logger.info("\n" + "=" * 80)
    logger.info("DOWNLOAD COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Files saved to: {output_dir}")


if __name__ == "__main__":
    main()
