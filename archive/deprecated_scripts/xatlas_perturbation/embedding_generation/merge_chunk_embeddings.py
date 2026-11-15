#!/usr/bin/env python3
"""
Merge chunked embeddings into single files per cell line.

Run after all chunks are processed:
    python merge_chunk_embeddings.py --model STATE
    python merge_chunk_embeddings.py --model UCE

Author: Project 302006
Date: December 2025
"""

import argparse
import logging
from pathlib import Path
import anndata as ad
import gc

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"


def merge_chunks(model: str, cell_line: str):
    """Merge all chunks for a cell line."""
    chunk_dir = EMB_DIR / model.lower() / "xatlas_chunks"
    pattern = f"{cell_line.lower()}_chunk*_{model.lower()}.h5ad"

    chunks = sorted(chunk_dir.glob(pattern))

    if not chunks:
        logger.warning(f"No chunks found for {cell_line} {model}")
        return None

    logger.info(f"Found {len(chunks)} chunks for {cell_line}")

    # Load all chunks
    adatas = []
    for chunk_path in chunks:
        logger.info(f"  Loading: {chunk_path.name}")
        adata = ad.read_h5ad(chunk_path)
        logger.info(f"    Cells: {adata.n_obs:,}")
        adatas.append(adata)
        gc.collect()

    # Concatenate
    logger.info("Concatenating...")
    if len(adatas) == 1:
        merged = adatas[0]
    else:
        merged = ad.concat(adatas, join='outer')

    logger.info(f"  Total cells: {merged.n_obs:,}")

    # Save
    output_path = EMB_DIR / model.lower() / f"{cell_line.lower()}_xatlas_full_{model.lower()}.h5ad"
    logger.info(f"Saving to: {output_path}")
    merged.write_h5ad(output_path)

    del adatas, merged
    gc.collect()

    return output_path


def main():
    parser = argparse.ArgumentParser(description='Merge chunk embeddings')
    parser.add_argument('--model', type=str, required=True, choices=['STATE', 'UCE'])
    parser.add_argument('--cell-lines', type=str, nargs='+', default=['HCT116', 'HEK293T'])
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info(f"MERGING {args.model} CHUNK EMBEDDINGS")
    logger.info("=" * 60)

    for cell_line in args.cell_lines:
        logger.info(f"\n{cell_line}:")
        merge_chunks(args.model, cell_line)

    logger.info("\nDone!")


if __name__ == "__main__":
    main()
