#!/usr/bin/env python3
"""
Process a SINGLE X-Atlas chunk for embedding generation.

Designed for short_gpu partition (2-hour limit):
- Processes ONE chunk at a time
- Saves embeddings to separate file per chunk
- Can be resubmitted for next chunk

Usage:
    python process_single_chunk.py --model STATE --cell-line HCT116 --chunk 0
    python process_single_chunk.py --model STATE --cell-line HEK293T --chunk 5
    python process_single_chunk.py --model STATE --auto  # Auto-find next unprocessed chunk

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
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
DATA_DIR = BASE_DIR / "01_DATA" / "xatlas"
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"


def get_all_chunks(cell_line: str):
    """Get list of all chunk files for a cell line."""
    pattern = f"{cell_line.lower()}_xatlas_chunk*.h5ad"
    chunks = sorted(DATA_DIR.glob(pattern))
    return chunks


def get_output_path(cell_line: str, chunk_idx: int, model: str):
    """Get output path for a chunk's embeddings."""
    model_dir = EMB_DIR / model.lower() / "xatlas_chunks"
    model_dir.mkdir(parents=True, exist_ok=True)
    return model_dir / f"{cell_line.lower()}_chunk{chunk_idx:02d}_{model.lower()}.h5ad"


def find_next_chunk(model: str):
    """Find the next unprocessed chunk across all cell lines."""
    for cell_line in ['HCT116', 'HEK293T']:
        chunks = get_all_chunks(cell_line)
        for i, chunk_path in enumerate(chunks):
            output_path = get_output_path(cell_line, i, model)
            if not output_path.exists():
                return cell_line, i, chunk_path
    return None, None, None


def get_progress_summary(model: str):
    """Get summary of processed/remaining chunks."""
    summary = {}
    for cell_line in ['HCT116', 'HEK293T']:
        chunks = get_all_chunks(cell_line)
        processed = 0
        for i in range(len(chunks)):
            output_path = get_output_path(cell_line, i, model)
            if output_path.exists():
                processed += 1
        summary[cell_line] = {'total': len(chunks), 'processed': processed}
    return summary


def generate_state_embeddings(input_path: Path, output_path: Path):
    """Generate STATE embeddings for a single chunk."""
    import torch
    import scanpy as sc

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
        return False

    # STATE model paths
    model_folder = BASE_DIR / "state_models" / "SE-600M"
    checkpoint = model_folder / "se600m_epoch4.ckpt"
    config = model_folder / "config.yaml"
    protein_embeddings = model_folder / "protein_embeddings.pt"

    if not all(p.exists() for p in [checkpoint, config, protein_embeddings]):
        logger.error("STATE model files not found")
        return False

    # Load model
    logger.info("Loading STATE model...")
    cfg = OmegaConf.load(config)
    protein_embeds = torch.load(str(protein_embeddings), map_location='cuda', weights_only=False)
    inferer = Inference(cfg=cfg, protein_embeds=protein_embeds)
    inferer.load_model(str(checkpoint))

    # Move to GPU if available
    if torch.cuda.is_available():
        logger.info(f"Using GPU: {torch.cuda.get_device_name(0)}")
        inferer.model = inferer.model.cuda()

    # Generate embeddings
    logger.info(f"Generating STATE embeddings for: {input_path.name}")
    logger.info(f"Output: {output_path}")

    inferer.encode_adata(str(input_path), str(output_path))

    logger.info(f"Saved embeddings to: {output_path}")
    return True


def generate_uce_embeddings(input_path: Path, output_path: Path):
    """Generate UCE embeddings for a single chunk."""
    import subprocess
    import shutil

    # UCE model path
    uce_dir = Path("/bigdata/godziklab/shared/Xinru/UCE")

    if not uce_dir.exists():
        logger.error(f"UCE directory not found: {uce_dir}")
        return False

    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Run UCE
    logger.info(f"Generating UCE embeddings for: {input_path.name}")
    cmd = [
        'python', str(uce_dir / 'eval_single_anndata.py'),
        '--adata_path', str(input_path),
        '--dir', str(output_path.parent),
        '--species', 'human',
        '--batch_size', '200',
    ]

    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(uce_dir))

    if result.returncode != 0:
        logger.error(f"UCE failed: {result.stderr}")
        return False

    # Find and rename output file
    temp_output = output_path.parent / input_path.name.replace('.h5ad', '_uce_adata.h5ad')
    if temp_output.exists():
        shutil.move(str(temp_output), str(output_path))

    logger.info(f"Saved embeddings to: {output_path}")
    return True


def main():
    parser = argparse.ArgumentParser(description='Process single X-Atlas chunk')
    parser.add_argument('--model', type=str, required=True, choices=['STATE', 'UCE'],
                        help='Foundation model to use')
    parser.add_argument('--cell-line', type=str, choices=['HCT116', 'HEK293T'],
                        help='Cell line to process')
    parser.add_argument('--chunk', type=int,
                        help='Chunk index to process (0-based)')
    parser.add_argument('--auto', action='store_true',
                        help='Auto-find next unprocessed chunk')
    parser.add_argument('--status', action='store_true',
                        help='Show progress status only')
    args = parser.parse_args()

    logger.info("=" * 70)
    logger.info(f"X-ATLAS CHUNK PROCESSOR - {args.model}")
    logger.info("=" * 70)

    # Show status
    if args.status:
        summary = get_progress_summary(args.model)
        logger.info("\nProgress Summary:")
        for cell_line, stats in summary.items():
            pct = 100 * stats['processed'] / stats['total'] if stats['total'] > 0 else 0
            logger.info(f"  {cell_line}: {stats['processed']}/{stats['total']} ({pct:.0f}%)")

        cell_line, chunk_idx, chunk_path = find_next_chunk(args.model)
        if cell_line:
            logger.info(f"\nNext chunk: {cell_line} chunk {chunk_idx}")
        else:
            logger.info("\nAll chunks processed!")
        return

    # Determine which chunk to process
    if args.auto:
        cell_line, chunk_idx, chunk_path = find_next_chunk(args.model)
        if cell_line is None:
            logger.info("All chunks already processed!")
            summary = get_progress_summary(args.model)
            for cl, stats in summary.items():
                logger.info(f"  {cl}: {stats['processed']}/{stats['total']}")
            return
    else:
        if not args.cell_line or args.chunk is None:
            parser.error("Either use --auto or specify --cell-line and --chunk")
        cell_line = args.cell_line
        chunk_idx = args.chunk
        chunks = get_all_chunks(cell_line)
        if chunk_idx >= len(chunks):
            logger.error(f"Chunk {chunk_idx} not found. Available: 0-{len(chunks)-1}")
            return
        chunk_path = chunks[chunk_idx]

    # Show progress
    summary = get_progress_summary(args.model)
    logger.info("\nProgress:")
    for cl, stats in summary.items():
        marker = " <--" if cl == cell_line else ""
        logger.info(f"  {cl}: {stats['processed']}/{stats['total']}{marker}")

    logger.info(f"\nProcessing: {cell_line} chunk {chunk_idx}")
    logger.info(f"Input: {chunk_path}")

    output_path = get_output_path(cell_line, chunk_idx, args.model)

    if output_path.exists():
        logger.info(f"Output already exists: {output_path}")
        logger.info("Skipping... (delete file to reprocess)")
        return

    # Process chunk
    start_time = datetime.now()

    if args.model == 'STATE':
        success = generate_state_embeddings(chunk_path, output_path)
    else:
        success = generate_uce_embeddings(chunk_path, output_path)

    elapsed = datetime.now() - start_time

    if success:
        logger.info(f"\nCompleted in {elapsed}")

        # Show updated progress
        summary = get_progress_summary(args.model)
        logger.info("\nUpdated Progress:")
        for cl, stats in summary.items():
            pct = 100 * stats['processed'] / stats['total'] if stats['total'] > 0 else 0
            logger.info(f"  {cl}: {stats['processed']}/{stats['total']} ({pct:.0f}%)")

        # Show next chunk
        next_cl, next_idx, _ = find_next_chunk(args.model)
        if next_cl:
            logger.info(f"\nNext chunk: {next_cl} chunk {next_idx}")
            logger.info("Resubmit job to continue processing.")
        else:
            logger.info("\nAll chunks processed!")
    else:
        logger.error("Processing failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
