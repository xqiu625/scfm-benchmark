#!/usr/bin/env python3
"""
scGPT Embedding Generation for Single-Cell Data Chunks

This script processes single-cell RNA-seq data chunks and generates scGPT embeddings
using batched inference for improved performance.

Usage:
    python run_scgpt_chunk.py --input data.h5ad --output embeddings.h5ad --model-dir /path/to/model

Author: Project 302006
Date: January 2026
"""

import argparse
import gc
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
import torch

warnings.filterwarnings('ignore')


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate scGPT embeddings for a data chunk')
    parser.add_argument('--input', '-i', required=True, help='Input h5ad file')
    parser.add_argument('--output', '-o', required=True, help='Output h5ad file')
    parser.add_argument('--model-dir', '-m', required=True, help='Path to scGPT model directory')
    parser.add_argument('--batch-size', '-b', type=int, default=64,
                        help='Batch size for inference (default: 64)')
    parser.add_argument('--n-bins', type=int, default=51,
                        help='Number of bins for expression discretization (default: 51)')
    return parser.parse_args()


def load_scgpt_components(model_dir: Path):
    """Load scGPT vocabulary and model."""
    from scgpt.model import TransformerModel
    from scgpt.tokenizer.gene_tokenizer import GeneVocab

    vocab_file = model_dir / "vocab.json"
    model_file = model_dir / "best_model.pt"
    args_file = model_dir / "args.json"

    # Validate files exist
    for f in [vocab_file, model_file, args_file]:
        if not f.exists():
            raise FileNotFoundError(f"Required file not found: {f}")

    # Load vocabulary
    print(f"Loading vocabulary from {vocab_file}")
    vocab = GeneVocab.from_file(vocab_file)
    special_tokens = ["<pad>", "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])

    # Load model configuration
    with open(args_file, "r") as f:
        model_args = json.load(f)

    # Initialize model
    print(f"Loading model from {model_file}")
    model = TransformerModel(
        len(vocab),
        model_args.get("embsize", 512),
        model_args.get("nheads", 8),
        model_args.get("d_hid", 512),
        model_args.get("nlayers", 12),
        vocab=vocab,
        dropout=0.0,
        pad_token="<pad>",
        pad_value=0
    )

    state_dict = torch.load(model_file, map_location='cpu')
    model.load_state_dict(state_dict, strict=False)
    model.eval()

    return vocab, model


def filter_genes_to_vocabulary(adata: sc.AnnData, vocab) -> tuple:
    """Filter genes to those present in the scGPT vocabulary.

    Returns:
        tuple: (gene_ids, valid_indices, filtered_adata)
    """
    gene_names = adata.var.index.tolist()
    gene_ids = []
    valid_indices = []

    for idx, gene in enumerate(gene_names):
        if gene in vocab:
            gene_ids.append(vocab[gene])
            valid_indices.append(idx)

    print(f"Genes in vocabulary: {len(gene_ids):,} / {len(gene_names):,} ({100*len(gene_ids)/len(gene_names):.1f}%)")

    # Filter AnnData to vocabulary genes
    adata_filtered = adata[:, valid_indices].copy()

    return gene_ids, valid_indices, adata_filtered


def bin_expression_values(X: np.ndarray, n_bins: int = 51) -> np.ndarray:
    """Discretize expression values into bins using quantile binning.

    Args:
        X: Expression matrix (cells x genes)
        n_bins: Number of bins (default: 51, giving values 0-50)

    Returns:
        np.ndarray: Binned expression matrix
    """
    n_cells, n_genes = X.shape
    X_binned = np.zeros_like(X, dtype=np.int32)

    print(f"Binning {n_genes:,} genes into {n_bins} bins...")

    for i in range(n_genes):
        if i % 2000 == 0 and i > 0:
            print(f"  Binning progress: {i:,}/{n_genes:,} genes")

        gene_expr = X[:, i]
        max_expr = np.max(gene_expr)

        if max_expr > 0:
            non_zero = gene_expr[gene_expr > 0]
            if len(non_zero) > n_bins:
                # Quantile binning for genes with enough non-zero values
                bins = np.quantile(non_zero, np.linspace(0, 1, n_bins))
                bins = np.unique(bins)
                X_binned[:, i] = np.digitize(gene_expr, bins)
            else:
                # Linear scaling for sparse genes
                X_binned[:, i] = np.clip(
                    (gene_expr / max_expr * (n_bins - 1)).astype(int),
                    0,
                    n_bins - 1
                )

    return X_binned


def generate_embeddings_batched(
    model,
    gene_ids: list,
    X_binned: np.ndarray,
    batch_size: int = 64
) -> np.ndarray:
    """Generate embeddings using batched inference.

    Args:
        model: Loaded scGPT model
        gene_ids: List of gene vocabulary IDs
        X_binned: Binned expression matrix (cells x genes)
        batch_size: Number of cells per batch

    Returns:
        np.ndarray: Cell embeddings matrix
    """
    n_cells = X_binned.shape[0]
    n_batches = (n_cells + batch_size - 1) // batch_size
    embeddings = []

    print(f"Generating embeddings for {n_cells:,} cells in {n_batches} batches (batch_size={batch_size})")

    start_time = time.time()
    gene_tensor = torch.tensor(gene_ids, dtype=torch.long)

    with torch.no_grad():
        for batch_idx in range(n_batches):
            start_idx = batch_idx * batch_size
            end_idx = min(start_idx + batch_size, n_cells)
            batch_cells = end_idx - start_idx

            # Progress reporting
            if batch_idx % 10 == 0 or batch_idx == n_batches - 1:
                elapsed = time.time() - start_time
                cells_done = end_idx
                rate = cells_done / max(elapsed, 0.001)
                remaining = (n_cells - cells_done) / max(rate, 1)
                print(f"  Batch {batch_idx+1}/{n_batches} | "
                      f"Cells: {cells_done:,}/{n_cells:,} | "
                      f"Rate: {rate:.0f} cells/s | "
                      f"ETA: {remaining/60:.1f} min")

            try:
                # Prepare batch tensors
                # Gene IDs repeated for each cell in batch
                gene_batch = gene_tensor.unsqueeze(0).expand(batch_cells, -1)

                # Expression values for this batch
                value_batch = torch.tensor(
                    X_binned[start_idx:end_idx],
                    dtype=torch.float
                )

                # Padding mask (no padding needed)
                padding_mask = torch.zeros(
                    batch_cells,
                    len(gene_ids),
                    dtype=torch.bool
                )

                # Forward pass
                output = model(gene_batch, value_batch, src_key_padding_mask=padding_mask)

                # Extract embeddings
                if hasattr(output, 'cell_emb') and output.cell_emb is not None:
                    batch_emb = output.cell_emb.cpu().numpy()
                elif isinstance(output, dict) and 'cell_emb' in output:
                    batch_emb = output['cell_emb'].cpu().numpy()
                else:
                    # Mean pooling over sequence dimension
                    batch_emb = output.mean(dim=1).cpu().numpy()

                embeddings.append(batch_emb)

                # Cleanup
                del gene_batch, value_batch, padding_mask, output

            except Exception as e:
                print(f"ERROR in batch {batch_idx}: {e}")
                raise

            # Periodic memory cleanup
            if batch_idx % 50 == 0:
                gc.collect()

    # Combine all batch embeddings
    final_embeddings = np.vstack(embeddings)

    elapsed_total = time.time() - start_time
    print(f"Embedding generation complete: {final_embeddings.shape} in {elapsed_total/60:.1f} min")
    print(f"Average rate: {n_cells/elapsed_total:.0f} cells/s")

    return final_embeddings


def main():
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    model_dir = Path(args.model_dir)

    print("=" * 60)
    print("scGPT EMBEDDING GENERATION")
    print("=" * 60)
    print(f"Input:      {input_path}")
    print(f"Output:     {output_path}")
    print(f"Model:      {model_dir}")
    print(f"Batch size: {args.batch_size}")
    print("=" * 60)

    # Check if output already exists
    if output_path.exists():
        print(f"Output already exists: {output_path}")
        print("Skipping processing.")
        return 0

    # Validate input
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}")
        return 1

    start_time = time.time()

    # Load data
    print("\n[1/5] Loading input data...")
    adata = sc.read_h5ad(input_path)
    print(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Load model components
    print("\n[2/5] Loading scGPT model...")
    vocab, model = load_scgpt_components(model_dir)

    # Filter genes
    print("\n[3/5] Filtering genes to vocabulary...")
    gene_ids, valid_indices, adata_filtered = filter_genes_to_vocabulary(adata, vocab)

    # Get expression matrix
    X = adata_filtered.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    # Bin expression values
    print("\n[4/5] Binning expression values...")
    X_binned = bin_expression_values(X, n_bins=args.n_bins)

    # Generate embeddings
    print("\n[5/5] Generating embeddings...")
    embeddings = generate_embeddings_batched(
        model,
        gene_ids,
        X_binned,
        batch_size=args.batch_size
    )

    # Create output AnnData
    print("\nSaving results...")
    result_adata = sc.AnnData(
        X=embeddings,
        obs=adata.obs.copy()
    )
    result_adata.obsm['X_scgpt'] = embeddings

    # Store metadata
    total_time = time.time() - start_time
    result_adata.uns['scgpt_processing'] = {
        'model': 'scGPT_human',
        'embedding_dim': embeddings.shape[1],
        'genes_used': len(gene_ids),
        'genes_total': adata.shape[1],
        'processing_time_seconds': total_time,
        'batch_size': args.batch_size,
        'n_bins': args.n_bins,
        'input_file': str(input_path),
    }

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_adata.write_h5ad(output_path)

    print(f"\nSaved: {output_path}")
    print(f"Embedding shape: {embeddings.shape}")
    print(f"Total time: {total_time/60:.1f} min ({total_time/3600:.2f} h)")

    # Verify output
    print("\nVerifying output...")
    result = sc.read_h5ad(output_path)
    print(f"Verified: {result.shape}")
    print(f"Embedding range: [{np.min(result.X):.4f}, {np.max(result.X):.4f}]")

    # Cleanup
    del model, adata, embeddings, result_adata
    gc.collect()

    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
