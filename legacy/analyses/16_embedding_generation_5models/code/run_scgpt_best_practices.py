#!/usr/bin/env python3
"""
Generate scGPT Embeddings Following Official Best Practices

This script implements scGPT embedding generation following the official
documentation recommendations:
1. Proper preprocessing with Preprocessor class
2. Depth normalization (normalize_total=1e4)
3. Log1p transformation
4. HVG selection (1200 genes)
5. Value binning (51 bins)
6. Batch-aware processing

Author: Project 302006
Date: January 2026
"""

import os
import sys
import gc
import json
import time
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import scanpy as sc
import torch

warnings.filterwarnings('ignore')

# Setup logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class scGPTBestPractices:
    """Generate scGPT embeddings following official best practices."""

    def __init__(
        self,
        model_dir: str,
        n_hvg: int = 1200,
        n_bins: int = 51,
        batch_size: int = 64,
        batch_key: str = "Data_NO"
    ):
        """
        Initialize scGPT embedding generator.

        Args:
            model_dir: Path to scGPT model directory
            n_hvg: Number of highly variable genes (default: 1200)
            n_bins: Number of bins for value binning (default: 51)
            batch_size: Batch size for processing (default: 64)
            batch_key: Column name for batch information
        """
        self.model_dir = Path(model_dir)
        self.n_hvg = n_hvg
        self.n_bins = n_bins
        self.batch_size = batch_size
        self.batch_key = batch_key

        # Model files
        self.vocab_file = self.model_dir / "vocab.json"
        self.model_file = self.model_dir / "best_model.pt"
        self.args_file = self.model_dir / "args.json"

        # Verify model files exist
        for f in [self.vocab_file, self.model_file, self.args_file]:
            if not f.exists():
                raise FileNotFoundError(f"Model file not found: {f}")

        logger.info(f"Model directory: {self.model_dir}")
        logger.info(f"HVG count: {n_hvg}")
        logger.info(f"Bins: {n_bins}")
        logger.info(f"Batch key: {batch_key}")

    def preprocess_data(self, adata: sc.AnnData) -> sc.AnnData:
        """
        Preprocess data following scGPT best practices.

        Steps:
        1. Filter genes by counts
        2. Depth normalization (1e4)
        3. Log1p transformation
        4. HVG selection
        5. Value binning
        """
        logger.info("=" * 60)
        logger.info("PREPROCESSING (Best Practices)")
        logger.info("=" * 60)

        adata = adata.copy()

        # Store raw counts
        adata.layers["counts"] = adata.X.copy()

        # Step 1: Filter genes with very low counts
        logger.info("Step 1: Filtering genes...")
        sc.pp.filter_genes(adata, min_counts=3)
        logger.info(f"  After gene filter: {adata.shape}")

        # Step 2: Depth normalization
        logger.info("Step 2: Depth normalization (target=1e4)...")
        sc.pp.normalize_total(adata, target_sum=1e4)

        # Step 3: Log1p transformation
        logger.info("Step 3: Log1p transformation...")
        sc.pp.log1p(adata)
        adata.layers["X_normed"] = adata.X.copy()

        # Step 4: HVG selection (batch-aware)
        logger.info(f"Step 4: Selecting {self.n_hvg} highly variable genes...")

        # Try seurat_v3 first, fall back to cell_ranger if it fails
        hvg_success = False
        for flavor in ["seurat_v3", "cell_ranger"]:
            try:
                if self.batch_key and self.batch_key in adata.obs.columns:
                    logger.info(f"  Trying {flavor} with batch_key={self.batch_key}...")
                    if flavor == "seurat_v3":
                        sc.pp.highly_variable_genes(
                            adata,
                            n_top_genes=self.n_hvg,
                            flavor=flavor,
                            layer="counts",
                            batch_key=self.batch_key
                        )
                    else:
                        # cell_ranger doesn't support layer parameter
                        sc.pp.highly_variable_genes(
                            adata,
                            n_top_genes=self.n_hvg,
                            flavor=flavor,
                            batch_key=self.batch_key
                        )
                else:
                    logger.info(f"  Trying {flavor} HVG selection...")
                    if flavor == "seurat_v3":
                        sc.pp.highly_variable_genes(
                            adata,
                            n_top_genes=self.n_hvg,
                            flavor=flavor,
                            layer="counts"
                        )
                    else:
                        sc.pp.highly_variable_genes(
                            adata,
                            n_top_genes=self.n_hvg,
                            flavor=flavor
                        )
                hvg_success = True
                logger.info(f"  HVG selection successful with {flavor}")
                break
            except Exception as e:
                logger.warning(f"  {flavor} failed: {e}")
                continue

        if not hvg_success:
            # Last resort: select top variable genes by variance
            logger.warning("  All HVG methods failed, using variance-based selection...")
            gene_vars = np.var(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X, axis=0)
            top_genes = np.argsort(gene_vars)[-self.n_hvg:]
            adata.var['highly_variable'] = False
            adata.var.iloc[top_genes, adata.var.columns.get_loc('highly_variable')] = True

        # Subset to HVGs
        adata = adata[:, adata.var.highly_variable].copy()
        logger.info(f"  After HVG selection: {adata.shape}")

        # Step 5: Value binning
        logger.info(f"Step 5: Value binning ({self.n_bins} bins)...")
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()

        X_binned = np.zeros(X.shape, dtype=np.int32)
        for i in range(X.shape[1]):
            gene_expr = X[:, i]
            non_zero_mask = gene_expr > 0
            if non_zero_mask.sum() > self.n_bins:
                # Use quantile binning for genes with enough non-zero values
                non_zero_vals = gene_expr[non_zero_mask]
                bins = np.quantile(non_zero_vals, np.linspace(0, 1, self.n_bins))
                bins = np.unique(bins)
                X_binned[:, i] = np.digitize(gene_expr, bins)
            elif non_zero_mask.sum() > 0:
                # Simple scaling for sparse genes
                max_val = gene_expr.max()
                if max_val > 0:
                    X_binned[:, i] = np.clip(
                        (gene_expr / max_val * (self.n_bins - 1)).astype(int),
                        0, self.n_bins - 1
                    )

        adata.layers["X_binned"] = X_binned
        logger.info(f"  Binned expression range: [{X_binned.min()}, {X_binned.max()}]")

        logger.info("Preprocessing complete!")
        return adata

    def load_model(self):
        """Load scGPT model and vocabulary."""
        logger.info("=" * 60)
        logger.info("LOADING MODEL")
        logger.info("=" * 60)

        # Import scGPT modules
        from scgpt.model import TransformerModel
        from scgpt.tokenizer.gene_tokenizer import GeneVocab

        # Load vocabulary
        logger.info(f"Loading vocabulary from {self.vocab_file}...")
        self.vocab = GeneVocab.from_file(str(self.vocab_file))

        # Add special tokens
        special_tokens = ["<pad>", "<cls>", "<eoc>"]
        for s in special_tokens:
            if s not in self.vocab:
                self.vocab.append_token(s)
        self.vocab.set_default_index(self.vocab["<pad>"])
        logger.info(f"  Vocabulary size: {len(self.vocab)}")

        # Load model args
        logger.info(f"Loading model config from {self.args_file}...")
        with open(self.args_file, "r") as f:
            args = json.load(f)

        # Create model
        logger.info("Creating TransformerModel...")
        self.model = TransformerModel(
            ntoken=len(self.vocab),
            d_model=args.get("embsize", 512),
            nhead=args.get("nheads", 8),
            d_hid=args.get("d_hid", 512),
            nlayers=args.get("nlayers", 12),
            vocab=self.vocab,
            dropout=0.0,
            pad_token="<pad>",
            pad_value=0
        )

        # Load weights
        logger.info(f"Loading model weights from {self.model_file}...")
        state_dict = torch.load(str(self.model_file), map_location='cpu')
        self.model.load_state_dict(state_dict, strict=False)
        self.model.eval()

        # Move to GPU if available
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = self.model.to(self.device)
        logger.info(f"  Model loaded on: {self.device}")

        return self.model

    def generate_embeddings(self, adata: sc.AnnData) -> np.ndarray:
        """
        Generate embeddings for preprocessed data.

        Args:
            adata: Preprocessed AnnData with X_binned layer

        Returns:
            Embedding matrix (n_cells x embedding_dim)
        """
        logger.info("=" * 60)
        logger.info("GENERATING EMBEDDINGS")
        logger.info("=" * 60)

        # Get binned expression
        X_binned = adata.layers["X_binned"]
        if hasattr(X_binned, 'toarray'):
            X_binned = X_binned.toarray()

        # Map genes to vocabulary
        gene_names = adata.var_names.tolist()
        gene_ids = []
        valid_indices = []

        for idx, gene in enumerate(gene_names):
            if gene in self.vocab:
                gene_ids.append(self.vocab[gene])
                valid_indices.append(idx)

        logger.info(f"Genes in vocabulary: {len(gene_ids)}/{len(gene_names)}")

        if len(gene_ids) == 0:
            raise ValueError("No genes found in vocabulary!")

        # Filter to vocab genes
        X_binned = X_binned[:, valid_indices]

        # Process in batches
        n_cells = X_binned.shape[0]
        embeddings = []

        logger.info(f"Processing {n_cells:,} cells in batches of {self.batch_size}...")
        start_time = time.time()

        with torch.no_grad():
            for batch_start in range(0, n_cells, self.batch_size):
                batch_end = min(batch_start + self.batch_size, n_cells)
                batch_X = X_binned[batch_start:batch_end]
                batch_size_actual = batch_end - batch_start

                # Create tensors
                gene_tensor = torch.tensor(
                    [gene_ids] * batch_size_actual,
                    dtype=torch.long
                ).to(self.device)

                value_tensor = torch.tensor(
                    batch_X,
                    dtype=torch.float
                ).to(self.device)

                padding_mask = torch.zeros(
                    batch_size_actual, len(gene_ids),
                    dtype=torch.bool
                ).to(self.device)

                # Forward pass
                output = self.model(
                    gene_tensor,
                    value_tensor,
                    src_key_padding_mask=padding_mask
                )

                # Extract embeddings (mean pooling)
                if hasattr(output, 'cell_emb') and output.cell_emb is not None:
                    batch_emb = output.cell_emb.cpu().numpy()
                elif isinstance(output, dict) and 'cell_emb' in output:
                    batch_emb = output['cell_emb'].cpu().numpy()
                else:
                    batch_emb = output.mean(dim=1).cpu().numpy()

                embeddings.append(batch_emb)

                # Progress
                if (batch_start // self.batch_size) % 100 == 0:
                    elapsed = time.time() - start_time
                    progress = batch_end / n_cells * 100
                    rate = batch_end / elapsed if elapsed > 0 else 0
                    eta = (n_cells - batch_end) / rate if rate > 0 else 0
                    logger.info(
                        f"  Progress: {batch_end:,}/{n_cells:,} ({progress:.1f}%) "
                        f"| Rate: {rate:.0f} cells/s | ETA: {eta/60:.1f} min"
                    )

                # Cleanup
                del gene_tensor, value_tensor, padding_mask, output
                if batch_start % (self.batch_size * 10) == 0:
                    gc.collect()
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()

        # Combine embeddings
        embeddings = np.vstack(embeddings)
        elapsed = time.time() - start_time

        logger.info(f"Embedding generation complete!")
        logger.info(f"  Shape: {embeddings.shape}")
        logger.info(f"  Time: {elapsed/60:.1f} minutes")
        logger.info(f"  Rate: {n_cells/elapsed:.0f} cells/second")

        return embeddings

    def process(
        self,
        input_path: str,
        output_path: str
    ) -> sc.AnnData:
        """
        Full pipeline: preprocess, load model, generate embeddings.

        Args:
            input_path: Path to input h5ad file
            output_path: Path to save output h5ad file

        Returns:
            AnnData with embeddings in obsm['X_scgpt_bp']
        """
        logger.info("=" * 60)
        logger.info("scGPT EMBEDDING GENERATION (Best Practices)")
        logger.info("=" * 60)
        logger.info(f"Input: {input_path}")
        logger.info(f"Output: {output_path}")
        logger.info(f"Start time: {datetime.now()}")

        # Load data
        logger.info("\nLoading input data...")
        adata_raw = sc.read_h5ad(input_path)
        logger.info(f"  Shape: {adata_raw.shape}")

        # Preprocess
        adata = self.preprocess_data(adata_raw)

        # Load model
        self.load_model()

        # Generate embeddings
        embeddings = self.generate_embeddings(adata)

        # Add embeddings to original data
        logger.info("\nSaving results...")
        adata_raw.obsm['X_scgpt_bp'] = embeddings  # bp = best practices

        # Also store preprocessing info
        adata_raw.uns['scgpt_best_practices'] = {
            'n_hvg': self.n_hvg,
            'n_bins': self.n_bins,
            'batch_key': self.batch_key,
            'embedding_dim': embeddings.shape[1],
            'genes_used': adata.shape[1],
            'date': str(datetime.now())
        }

        # Save
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        adata_raw.write_h5ad(output_path)
        logger.info(f"  Saved to: {output_path}")

        logger.info("\n" + "=" * 60)
        logger.info("COMPLETE")
        logger.info("=" * 60)

        return adata_raw


def main():
    """Main execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate scGPT embeddings following best practices'
    )
    parser.add_argument(
        '--input', type=str, required=True,
        help='Path to input h5ad file'
    )
    parser.add_argument(
        '--output', type=str, required=True,
        help='Path to output h5ad file'
    )
    parser.add_argument(
        '--model-dir', type=str,
        default="/bigdata/godziklab/shared/Xinru/302006/models/scGPT_human/scGPT_human",
        help='Path to scGPT model directory'
    )
    parser.add_argument(
        '--n-hvg', type=int, default=1200,
        help='Number of highly variable genes (default: 1200)'
    )
    parser.add_argument(
        '--n-bins', type=int, default=51,
        help='Number of bins for value binning (default: 51)'
    )
    parser.add_argument(
        '--batch-size', type=int, default=64,
        help='Batch size for processing (default: 64)'
    )
    parser.add_argument(
        '--batch-key', type=str, default="Data_NO",
        help='Column name for batch information (default: Data_NO)'
    )

    args = parser.parse_args()

    # Check CUDA
    logger.info(f"PyTorch version: {torch.__version__}")
    logger.info(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        logger.info(f"CUDA device: {torch.cuda.get_device_name(0)}")

    # Run pipeline
    processor = scGPTBestPractices(
        model_dir=args.model_dir,
        n_hvg=args.n_hvg,
        n_bins=args.n_bins,
        batch_size=args.batch_size,
        batch_key=args.batch_key
    )

    processor.process(args.input, args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
