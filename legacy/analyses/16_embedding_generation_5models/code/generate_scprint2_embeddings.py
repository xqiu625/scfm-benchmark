#!/usr/bin/env python3
"""
Generate scPRINT-2 Embeddings for COVID/Mouse datasets

Purpose: Generate cell embeddings using scPRINT-2 foundation model.
Key features:
  - Uses scPRINT-2 Embedder class for zero-shot cell embeddings
  - Supports human (NCBITaxon:9606) and mouse (NCBITaxon:10090)
  - Raw counts input (no normalization/integration needed)
  - Gene names as Ensembl IDs or HUGO symbols
  - Outputs embeddings in obsm['X_scprint2'] for benchmark compatibility

Requirements:
  - scprint2 package: uv pip install scprint2
  - lamin initialized: lamin init --storage <path> --name <name> --modules bionty
  - Model checkpoint: hf download jkobject/scPRINT small-v2.ckpt --local-dir .

Author: Project 302006
Date: April 2026
"""

import scanpy as sc
import numpy as np
import torch
from pathlib import Path
import logging
from datetime import datetime
import sys
import argparse

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(
            f'scprint2_embeddings_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
        ),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


# Species ontology IDs
SPECIES_MAP = {
    'human': 'NCBITaxon:9606',
    'mouse': 'NCBITaxon:10090',
    'NCBITaxon:9606': 'NCBITaxon:9606',
    'NCBITaxon:10090': 'NCBITaxon:10090',
}


def setup_lamin():
    """Initialize lamin if not already set up."""
    try:
        import lamindb as ln
        ln.connect("scprint2_db")
        logger.info("Connected to existing lamin instance")
    except Exception:
        logger.warning("Could not connect to lamin. If this is first run, "
                       "run: lamin init --storage <path> --name scprint2_db --modules bionty")
        logger.warning("Then: scprint2 easy_setup")
        raise


def load_model(ckpt_path):
    """Load scPRINT-2 model from checkpoint."""
    from scprint2 import scPRINT2
    from scdataloader.utils import load_genes

    logger.info(f"Loading scPRINT-2 from: {ckpt_path}")

    model = scPRINT2.load_from_checkpoint(
        str(ckpt_path),
        precpt_gene_emb=None,
        gene_pos_file=None,
    )

    # Handle CPU-only environments
    if not torch.cuda.is_available():
        logger.warning("No CUDA GPU detected, using CPU (will be slow)")
        model = model.to(torch.float32)
        device = "cpu"
    else:
        device = "cuda"
        logger.info(f"Using GPU: {torch.cuda.get_device_name(0)}")

    model = model.to(device)

    # Remove genes that no longer exist in lamin ontology
    missing = set(model.genes) - set(load_genes(model.organisms).index)
    if len(missing) > 0:
        logger.warning(
            f"{len(missing)} genes mismatch between model and ontology, removing..."
        )
        model._rm_genes(missing)

    logger.info(f"Model loaded: {len(model.genes)} genes, device={device}")
    return model


def prepare_adata(adata, species_id):
    """Prepare AnnData for scPRINT-2 inference.

    scPRINT-2 requires:
      - Raw counts (not normalized/integrated)
      - obs['organism_ontology_term_id'] with species ontology ID
      - var_names as Ensembl IDs or HUGO gene symbols
    """
    logger.info(f"Preparing AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Add species annotation if missing
    if 'organism_ontology_term_id' not in adata.obs.columns:
        logger.info(f"Adding organism_ontology_term_id: {species_id}")
        adata.obs['organism_ontology_term_id'] = species_id
    else:
        existing = adata.obs['organism_ontology_term_id'].unique()
        logger.info(f"Existing organism IDs: {existing}")

    # Ensure raw counts (warn if data looks normalized)
    if adata.X is not None:
        x_sample = adata.X[:100] if hasattr(adata.X, 'toarray') else adata.X[:100]
        if hasattr(x_sample, 'toarray'):
            x_sample = x_sample.toarray()
        max_val = np.max(x_sample)
        if max_val < 20:
            logger.warning(
                f"Max expression value is {max_val:.2f} — data may be normalized. "
                "scPRINT-2 expects raw counts!"
            )
        else:
            logger.info(f"Data looks like raw counts (max={max_val:.0f})")

    # Check gene naming
    sample_genes = list(adata.var_names[:5])
    if any(g.startswith('ENSG') or g.startswith('ENSMUSG') for g in sample_genes):
        logger.info("Gene names appear to be Ensembl IDs")
    else:
        logger.info(f"Gene names appear to be symbols (e.g., {sample_genes[:3]})")

    logger.info(f"AnnData prepared: {adata.shape}")
    return adata


def generate_embeddings(model, adata, batch_size=64, max_len=2000, how='random expr'):
    """Generate cell embeddings using scPRINT-2 Embedder."""
    from scprint2.tasks import Embedder

    logger.info(f"Generating embeddings (batch_size={batch_size}, "
                f"max_len={max_len}, how='{how}')")

    embedder = Embedder(
        batch_size=batch_size,
        num_workers=8,
        how=how,
        max_len=max_len,
        doclass=False,   # skip classification — we only need embeddings
        doplot=False,     # skip UMAP/plots
        pred_embedding=['cell_type_ontology_term_id'],
        use_knn=False,
    )

    # Run embedding
    logger.info("Running scPRINT-2 inference...")
    start = datetime.now()
    embedder(model, adata=adata)
    elapsed = (datetime.now() - start).total_seconds()
    logger.info(f"Inference completed in {elapsed:.1f}s "
                f"({adata.shape[0] / elapsed:.0f} cells/sec)")

    # Check which obsm key was created
    scprint_keys = [k for k in adata.obsm.keys() if 'scprint' in k.lower()]
    logger.info(f"scPRINT-2 obsm keys created: {scprint_keys}")

    # Rename to X_scprint2 for benchmark compatibility
    if 'scprint_emb' in adata.obsm:
        emb = adata.obsm['scprint_emb']
        adata.obsm['X_scprint2'] = emb
        logger.info(f"Embeddings shape: {emb.shape} -> obsm['X_scprint2']")
    else:
        # Try other possible keys
        for key in scprint_keys:
            if 'umap' not in key:
                emb = adata.obsm[key]
                adata.obsm['X_scprint2'] = emb
                logger.info(f"Embeddings from '{key}': {emb.shape} -> obsm['X_scprint2']")
                break
        else:
            logger.error("No scPRINT-2 embeddings found in obsm!")
            logger.error(f"Available keys: {list(adata.obsm.keys())}")
            raise RuntimeError("scPRINT-2 embedding generation failed")

    return adata


def main():
    parser = argparse.ArgumentParser(
        description='Generate scPRINT-2 cell embeddings'
    )
    parser.add_argument('--input', type=str, required=True,
                        help='Input h5ad file (raw counts)')
    parser.add_argument('--output', type=str, required=True,
                        help='Output h5ad file with embeddings')
    parser.add_argument('--ckpt', type=str, required=True,
                        help='Path to scPRINT-2 checkpoint (e.g., small-v2.ckpt)')
    parser.add_argument('--species', type=str, default='human',
                        help='Species: "human", "mouse", or ontology ID '
                             '(default: human)')
    parser.add_argument('--batch-size', type=int, default=64,
                        help='Batch size for inference (default: 64)')
    parser.add_argument('--max-len', type=int, default=2000,
                        help='Max gene sequence length (default: 2000)')
    parser.add_argument('--how', type=str, default='random expr',
                        choices=['random expr', 'most var', 'some', 'most expr'],
                        help='Gene selection method (default: "random expr")')
    args = parser.parse_args()

    # Resolve species
    species_id = SPECIES_MAP.get(args.species, args.species)
    if not species_id.startswith('NCBITaxon:'):
        logger.error(f"Unknown species: {args.species}. "
                     "Use 'human', 'mouse', or NCBITaxon ID.")
        return 1

    logger.info("=" * 80)
    logger.info("scPRINT-2 EMBEDDING GENERATION")
    logger.info("=" * 80)
    logger.info(f"Input:   {args.input}")
    logger.info(f"Output:  {args.output}")
    logger.info(f"Ckpt:    {args.ckpt}")
    logger.info(f"Species: {species_id}")
    logger.info(f"Batch:   {args.batch_size}")

    # Validate paths
    input_path = Path(args.input)
    ckpt_path = Path(args.ckpt)
    output_path = Path(args.output)

    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return 1
    if not ckpt_path.exists():
        logger.error(f"Checkpoint not found: {ckpt_path}")
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        # Step 1: Setup lamin
        setup_lamin()

        # Step 2: Load model
        model = load_model(ckpt_path)

        # Step 3: Load and prepare data
        logger.info(f"\nLoading input: {input_path}")
        adata = sc.read_h5ad(str(input_path))
        logger.info(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")
        adata = prepare_adata(adata, species_id)

        # Step 4: Generate embeddings
        adata = generate_embeddings(
            model, adata,
            batch_size=args.batch_size,
            max_len=args.max_len,
            how=args.how
        )

        # Step 5: Save
        logger.info(f"\nSaving to: {output_path}")
        adata.write_h5ad(str(output_path))

        emb_shape = adata.obsm['X_scprint2'].shape
        logger.info("=" * 80)
        logger.info("COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Cells:      {emb_shape[0]:,}")
        logger.info(f"Dimensions: {emb_shape[1]}")
        logger.info(f"Output:     {output_path}")
        logger.info(f"Obsm key:   X_scprint2")

        return 0

    except Exception as e:
        logger.error(f"\nFATAL ERROR: {str(e)}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
