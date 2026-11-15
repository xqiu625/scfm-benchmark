#!/usr/bin/env python3
"""
Prepare Data for Geneformer Embedding Generation

Purpose: Convert COVID/HEK293T/HCT116 datasets to Geneformer-compatible format
Requirements:
  - Raw counts (no normalization)
  - Ensembl gene IDs in adata.var['ensembl_id']
  - n_counts in adata.obs
  - H5AD format

Uses the TranscriptFormer-prepared data which already has Ensembl IDs.

Author: Project 302006
Date: December 2025
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from datetime import datetime
import sys

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'prepare_geneformer_data_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class GeneformerDataPrep:
    """Prepare data for Geneformer processing"""

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_dir = Path(base_dir)
        # Use TranscriptFormer-ready data which has Ensembl IDs
        self.input_dir = self.base_dir / "01_DATA/transcriptformer_ready"
        self.output_dir = self.base_dir / "01_DATA/geneformer_ready"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def prepare_dataset(self, dataset_name):
        """Prepare a single dataset for Geneformer"""
        logger.info(f"\n{'='*80}")
        logger.info(f"Preparing for Geneformer: {dataset_name}")
        logger.info(f"{'='*80}")

        # Input file (from TranscriptFormer prep - has Ensembl IDs)
        input_file = self.input_dir / f"{dataset_name}_transcriptformer_ready.h5ad"

        if not input_file.exists():
            logger.error(f"Input file not found: {input_file}")
            logger.error("Please run TranscriptFormer data preparation first")
            return None

        # Load data
        logger.info(f"Loading {input_file.name}...")
        adata = sc.read_h5ad(input_file)
        logger.info(f"  Shape: {adata.shape}")
        logger.info(f"  var columns: {list(adata.var.columns)}")
        logger.info(f"  obs columns: {list(adata.obs.columns)}")

        # Geneformer requires ensembl_id in var
        if 'ensembl_id' not in adata.var.columns:
            logger.error("Missing 'ensembl_id' column in adata.var!")
            return None

        # Check if var_names are Ensembl IDs or gene symbols
        first_var = str(adata.var_names[0])
        if first_var.startswith('ENSG'):
            logger.info("  var_names are already Ensembl IDs ✓")
        else:
            # Set var_names to ensembl_id for Geneformer
            logger.info("  Setting var_names to Ensembl IDs...")
            if 'gene_symbol' not in adata.var.columns:
                adata.var['gene_symbol'] = adata.var_names.values
            adata.var_names = pd.Index(adata.var['ensembl_id'].values)

        # Ensure n_counts is present
        logger.info("\nChecking required columns...")
        if 'n_counts' not in adata.obs.columns:
            if 'total_counts' in adata.obs.columns:
                adata.obs['n_counts'] = adata.obs['total_counts']
                logger.info("  Copied 'total_counts' to 'n_counts'")
            else:
                # Calculate total counts per cell
                logger.info("  Calculating n_counts from expression matrix...")
                if hasattr(adata.X, 'sum'):
                    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
                else:
                    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()

        logger.info(f"  n_counts range: {adata.obs['n_counts'].min():.0f} - {adata.obs['n_counts'].max():.0f}")

        # Add filter_pass column (Geneformer uses this)
        if 'filter_pass' not in adata.obs.columns:
            adata.obs['filter_pass'] = True
            logger.info("  Added 'filter_pass' column (all True)")

        # Add cell_type if not present (for Geneformer metadata)
        if 'cell_type' not in adata.obs.columns:
            if dataset_name == 'covid':
                # Use severity as cell_type proxy
                if 'severity' in adata.obs.columns:
                    adata.obs['cell_type'] = adata.obs['severity'].astype(str)
                elif 'condition' in adata.obs.columns:
                    adata.obs['cell_type'] = adata.obs['condition'].astype(str)
                else:
                    adata.obs['cell_type'] = 'platelet'
            else:
                adata.obs['cell_type'] = dataset_name.upper()
            logger.info(f"  Added 'cell_type' column: {adata.obs['cell_type'].unique()[:5].tolist()}")

        # Remove adata.raw if present (can cause issues)
        if adata.raw is not None:
            logger.info("  Removing adata.raw to avoid conflicts")
            adata.raw = None

        # Convert categorical columns to strings for serialization
        for col in adata.obs.columns:
            if adata.obs[col].dtype.name == 'category':
                adata.obs[col] = adata.obs[col].astype(str)

        for col in adata.var.columns:
            if adata.var[col].dtype.name == 'category':
                adata.var[col] = adata.var[col].astype(str)

        # Verify data is ready
        logger.info("\nVerification:")
        logger.info(f"  Cells: {adata.n_obs}")
        logger.info(f"  Genes: {adata.n_vars}")
        logger.info(f"  var_names (first 5): {adata.var_names[:5].tolist()}")
        logger.info(f"  ensembl_id (first 5): {adata.var['ensembl_id'][:5].tolist()}")
        logger.info(f"  n_counts present: {'n_counts' in adata.obs.columns}")

        # Save prepared data
        output_file = self.output_dir / f"{dataset_name}_geneformer_ready.h5ad"
        logger.info(f"\nSaving to: {output_file}")
        adata.write_h5ad(output_file)

        logger.info(f"✓ {dataset_name} prepared for Geneformer")

        return adata

    def prepare_all_datasets(self):
        """Prepare all three datasets"""
        logger.info("=" * 80)
        logger.info("GENEFORMER DATA PREPARATION")
        logger.info("=" * 80)

        datasets = ['covid', 'hek293t', 'hct116']
        results = {}

        for dataset_name in datasets:
            try:
                adata = self.prepare_dataset(dataset_name)
                if adata is not None:
                    results[dataset_name] = adata
            except Exception as e:
                logger.error(f"Failed to prepare {dataset_name}: {str(e)}")
                import traceback
                traceback.print_exc()
                continue

        # Summary
        logger.info("\n" + "=" * 80)
        logger.info("PREPARATION COMPLETE")
        logger.info("=" * 80)
        logger.info(f"\nOutput directory: {self.output_dir}")
        logger.info("\nPrepared datasets:")
        for dataset_name, adata in results.items():
            logger.info(f"  {dataset_name}: {adata.shape[0]} cells x {adata.shape[1]} genes")

        logger.info("\nNext steps:")
        logger.info("  1. Run Geneformer tokenization")
        logger.info("  2. Generate Geneformer embeddings")

        return results


def main():
    """Main execution"""
    try:
        prep = GeneformerDataPrep()
        results = prep.prepare_all_datasets()
        return 0

    except Exception as e:
        logger.error(f"\nERROR: {str(e)}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
