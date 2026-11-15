#!/usr/bin/env python3
"""
Clean scGPT Embeddings: Remove Orion Data (HCT116, HEK293T)

The scGPT embedding file contains mixed data:
- Platelet COVID data (~47,000 cells) - KEEP
- HCT116 cell line (~90,000 cells) - REMOVE
- HEK293T cell line (~90,000 cells) - REMOVE

This script filters to keep only platelet data.

Author: Project 302006
Date: January 2026
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def clean_scgpt_embeddings(input_path, output_path):
    """Remove Orion data (HCT116, HEK293T) from scGPT embeddings"""

    logger.info(f"Loading scGPT embeddings from: {input_path}")
    adata = sc.read_h5ad(input_path)

    logger.info(f"Original shape: {adata.shape[0]} cells x {adata.shape[1]} features")

    # Check source_dataset column
    if 'source_dataset' in adata.obs.columns:
        logger.info("\nSource dataset distribution (before filtering):")
        source_counts = adata.obs['source_dataset'].value_counts(dropna=False)
        for src, count in source_counts.items():
            logger.info(f"  {src}: {count:,} cells")

        # Filter: keep only cells where source_dataset is NaN (platelet data)
        mask = adata.obs['source_dataset'].isna()
        adata_clean = adata[mask].copy()

        logger.info(f"\nAfter filtering (keeping source_dataset=NaN):")
        logger.info(f"  Kept: {adata_clean.shape[0]:,} cells")
        logger.info(f"  Removed: {adata.shape[0] - adata_clean.shape[0]:,} cells")

    else:
        logger.warning("'source_dataset' column not found. Trying alternative filtering...")

        # Alternative: filter by Disease column (platelet data has Disease values)
        if 'Disease' in adata.obs.columns:
            mask = ~adata.obs['Disease'].isna()
            adata_clean = adata[mask].copy()
            logger.info(f"Filtered by Disease column: {adata_clean.shape[0]:,} cells")
        else:
            logger.error("Cannot determine filtering criteria!")
            return None

    # Verify the cleaned data
    logger.info(f"\nCleaned dataset shape: {adata_clean.shape[0]} cells x {adata_clean.shape[1]} features")

    if 'Disease' in adata_clean.obs.columns:
        logger.info("\nDisease distribution (after filtering):")
        disease_counts = adata_clean.obs['Disease'].value_counts(dropna=False)
        for disease, count in disease_counts.items():
            logger.info(f"  {disease}: {count:,} cells")

    if 'Data_NO' in adata_clean.obs.columns:
        logger.info(f"\nData_NO unique values: {adata_clean.obs['Data_NO'].nunique()}")
        logger.info(f"  {adata_clean.obs['Data_NO'].unique().tolist()}")

    # Save cleaned embeddings
    logger.info(f"\nSaving cleaned embeddings to: {output_path}")
    adata_clean.write_h5ad(output_path)

    logger.info("Done!")
    return adata_clean


def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")

    input_path = base_dir / "02_EMBEDDINGS/scgpt/covid_scgpt.h5ad"
    output_path = base_dir / "02_EMBEDDINGS/scgpt/covid_scgpt_clean.h5ad"

    # Also check alternative path
    if not input_path.exists():
        input_path = base_dir / "02_EMBEDDINGS/scgpt/platelet_scgpt_combined.h5ad"
        output_path = base_dir / "02_EMBEDDINGS/scgpt/platelet_scgpt_clean.h5ad"

    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return 1

    clean_scgpt_embeddings(input_path, output_path)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
