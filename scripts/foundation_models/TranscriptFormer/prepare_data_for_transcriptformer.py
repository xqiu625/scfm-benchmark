#!/usr/bin/env python3
"""
Prepare Data for TranscriptFormer

Purpose: Convert raw count data to TranscriptFormer-compatible format
Requirements:
  - Raw counts (no normalization)
  - Ensembl gene IDs in adata.var
  - H5AD format

Author: Project 302006
Date: November 2025
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
        logging.FileHandler(f'prepare_transcriptformer_data_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class TranscriptFormerDataPrep:
    """Prepare data for TranscriptFormer processing"""

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_dir = Path(base_dir)
        self.raw_dir = self.base_dir / "01_DATA/raw"
        self.output_dir = self.base_dir / "01_DATA/transcriptformer_ready"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_gene_mapping(self):
        """Load or create gene symbol to Ensembl ID mapping"""
        logger.info("Loading gene symbol to Ensembl ID mapping...")

        # Try to load existing mapping
        mapping_file = self.base_dir / "reference/gene_mapping_ensembl.csv"

        if mapping_file.exists():
            logger.info(f"  Loading existing mapping from {mapping_file}")
            mapping = pd.read_csv(mapping_file, index_col=0)
            return mapping

        # If no mapping exists, we'll need to create one using biomaRt or similar
        logger.warning("  No gene mapping file found!")
        logger.info("  TranscriptFormer requires Ensembl IDs")
        logger.info("  Options:")
        logger.info("    1. Use mygene.info API to map gene symbols")
        logger.info("    2. Use existing Ensembl annotations")
        logger.info("    3. Skip genes without Ensembl IDs")

        return None

    def add_ensembl_ids(self, adata):
        """Add Ensembl IDs to gene metadata"""
        logger.info("  Adding Ensembl IDs...")

        # Check if already has ensembl_id column
        if 'ensembl_id' in adata.var.columns:
            logger.info("    adata.var already has 'ensembl_id' column")
            n_missing = adata.var['ensembl_id'].isna().sum()
            logger.info(f"    Genes with Ensembl IDs: {len(adata.var) - n_missing}/{len(adata.var)}")
            return adata

        # Try to map using mygene if possible
        try:
            import mygene
            mg = mygene.MyGeneInfo()

            logger.info("    Querying mygene.info for Ensembl IDs...")
            gene_symbols = adata.var_names.tolist()

            # Query in batches to avoid timeout
            batch_size = 1000
            all_results = []

            for i in range(0, len(gene_symbols), batch_size):
                batch = gene_symbols[i:i+batch_size]
                logger.info(f"    Querying batch {i//batch_size + 1}/{(len(gene_symbols)-1)//batch_size + 1}")
                results = mg.querymany(batch, scopes='symbol', fields='ensembl.gene', species='human')
                all_results.extend(results)

            # Extract Ensembl IDs (take first match for duplicates)
            ensembl_ids = []
            for i, result in enumerate(all_results):
                if 'ensembl' in result:
                    if isinstance(result['ensembl'], list):
                        # Multiple matches - take first one
                        ensembl_ids.append(result['ensembl'][0]['gene'])
                    else:
                        ensembl_ids.append(result['ensembl']['gene'])
                else:
                    ensembl_ids.append(None)

            # Ensure length matches
            if len(ensembl_ids) != len(gene_symbols):
                logger.warning(f"    Ensembl ID count mismatch: {len(ensembl_ids)} vs {len(gene_symbols)} genes")
                logger.warning(f"    Truncating to match gene count")
                ensembl_ids = ensembl_ids[:len(gene_symbols)]

            adata.var['ensembl_id'] = ensembl_ids

            n_mapped = sum(1 for eid in ensembl_ids if eid is not None)
            logger.info(f"    Mapped {n_mapped}/{len(gene_symbols)} genes to Ensembl IDs")

        except ImportError:
            logger.warning("    mygene package not installed")
            logger.warning("    Install with: pip install mygene")
            logger.warning("    Proceeding without Ensembl ID mapping")
            adata.var['ensembl_id'] = None

        return adata

    def prepare_dataset(self, dataset_name, input_file):
        """Prepare a single dataset for TranscriptFormer"""
        logger.info(f"\n{'='*80}")
        logger.info(f"Preparing: {dataset_name}")
        logger.info(f"{'='*80}")

        # Load data
        logger.info(f"Loading {input_file.name}...")
        adata = sc.read_h5ad(input_file)
        logger.info(f"  Shape: {adata.shape}")
        logger.info(f"  Current .var columns: {list(adata.var.columns)}")

        # Basic QC (TranscriptFormer expects filtered data)
        logger.info("\nBasic quality control...")
        logger.info(f"  Cells before QC: {adata.shape[0]}")

        # Calculate QC metrics if not present
        if 'n_genes' not in adata.obs.columns:
            sc.pp.calculate_qc_metrics(adata, inplace=True)

        # Filter cells (standard thresholds from guide)
        # - Minimum 200 genes per cell
        # - Maximum 5000 genes per cell
        # - Maximum 20% mitochondrial reads
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_cells(adata, max_genes=5000)

        if 'pct_counts_mt' in adata.obs.columns:
            adata = adata[adata.obs['pct_counts_mt'] < 20].copy()

        logger.info(f"  Cells after QC: {adata.shape[0]}")

        # Ensure we have raw counts
        logger.info("\nVerifying raw counts...")
        if hasattr(adata, 'raw') and adata.raw is not None:
            logger.info("  Using adata.raw for counts")
            # Create new adata with raw counts
            adata = sc.AnnData(
                X=adata.raw.X,
                obs=adata.obs,
                var=adata.raw.var
            )

        # Check if counts are integers (raw counts)
        if not np.allclose(adata.X.data if hasattr(adata.X, 'data') else adata.X,
                          np.round(adata.X.data if hasattr(adata.X, 'data') else adata.X)):
            logger.warning("  Warning: Data may be normalized (not integers)")
            logger.warning("  TranscriptFormer expects raw counts!")
        else:
            logger.info("  ✓ Data appears to be raw counts")

        # Add Ensembl IDs
        adata = self.add_ensembl_ids(adata)

        # Filter out genes without VALID Ensembl IDs (TranscriptFormer requirement)
        # IMPORTANT: Only keep IDs that start with "ENSG" (proper Ensembl gene IDs)
        # Reject Entrez IDs (plain numbers) and other non-Ensembl formats
        logger.info("\nFiltering genes...")
        n_before = adata.shape[1]

        # Check which genes have valid Ensembl IDs (start with "ENSG")
        def is_valid_ensembl(eid):
            if pd.isna(eid):
                return False
            eid_str = str(eid)
            return eid_str.startswith('ENSG')

        genes_with_valid_ensembl = adata.var['ensembl_id'].apply(is_valid_ensembl)
        n_valid = genes_with_valid_ensembl.sum()
        n_invalid = (~genes_with_valid_ensembl).sum()

        logger.info(f"  Genes before filtering: {n_before}")
        logger.info(f"  Genes with valid Ensembl IDs (ENSG*): {n_valid}")
        logger.info(f"  Genes with invalid/missing IDs (will be removed): {n_invalid}")

        # Keep only genes with valid Ensembl IDs
        adata = adata[:, genes_with_valid_ensembl].copy()
        logger.info(f"  Genes after filtering: {adata.shape[1]}")

        # Check for and remove duplicate Ensembl IDs
        # IMPORTANT: Strip version numbers first (e.g., ENSG00000123456.1 → ENSG00000123456)
        # This matches what TranscriptFormer does internally
        logger.info("\nChecking for duplicate Ensembl IDs...")
        ensembl_ids = adata.var['ensembl_id'].values

        # Strip version numbers from Ensembl IDs
        ensembl_ids_no_version = np.array([str(eid).split('.')[0] if pd.notna(eid) else eid for eid in ensembl_ids])

        # Check for duplicates after version removal
        unique_ids, unique_indices = np.unique(ensembl_ids_no_version, return_index=True)

        n_duplicates = len(ensembl_ids_no_version) - len(unique_ids)
        if n_duplicates > 0:
            logger.warning(f"  Found {n_duplicates} duplicate Ensembl IDs (after removing version numbers)")
            logger.warning(f"  Keeping first occurrence of each duplicate")

            # Keep only first occurrence of each Ensembl ID
            adata = adata[:, unique_indices].copy()

            # Update ensembl_id column with version-stripped IDs
            adata.var['ensembl_id'] = ensembl_ids_no_version[unique_indices]
            logger.info(f"  Genes after removing duplicates: {adata.shape[1]}")
        else:
            logger.info(f"  No duplicate Ensembl IDs found")
            # Still strip version numbers even if no duplicates
            adata.var['ensembl_id'] = ensembl_ids_no_version

        # CRITICAL: Set var_names (index) to ensembl_id for TranscriptFormer
        # TranscriptFormer uses var_names as gene identifiers when checking for duplicates
        # Keeping original gene symbols as a separate column
        logger.info("\nSetting var index to Ensembl IDs...")
        adata.var['gene_symbol'] = adata.var_names.values  # Preserve original gene symbols
        adata.var_names = adata.var['ensembl_id'].values
        adata.var_names_make_unique()  # Ensure uniqueness (adds suffix if needed)
        logger.info(f"  var_names now set to Ensembl IDs")

        # Store raw data
        adata.raw = adata

        # Add required TranscriptFormer metadata columns
        logger.info("\nAdding required TranscriptFormer metadata...")

        # Add 'assay' column (required by TranscriptFormer tokenizer)
        if 'assay' not in adata.obs.columns:
            # Use '10x 3' v3' as default (most common scRNA-seq platform)
            # Or 'unknown' if assay type is truly unknown
            adata.obs['assay'] = '10x 3\' v3'
            logger.info("  Added 'assay' column: 10x 3' v3")

        # Add 'organism' column if not present (TranscriptFormer may need this)
        if 'organism' not in adata.obs.columns:
            adata.obs['organism'] = 'Homo sapiens'
            logger.info("  Added 'organism' column: Homo sapiens")

        # Add 'tissue' column if not present
        if 'tissue' not in adata.obs.columns:
            # For platelets data, use appropriate tissue type
            if dataset_name == 'covid':
                adata.obs['tissue'] = 'blood'
            else:
                adata.obs['tissue'] = 'cell line'
            logger.info(f"  Added 'tissue' column: {adata.obs['tissue'].iloc[0]}")

        # Convert all .obs columns to strings to avoid H5AD serialization errors
        logger.info("\nConverting metadata to H5AD-compatible format...")
        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                adata.obs[col] = adata.obs[col].astype(str)

        # Convert all .var columns to strings
        for col in adata.var.columns:
            if adata.var[col].dtype == 'object':
                adata.var[col] = adata.var[col].astype(str)

        # Save prepared data
        output_file = self.output_dir / f"{dataset_name}_transcriptformer_ready.h5ad"
        logger.info(f"\nSaving to: {output_file}")
        adata.write_h5ad(output_file)

        logger.info(f"✓ {dataset_name} prepared successfully")
        logger.info(f"  Final shape: {adata.shape}")
        if 'ensembl_id' in adata.var.columns:
            n_with_ensembl = adata.var['ensembl_id'].notna().sum()
            logger.info(f"  Genes with Ensembl IDs: {n_with_ensembl}/{adata.shape[1]}")

        return adata

    def prepare_all_datasets(self):
        """Prepare all three datasets"""
        logger.info("="*80)
        logger.info("TRANSCRIPTFORMER DATA PREPARATION")
        logger.info("="*80)

        datasets = {
            'covid': 'human_platelet_covid_severity_raw_counts.h5ad',
            'hek293t': 'HEK293T_platelets_raw_counts.h5ad',
            'hct116': 'HCT116_platelets_raw_counts.h5ad'
        }

        results = {}
        for dataset_name, filename in datasets.items():
            input_file = self.raw_dir / filename

            if not input_file.exists():
                logger.error(f"File not found: {input_file}")
                continue

            try:
                adata = self.prepare_dataset(dataset_name, input_file)
                results[dataset_name] = adata
            except Exception as e:
                logger.error(f"Failed to prepare {dataset_name}: {str(e)}")
                continue

        # Summary
        logger.info("\n" + "="*80)
        logger.info("PREPARATION COMPLETE")
        logger.info("="*80)
        logger.info(f"\nOutput directory: {self.output_dir}")
        logger.info("\nPrepared datasets:")
        for dataset_name, adata in results.items():
            logger.info(f"  {dataset_name}: {adata.shape[0]} cells x {adata.shape[1]} genes")

        logger.info("\nNext steps:")
        logger.info("  1. Verify Ensembl ID coverage")
        logger.info("  2. Run TranscriptFormer embedding generation")
        logger.info("  3. Compare with UCE and scGPT embeddings")

        return results


def main():
    """Main execution"""
    try:
        prep = TranscriptFormerDataPrep()
        results = prep.prepare_all_datasets()
        return 0

    except Exception as e:
        logger.error(f"\nERROR: {str(e)}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
