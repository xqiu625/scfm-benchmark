#!/usr/bin/env python3
"""
Generate TranscriptFormer Embeddings using CLI

Purpose: Generate embeddings using TF-Sapiens model via TranscriptFormer CLI
Key features:
  - Uses transcriptformer CLI for inference
  - Automatic preprocessing (no manual normalization needed)
  - Mixed precision (fp16) for memory efficiency
  - Batch processing for 47K-89K cells

Author: Project 302006
Date: November 2025
"""

import scanpy as sc
import subprocess
from pathlib import Path
import logging
from datetime import datetime
import sys
import shutil

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'transcriptformer_embeddings_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class TranscriptFormerEmbedding:
    """Generate TranscriptFormer embeddings using CLI"""

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_dir = Path(base_dir)
        self.input_dir = self.base_dir / "01_DATA/transcriptformer_ready"
        self.output_dir = self.base_dir / "02_EMBEDDINGS/transcriptformer"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Model checkpoint path
        self.checkpoint_path = Path.cwd() / "checkpoints/tf_sapiens"

        # Check if model exists
        if not self.checkpoint_path.exists():
            logger.error(f"Model not found at {self.checkpoint_path}")
            logger.error("Run: transcriptformer download tf-sapiens")
            raise FileNotFoundError(f"Model checkpoint not found: {self.checkpoint_path}")

        logger.info(f"Model checkpoint: {self.checkpoint_path}")
        logger.info(f"Input directory: {self.input_dir}")
        logger.info(f"Output directory: {self.output_dir}")

    def run_inference(self, input_file, output_path, batch_size=16, emb_type="cell"):
        """Run TranscriptFormer inference using CLI"""

        cmd = [
            "transcriptformer", "inference",
            "--checkpoint-path", str(self.checkpoint_path),
            "--data-file", str(input_file),
            "--output-path", str(output_path),
            "--batch-size", str(batch_size),
            "--emb-type", emb_type,
            "--gene-col-name", "ensembl_id",
            "--precision", "16-mixed"
            # NOTE: --remove-duplicate-genes has a bug in TranscriptFormer
            # Duplicates are handled in data preparation instead
        ]

        logger.info(f"Running command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            logger.info("Inference completed successfully")
            if result.stdout:
                logger.info(f"Output:\n{result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            logger.error(f"Inference failed with exit code {e.returncode}")
            if e.stdout:
                logger.error(f"STDOUT:\n{e.stdout}")
            if e.stderr:
                logger.error(f"STDERR:\n{e.stderr}")
            raise

    def process_dataset(self, dataset_name, batch_size=16):
        """Process a single dataset"""
        logger.info(f"\n{'='*80}")
        logger.info(f"Processing: {dataset_name}")
        logger.info(f"{'='*80}")

        # Input file
        input_file = self.input_dir / f"{dataset_name}_transcriptformer_ready.h5ad"
        if not input_file.exists():
            logger.error(f"Input file not found: {input_file}")
            return None

        logger.info(f"Loading {input_file.name}...")
        adata = sc.read_h5ad(input_file)
        logger.info(f"  Shape: {adata.shape}")

        # Check for Ensembl IDs
        if 'ensembl_id' not in adata.var.columns:
            logger.error("  Missing 'ensembl_id' column in adata.var!")
            logger.error("  Run prepare_data_for_transcriptformer.py first")
            return None

        n_with_ensembl = adata.var['ensembl_id'].notna().sum()
        logger.info(f"  Genes with Ensembl IDs: {n_with_ensembl}/{adata.shape[1]}")

        # Create temporary output directory for this dataset
        temp_output_dir = self.output_dir / f"temp_{dataset_name}"
        temp_output_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Run inference
            logger.info(f"\nGenerating embeddings (batch_size={batch_size})...")
            self.run_inference(
                input_file=input_file,
                output_path=temp_output_dir,
                batch_size=batch_size,
                emb_type="cell"
            )

            # Load the generated embeddings
            embeddings_file = temp_output_dir / "embeddings.h5ad"
            if not embeddings_file.exists():
                logger.error(f"Embeddings file not found: {embeddings_file}")
                return None

            logger.info(f"Loading embeddings from {embeddings_file}")
            adata_emb = sc.read_h5ad(embeddings_file)

            # Extract embeddings and add to original adata
            if 'embeddings' in adata_emb.obsm:
                adata.obsm['X_transcriptformer'] = adata_emb.obsm['embeddings']
            elif 'X_transcriptformer' in adata_emb.obsm:
                adata.obsm['X_transcriptformer'] = adata_emb.obsm['X_transcriptformer']
            else:
                logger.error(f"No embeddings found in {embeddings_file}")
                logger.error(f"Available keys: {list(adata_emb.obsm.keys())}")
                return None

            # Save with embeddings
            output_file = self.output_dir / f"{dataset_name}_transcriptformer.h5ad"
            logger.info(f"\nSaving to: {output_file}")
            adata.write_h5ad(output_file)

            emb_shape = adata.obsm['X_transcriptformer'].shape
            logger.info(f"✓ {dataset_name} processed successfully")
            logger.info(f"  Cells: {emb_shape[0]:,}")
            logger.info(f"  Embedding dimensions: {emb_shape[1]}")

            # Clean up temporary directory
            shutil.rmtree(temp_output_dir)

            return adata

        except Exception as e:
            logger.error(f"Failed to process {dataset_name}: {str(e)}")
            # Clean up temporary directory on error
            if temp_output_dir.exists():
                shutil.rmtree(temp_output_dir)
            raise

    def process_all_datasets(self, batch_size=16):
        """Process all datasets"""
        logger.info("="*80)
        logger.info("TRANSCRIPTFORMER EMBEDDING GENERATION")
        logger.info("="*80)

        datasets = ['covid', 'hek293t', 'hct116']
        results = {}

        for dataset_name in datasets:
            try:
                adata = self.process_dataset(dataset_name, batch_size)
                if adata is not None:
                    results[dataset_name] = adata
            except Exception as e:
                logger.error(f"Failed to process {dataset_name}: {str(e)}")
                continue

        # Summary
        logger.info("\n" + "="*80)
        logger.info("EMBEDDING GENERATION COMPLETE")
        logger.info("="*80)
        logger.info(f"\nOutput directory: {self.output_dir}")
        logger.info("\nGenerated embeddings:")
        for dataset_name, adata in results.items():
            emb_shape = adata.obsm['X_transcriptformer'].shape
            logger.info(f"  {dataset_name}: {emb_shape[0]} cells x {emb_shape[1]} dimensions")

        logger.info("\nNext steps:")
        logger.info("  1. Merge with UCE and scGPT embeddings")
        logger.info("  2. Run three-model comparison analysis")
        logger.info("  3. Evaluate TranscriptFormer performance")

        return results


def main():
    """Main execution"""
    try:
        # Parse command line arguments
        import argparse
        parser = argparse.ArgumentParser(
            description='Generate TranscriptFormer embeddings using CLI'
        )
        parser.add_argument('--batch-size', type=int, default=16,
                          help='Batch size for inference (default: 16)')
        parser.add_argument('--dataset', type=str, default=None,
                          help='Process specific dataset only (covid/hek293t/hct116)')
        args = parser.parse_args()

        embedder = TranscriptFormerEmbedding()

        if args.dataset:
            # Process single dataset
            embedder.process_dataset(args.dataset, args.batch_size)
        else:
            # Process all datasets
            embedder.process_all_datasets(args.batch_size)

        return 0

    except Exception as e:
        logger.error(f"\nERROR: {str(e)}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
