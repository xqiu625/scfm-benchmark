#!/usr/bin/env python3
"""
UCE (Universal Cell Embeddings) Generation Script

Generates UCE embeddings for single-cell RNA-seq datasets using the 33-layer
transformer model with ESM2 protein embeddings for gene tokenization.

Usage:
    python generate_uce_embeddings.py --base-dir /path/to/project --dataset covid
    python generate_uce_embeddings.py --base-dir /path/to/project --dataset all
    python generate_uce_embeddings.py --input data.h5ad --output embeddings/

Requirements:
    - UCE repository cloned with eval_single_anndata.py
    - GPU with sufficient memory (recommended: 24GB+)
    - Conda environment with UCE dependencies

Author: Project 302006
Date: January 2026
"""

import argparse
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Default dataset configurations
DATASET_CONFIGS = {
    'covid': {
        'input': 'data/raw/human_platelet_covid_severity.h5ad',
        'output_dir': 'embeddings/uce_embeddings/covid',
        'species': 'human',
        'description': 'COVID-19 platelet dataset'
    },
    'hek293t': {
        'input': 'data/raw/hek293t_perturb_seq.h5ad',
        'output_dir': 'embeddings/uce_embeddings/hek293t',
        'species': 'human',
        'description': 'HEK293T Perturb-seq dataset'
    },
    'hct116': {
        'input': 'data/raw/hct116_perturb_seq.h5ad',
        'output_dir': 'embeddings/uce_embeddings/hct116',
        'species': 'human',
        'description': 'HCT116 Perturb-seq dataset'
    }
}


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate UCE embeddings for single-cell data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a specific dataset
  python generate_uce_embeddings.py --base-dir /path/to/project --dataset covid

  # Process all datasets
  python generate_uce_embeddings.py --base-dir /path/to/project --dataset all

  # Process a custom file
  python generate_uce_embeddings.py --input data.h5ad --output embeddings/ --species human
        """
    )

    # Mode 1: Use predefined dataset configurations
    parser.add_argument('--base-dir', '-b', type=str,
                        help='Base directory containing project structure')
    parser.add_argument('--dataset', '-d', type=str,
                        choices=['covid', 'hek293t', 'hct116', 'all'],
                        help='Dataset to process (or "all" for all datasets)')

    # Mode 2: Direct file specification
    parser.add_argument('--input', '-i', type=str,
                        help='Input h5ad file path (alternative to --dataset)')
    parser.add_argument('--output', '-o', type=str,
                        help='Output directory for embeddings')

    # Common parameters
    parser.add_argument('--species', '-s', type=str, default='human',
                        help='Species for gene mapping (default: human)')
    parser.add_argument('--batch-size', type=int, default=200,
                        help='Batch size for inference (default: 200)')
    parser.add_argument('--uce-dir', type=str,
                        help='Path to UCE repository (default: {base_dir}/scfm/UCE)')
    parser.add_argument('--conda-env', type=str,
                        help='Conda environment path or name')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print commands without executing')

    return parser.parse_args()


def validate_paths(input_path: Path, output_dir: Path, uce_dir: Path) -> bool:
    """Validate that required paths exist."""
    errors = []

    if not input_path.exists():
        errors.append(f"Input file not found: {input_path}")

    if not uce_dir.exists():
        errors.append(f"UCE directory not found: {uce_dir}")

    eval_script = uce_dir / "eval_single_anndata.py"
    if not eval_script.exists():
        errors.append(f"UCE eval script not found: {eval_script}")

    if errors:
        for error in errors:
            logger.error(error)
        return False

    return True


def run_uce_embedding(
    input_path: Path,
    output_dir: Path,
    uce_dir: Path,
    species: str = 'human',
    batch_size: int = 200,
    conda_env: str = None,
    dry_run: bool = False
) -> bool:
    """
    Run UCE embedding generation for a single dataset.

    Args:
        input_path: Path to input h5ad file
        output_dir: Directory for output embeddings
        uce_dir: Path to UCE repository
        species: Species name for gene mapping
        batch_size: Batch size for inference
        conda_env: Conda environment path or name
        dry_run: If True, print command without executing

    Returns:
        bool: True if successful, False otherwise
    """
    eval_script = uce_dir / "eval_single_anndata.py"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if output already exists
    expected_output = output_dir / f"{input_path.stem}_uce.h5ad"
    if expected_output.exists():
        logger.info(f"Output already exists: {expected_output}")
        return True

    # Build command
    cmd = [
        "python", str(eval_script),
        "--adata_path", str(input_path),
        "--dir", str(output_dir),
        "--species", species,
        "--batch_size", str(batch_size)
    ]

    logger.info("=" * 60)
    logger.info("UCE EMBEDDING GENERATION")
    logger.info("=" * 60)
    logger.info(f"Input:      {input_path}")
    logger.info(f"Output:     {output_dir}")
    logger.info(f"Species:    {species}")
    logger.info(f"Batch size: {batch_size}")
    logger.info(f"Command:    {' '.join(cmd)}")
    logger.info("=" * 60)

    if dry_run:
        logger.info("DRY RUN - Command not executed")
        return True

    try:
        # Change to UCE directory for relative imports
        start_time = datetime.now()
        logger.info(f"Starting at {start_time}")

        result = subprocess.run(
            cmd,
            cwd=str(uce_dir),
            check=True,
            capture_output=False,  # Let output stream to console
            text=True
        )

        end_time = datetime.now()
        duration = end_time - start_time
        logger.info(f"Completed in {duration}")

        # Verify output was created
        if expected_output.exists():
            logger.info(f"Output verified: {expected_output}")

            # Load and report basic stats
            try:
                import scanpy as sc
                adata = sc.read_h5ad(expected_output)
                logger.info(f"Embedding shape: {adata.shape}")
                if 'X_uce' in adata.obsm:
                    emb = adata.obsm['X_uce']
                    logger.info(f"Embedding dim: {emb.shape[1]}")
                    logger.info(f"Embedding range: [{np.min(emb):.4f}, {np.max(emb):.4f}]")
            except ImportError:
                pass

            return True
        else:
            logger.warning(f"Expected output not found: {expected_output}")
            # Check for alternative output patterns
            outputs = list(output_dir.glob("*.h5ad"))
            if outputs:
                logger.info(f"Found outputs: {[o.name for o in outputs]}")
                return True
            return False

    except subprocess.CalledProcessError as e:
        logger.error(f"UCE embedding generation failed with return code {e.returncode}")
        if e.stderr:
            logger.error(f"Error output: {e.stderr}")
        return False

    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return False


def main():
    args = parse_args()

    # Determine mode of operation
    if args.input:
        # Mode 2: Direct file specification
        if not args.output:
            logger.error("--output is required when using --input")
            return 1

        input_path = Path(args.input)
        output_dir = Path(args.output)

        # Determine UCE directory
        if args.uce_dir:
            uce_dir = Path(args.uce_dir)
        elif args.base_dir:
            uce_dir = Path(args.base_dir) / "scfm" / "UCE"
        else:
            logger.error("--uce-dir or --base-dir is required")
            return 1

        if not validate_paths(input_path, output_dir, uce_dir):
            return 1

        success = run_uce_embedding(
            input_path=input_path,
            output_dir=output_dir,
            uce_dir=uce_dir,
            species=args.species,
            batch_size=args.batch_size,
            conda_env=args.conda_env,
            dry_run=args.dry_run
        )
        return 0 if success else 1

    elif args.base_dir and args.dataset:
        # Mode 1: Use predefined dataset configurations
        base_dir = Path(args.base_dir)

        if not base_dir.exists():
            logger.error(f"Base directory not found: {base_dir}")
            return 1

        # Determine UCE directory
        uce_dir = Path(args.uce_dir) if args.uce_dir else base_dir / "scfm" / "UCE"

        # Determine which datasets to process
        if args.dataset == 'all':
            datasets = list(DATASET_CONFIGS.keys())
        else:
            datasets = [args.dataset]

        logger.info(f"Processing datasets: {datasets}")

        results = {}
        for dataset_name in datasets:
            config = DATASET_CONFIGS[dataset_name]
            logger.info(f"\n{'='*60}")
            logger.info(f"Processing: {config['description']}")
            logger.info(f"{'='*60}")

            input_path = base_dir / config['input']
            output_dir = base_dir / config['output_dir']
            species = config.get('species', args.species)

            if not validate_paths(input_path, output_dir, uce_dir):
                logger.error(f"Skipping {dataset_name} due to validation errors")
                results[dataset_name] = False
                continue

            success = run_uce_embedding(
                input_path=input_path,
                output_dir=output_dir,
                uce_dir=uce_dir,
                species=species,
                batch_size=args.batch_size,
                conda_env=args.conda_env,
                dry_run=args.dry_run
            )
            results[dataset_name] = success

        # Summary
        logger.info("\n" + "=" * 60)
        logger.info("SUMMARY")
        logger.info("=" * 60)
        for dataset_name, success in results.items():
            status = "SUCCESS" if success else "FAILED"
            logger.info(f"  {dataset_name}: {status}")

        failed = [k for k, v in results.items() if not v]
        if failed:
            logger.error(f"Failed datasets: {failed}")
            return 1
        return 0

    else:
        logger.error("Either (--base-dir and --dataset) or (--input and --output) are required")
        logger.error("Run with --help for usage information")
        return 1


if __name__ == "__main__":
    sys.exit(main())
