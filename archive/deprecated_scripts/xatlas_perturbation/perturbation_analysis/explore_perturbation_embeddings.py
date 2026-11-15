#!/usr/bin/env python3
"""
Explore available perturbation embeddings on the cluster.
Run this first to verify what data is available for perturbation analysis.

Usage:
    python scripts/analysis/explore_perturbation_embeddings.py

Author: Project 302006
Date: December 2025
"""

import os
from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def explore_embeddings(base_dir="/bigdata/godziklab/shared/Xinru/302006"):
    """Explore available embedding files"""
    embedding_dir = Path(base_dir) / "02_EMBEDDINGS"

    print("="*70)
    print("EXPLORING PERTURBATION EMBEDDINGS")
    print("="*70)
    print(f"\nBase directory: {base_dir}")
    print(f"Embedding directory: {embedding_dir}")
    print(f"Directory exists: {embedding_dir.exists()}")

    if not embedding_dir.exists():
        print("\nERROR: Embedding directory not found!")
        return

    # List all subdirectories
    print("\n" + "-"*50)
    print("AVAILABLE EMBEDDING DIRECTORIES:")
    print("-"*50)

    for subdir in sorted(embedding_dir.iterdir()):
        if subdir.is_dir():
            files = list(subdir.glob("*.h5ad"))
            print(f"\n📁 {subdir.name}/")
            if files:
                for f in sorted(files):
                    size_mb = f.stat().st_size / (1024*1024)
                    print(f"   - {f.name} ({size_mb:.1f} MB)")
            else:
                print("   (no .h5ad files)")

    # Check UCE files specifically (primary for perturbation analysis)
    print("\n" + "="*70)
    print("UCE EMBEDDING DETAILS (Primary for Perturbation Analysis)")
    print("="*70)

    uce_paths = {
        'UCE 33-layer COVID': embedding_dir / "uce_33layer/covid_uce_33layer.h5ad",
        'UCE 33-layer HEK293T': embedding_dir / "uce_33layer/hek293t_uce_33layer.h5ad",
        'UCE 33-layer HCT116': embedding_dir / "uce_33layer/hct116_uce_33layer.h5ad",
        'UCE 4-layer COVID': embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad",
        'UCE 4-layer HEK293T': embedding_dir / "uce/human_platelet_hek293t_uce_adata.h5ad",
        'UCE 4-layer HCT116': embedding_dir / "uce/human_platelet_hct116_uce_adata.h5ad",
    }

    for name, path in uce_paths.items():
        print(f"\n{name}:")
        print(f"  Path: {path}")

        if path.exists():
            print(f"  Status: ✅ EXISTS")
            try:
                adata = sc.read_h5ad(path)
                print(f"  Cells: {adata.n_obs}")
                print(f"  Genes: {adata.n_vars}")
                print(f"  Obs columns: {list(adata.obs.columns)[:10]}...")
                print(f"  Obsm keys: {list(adata.obsm.keys())}")

                # Check for severity or perturbation labels
                for col in ['covid_severity', 'severity', 'Category', 'guide_target', 'gene_target', 'perturbation']:
                    if col in adata.obs.columns:
                        unique = adata.obs[col].unique()
                        print(f"  {col}: {len(unique)} unique values")
                        if len(unique) <= 10:
                            print(f"    Values: {list(unique)[:10]}")

                del adata  # Free memory

            except Exception as e:
                print(f"  Error reading: {e}")
        else:
            print(f"  Status: ❌ NOT FOUND")

    # Check scGPT combined file
    print("\n" + "="*70)
    print("scGPT COMBINED FILE")
    print("="*70)

    scgpt_file = embedding_dir / "scgpt/platelet_scgpt_combined.h5ad"
    if scgpt_file.exists():
        print(f"  Path: {scgpt_file}")
        print(f"  Status: ✅ EXISTS")
        try:
            adata = sc.read_h5ad(scgpt_file)
            print(f"  Total cells: {adata.n_obs}")
            print(f"  Obsm keys: {list(adata.obsm.keys())}")

            if 'dataset' in adata.obs.columns:
                print(f"  Datasets: {adata.obs['dataset'].value_counts().to_dict()}")

            del adata
        except Exception as e:
            print(f"  Error: {e}")
    else:
        print(f"  Status: ❌ NOT FOUND")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY - PERTURBATION ANALYSIS READINESS")
    print("="*70)

    uce_33_hek = embedding_dir / "uce_33layer/hek293t_uce_33layer.h5ad"
    uce_33_hct = embedding_dir / "uce_33layer/hct116_uce_33layer.h5ad"
    uce_4_hek = embedding_dir / "uce/human_platelet_hek293t_uce_adata.h5ad"
    uce_4_hct = embedding_dir / "uce/human_platelet_hct116_uce_adata.h5ad"

    if uce_33_hek.exists() or uce_4_hek.exists():
        print("✅ UCE HEK293T embeddings: AVAILABLE")
        print("   - Ready for perturbation analysis")
    else:
        print("❌ UCE HEK293T embeddings: NOT FOUND")

    if uce_33_hct.exists() or uce_4_hct.exists():
        print("✅ UCE HCT116 embeddings: AVAILABLE")
        print("   - Ready for perturbation analysis")
    else:
        print("❌ UCE HCT116 embeddings: NOT FOUND")

    print("\n" + "-"*50)
    print("RECOMMENDED NEXT STEPS:")
    print("-"*50)

    if uce_33_hek.exists() or uce_4_hek.exists():
        print("1. Run perturbation analysis with UCE:")
        print("   python scripts/analysis/perturbation_severity_analysis.py --model UCE")
    else:
        print("1. UCE perturbation embeddings not found. Check:")
        print("   - If embeddings were generated")
        print("   - If paths are correct")

    print("\n2. If using scGPT combined file, ensure proper dataset filtering")
    print("3. TranscriptFormer only supports COVID (no perturbation data)")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Explore perturbation embeddings')
    parser.add_argument('--base-dir', type=str,
                        default="/bigdata/godziklab/shared/Xinru/302006",
                        help='Base directory')
    args = parser.parse_args()

    explore_embeddings(base_dir=args.base_dir)
