#!/usr/bin/env python3
"""
Quick UCE Files Inspector
Check structure and content of UCE embedding files
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# UCE embeddings path
uce_path = Path("/bigdata/godziklab/shared/Xinru/302006/embeddings/uce")

# Key files to inspect
files_to_check = {
    'HEK293T Epithelial': 'hek293t_epithelial_uce_adata.h5ad',
    'HCT116 Epithelial': 'hct116_epithelial_uce_adata.h5ad', 
    'Mouse RETNLA': 'mouse_epithelial_retnla_ko_wt_uce_adata.h5ad',
    'HEK293T Platelets': 'hek293t_platelets_uce_adata.h5ad',
    'HCT116 Platelets': 'hct116_platelets_uce_adata.h5ad',
    'COVID Platelets': 'human_platelet_covid_severity_uce_adata.h5ad'
}

print("🔍 UCE FILES INSPECTION REPORT")
print("=" * 50)

for name, filename in files_to_check.items():
    filepath = uce_path / filename
    
    if filepath.exists():
        print(f"\n📁 {name}")
        print(f"   File: {filename}")
        
        try:
            adata = sc.read_h5ad(filepath)
            
            print(f"   📊 Shape: {adata.shape} (cells × genes)")
            print(f"   📋 .obs columns: {list(adata.obs.columns)}")
            print(f"   🔢 .obsm keys: {list(adata.obsm.keys())}")
            print(f"   📈 .var columns: {list(adata.var.columns)}")
            
            # Check main data matrix
            if hasattr(adata, 'X') and adata.X is not None:
                print(f"   🧮 .X shape: {adata.X.shape}")
                print(f"   🧮 .X type: {type(adata.X)}")
            
            # Look for embeddings in .obsm
            embedding_keys = []
            for key in adata.obsm.keys():
                if any(term in key.lower() for term in ['uce', 'embed', 'x_']):
                    embedding_keys.append(key)
                    print(f"   ✨ {key}: {adata.obsm[key].shape}")
            
            # Check specific metadata columns
            if 'cell_type' in adata.obs.columns:
                unique_types = adata.obs['cell_type'].unique()
                print(f"   🧬 Cell types ({len(unique_types)}): {list(unique_types)[:5]}...")
            
            # Check for COVID severity if applicable
            if 'covid' in filename:
                severity_cols = [col for col in adata.obs.columns 
                               if any(term in col.lower() for term in ['sever', 'covid', 'stage', 'progress'])]
                if severity_cols:
                    print(f"   🦠 Severity columns: {severity_cols}")
                    for col in severity_cols[:2]:  # Show first 2
                        unique_vals = adata.obs[col].unique()
                        print(f"      {col}: {list(unique_vals)}")
            
            # Check for mouse-specific metadata
            if 'mouse' in filename:
                mouse_cols = [col for col in adata.obs.columns 
                             if any(term in col.lower() for term in ['ko', 'wt', 'genotype', 'condition'])]
                if mouse_cols:
                    print(f"   🐭 Mouse-specific columns: {mouse_cols}")
                    for col in mouse_cols[:2]:
                        unique_vals = adata.obs[col].unique()
                        print(f"      {col}: {list(unique_vals)}")
            
        except Exception as e:
            print(f"   ❌ Error reading file: {e}")
    else:
        print(f"\n❌ {name}: File not found - {filename}")

print(f"\n🎯 SUMMARY:")
print(f"UCE embeddings directory: {uce_path}")
print("Ready to run comprehensive analysis!")