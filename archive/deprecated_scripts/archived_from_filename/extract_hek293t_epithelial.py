#!/usr/bin/env python3
# Extract Epithelial Cells from HEK293T

import scanpy as sc
import numpy as np
from sklearn.preprocessing import StandardScaler

# Load HEK293T processed data
print("Loading HEK293T processed data...")
adata = sc.read_h5ad("/bigdata/godziklab/shared/Xinru/302006/data/processed/HEK293T_processed.h5ad")
print(f"Loaded: {adata.shape}")

# Epithelial cell markers
markers = [
    'EPCAM',    # Epithelial cell adhesion molecule (primary epithelial marker)
    'CDH1',     # E-cadherin (epithelial adhesion)
    'KRT8',     # Keratin 8 (epithelial cytoskeleton)
    'KRT18',    # Keratin 18 (epithelial cytoskeleton)
    'KRT19',    # Keratin 19 (epithelial marker)
    'CLDN1',    # Claudin-1 (tight junction protein)
    'TJP1',     # Tight junction protein 1 (ZO-1)
    'OCLN',     # Occludin (tight junction)
    'DSP',      # Desmoplakin (desmosome)
    'PKP1'      # Plakophilin 1 (desmosome)
]

available_markers = [m for m in markers if m in adata.var_names]
print(f"Available epithelial markers: {available_markers}")

if len(available_markers) >= 3:
    # Calculate epithelial score
    marker_expression = []
    for marker in available_markers:
        if hasattr(adata.X, 'toarray'):
            expr = adata[:, marker].X.toarray().flatten()
        else:
            expr = adata[:, marker].X.flatten()
        marker_expression.append(expr)
    
    marker_matrix = np.column_stack(marker_expression)
    scaler = StandardScaler()
    marker_matrix_norm = scaler.fit_transform(marker_matrix)
    epithelial_score = np.mean(marker_matrix_norm, axis=1)
    
    # Select top 20% as epithelial cells (HEK293T are epithelial-derived)
    threshold = np.percentile(epithelial_score, 80)
    is_epithelial = epithelial_score >= threshold
    
    # Quality filter
    if hasattr(adata.X, 'toarray'):
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    else:
        total_counts = adata.X.sum(axis=1)
    
    min_counts = np.percentile(total_counts, 10)
    high_quality = total_counts >= min_counts
    
    final_mask = is_epithelial & high_quality
    
    # Extract epithelial cells
    epithelial = adata[final_mask].copy()
    epithelial.obs['cell_type'] = 'epithelial'
    epithelial.obs['epithelial_score'] = epithelial_score[final_mask]
    epithelial.obs['source_dataset'] = 'HEK293T'
    
    # Save to raw directory
    output_path = "/bigdata/godziklab/shared/Xinru/302006/data/raw/hek293t_epithelial.h5ad"
    epithelial.write_h5ad(output_path)
    print(f"✅ Extracted {epithelial.n_obs:,} epithelial cells from HEK293T")
    print(f"Saved to: {output_path}")
    
    # Show score distribution
    print(f"Epithelial score range: {epithelial_score.min():.3f} to {epithelial_score.max():.3f}")
    print(f"Selected threshold: {threshold:.3f}")
    
else:
    print("❌ Not enough epithelial markers found")
    print(f"Found only: {available_markers}")