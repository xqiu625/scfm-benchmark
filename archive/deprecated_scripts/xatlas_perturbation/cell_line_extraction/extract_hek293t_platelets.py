#!/usr/bin/env python3
# Extract Platelets from HEK293T - Memory Efficient Version

import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import gc

def extract_platelets_chunked():
    """Extract platelets with memory-efficient processing"""
    
    print("Loading HEK293T processed data...")
    
    # Load the processed data
    adata = sc.read_h5ad("/bigdata/godziklab/shared/Xinru/302006/data/processed/HEK293T_processed.h5ad")
    print(f"Loaded: {adata.shape}")
    
    # Platelet markers
    markers = ['ITGA2B', 'ITGB3', 'GP1BA', 'GP1BB', 'THBS1', 'PDGFA', 'PF4', 'PPBP', 'GP9', 'TUBB1']
    available_markers = [m for m in markers if m in adata.var_names]
    print(f"Available platelet markers: {available_markers}")
    
    if len(available_markers) < 3:
        print(f"❌ Only {len(available_markers)} markers found. Need at least 3.")
        return
    
    print("Calculating platelet scores...")
    
    # Process in chunks to save memory
    chunk_size = 10000  # Process 10k cells at a time
    n_cells = adata.n_obs
    n_chunks = (n_cells + chunk_size - 1) // chunk_size
    
    platelet_scores = np.zeros(n_cells)
    
    # Calculate scores chunk by chunk
    for i in range(n_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, n_cells)
        
        print(f"Processing chunk {i+1}/{n_chunks} (cells {start_idx}-{end_idx})")
        
        # Get chunk data
        chunk_data = adata[start_idx:end_idx, available_markers]
        
        # Extract expression
        if hasattr(chunk_data.X, 'toarray'):
            chunk_expr = chunk_data.X.toarray()
        else:
            chunk_expr = chunk_data.X
        
        # Calculate marker scores for this chunk
        if chunk_expr.shape[0] > 1:  # Need multiple cells for scaling
            scaler = StandardScaler()
            chunk_expr_norm = scaler.fit_transform(chunk_expr)
            chunk_scores = np.mean(chunk_expr_norm, axis=1)
        else:
            chunk_scores = np.mean(chunk_expr, axis=1)
        
        platelet_scores[start_idx:end_idx] = chunk_scores
        
        # Clean up
        del chunk_data, chunk_expr
        gc.collect()
    
    print("Identifying platelets...")
    
    # Select top 10% as platelets (more selective than macrophages)
    threshold = np.percentile(platelet_scores, 90)
    is_platelet = platelet_scores >= threshold
    
    print(f"Platelet threshold: {threshold:.3f}")
    print(f"Cells above threshold: {np.sum(is_platelet):,}")
    
    # Quality filter based on total counts
    print("Applying quality filter...")
    
    if hasattr(adata.X, 'toarray'):
        # For sparse matrices, calculate row sums efficiently
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    else:
        total_counts = adata.X.sum(axis=1)
    
    min_counts = np.percentile(total_counts, 10)
    high_quality = total_counts >= min_counts
    
    print(f"Min count threshold: {min_counts:.0f}")
    print(f"High quality cells: {np.sum(high_quality):,}")
    
    # Combine filters
    final_mask = is_platelet & high_quality
    n_platelets = np.sum(final_mask)
    
    print(f"Final platelets identified: {n_platelets:,}")
    
    if n_platelets == 0:
        print("❌ No platelets identified")
        return
    
    # Extract platelets
    print("Extracting platelet subset...")
    platelets = adata[final_mask].copy()
    
    # Add metadata
    platelets.obs['cell_type'] = 'platelet'
    platelets.obs['platelet_score'] = platelet_scores[final_mask]
    platelets.obs['source_dataset'] = 'HEK293T'
    
    # Clean up original data
    del adata
    gc.collect()
    
    print("Saving platelets...")
    
    # Save to specified location
    output_path = "/bigdata/godziklab/shared/Xinru/302006/data/raw/hek293t_platelets.h5ad"
    platelets.write_h5ad(output_path)
    
    print(f"✅ Successfully extracted {platelets.n_obs:,} platelets from HEK293T")
    print(f"✅ Saved to: {output_path}")
    print(f"Platelet score range: {platelets.obs['platelet_score'].min():.3f} - {platelets.obs['platelet_score'].max():.3f}")
    
    # Show perturbation info if available
    if 'guide_target' in platelets.obs.columns:
        print("\nPerturbation distribution:")
        print(platelets.obs['guide_target'].value_counts().head(10))

if __name__ == "__main__":
    extract_platelets_chunked()