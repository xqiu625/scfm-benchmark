#!/usr/bin/env python3
# Extract HCT116 RETNLB perturbations and identify cell types

import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Load HCT116 processed data
print("Loading HCT116 processed data...")
adata = sc.read_h5ad("/bigdata/godziklab/shared/Xinru/302006/data/processed/HCT116_processed.h5ad")
print(f"Loaded: {adata.shape}")

# Check metadata columns
print(f"\nMetadata columns: {list(adata.obs.columns)}")

# Look for perturbation/guide columns
perturbation_cols = []
for col in adata.obs.columns:
    if any(word in col.lower() for word in ['guide', 'target', 'perturbation', 'treatment']):
        perturbation_cols.append(col)

print(f"Potential perturbation columns: {perturbation_cols}")

# Check each perturbation column for RETNLB
retnlb_found = False
perturbation_col = None

for col in perturbation_cols:
    unique_values = adata.obs[col].value_counts()
    print(f"\n{col} unique values:")
    for val, count in unique_values.head(20).items():
        print(f"  {val}: {count:,} cells")
    
    # Check for RETNLB
    retnlb_mask = adata.obs[col].str.contains('RETNLB', case=False, na=False)
    if retnlb_mask.sum() > 0:
        print(f"🎯 Found RETNLB in {col}: {retnlb_mask.sum():,} cells")
        retnlb_found = True
        perturbation_col = col
        break

if not retnlb_found:
    print("❌ RETNLB not found in obvious columns. Checking all columns...")
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            retnlb_mask = adata.obs[col].str.contains('RETNLB', case=False, na=False)
            if retnlb_mask.sum() > 0:
                print(f"🎯 Found RETNLB in {col}: {retnlb_mask.sum():,} cells")
                perturbation_col = col
                retnlb_found = True
                break

if retnlb_found:
    # Extract RETNLB perturbations
    print(f"\n📊 Extracting RETNLB perturbations from {perturbation_col}...")
    
    retnlb_mask = adata.obs[perturbation_col].str.contains('RETNLB', case=False, na=False)
    retnlb_cells = adata[retnlb_mask].copy()
    
    print(f"✅ Extracted {retnlb_cells.n_obs:,} RETNLB perturbed cells")
    
    # Show RETNLB perturbation details
    print(f"\nRETNLB perturbation details:")
    retnlb_values = retnlb_cells.obs[perturbation_col].value_counts()
    for val, count in retnlb_values.items():
        print(f"  {val}: {count:,} cells")
    
    # Identify cell types using multiple marker sets
    print(f"\n🔬 IDENTIFYING CELL TYPES...")
    
    # Define cell type markers
    cell_type_markers = {
        'macrophage': ['CD68', 'CD14', 'CSF1R', 'ADGRE1', 'CD163', 'MRC1', 'FCGR1A'],
        'platelet': ['ITGA2B', 'ITGB3', 'GP1BA', 'GP1BB', 'THBS1', 'PDGFA', 'PF4', 'PPBP'],
        't_cell': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B'],
        'b_cell': ['CD19', 'CD20', 'MS4A1', 'CD79A', 'CD79B'],
        'nk_cell': ['KLRB1', 'KLRD1', 'NCAM1', 'NKG7'],
        'neutrophil': ['ELANE', 'MPO', 'DEFA1', 'DEFA4', 'LTF'],
        'monocyte': ['CD14', 'FCGR3A', 'LYZ', 'S100A8', 'S100A9'],
        'epithelial': ['EPCAM', 'KRT18', 'KRT19', 'CDH1'],
        'endothelial': ['PECAM1', 'VWF', 'CDH5', 'FLT1'],
        'fibroblast': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM']
    }
    
    # Calculate scores for each cell type
    cell_type_scores = {}
    
    for cell_type, markers in cell_type_markers.items():
        available_markers = [m for m in markers if m in retnlb_cells.var_names]
        
        if len(available_markers) >= 2:  # Need at least 2 markers
            print(f"\n{cell_type.upper()}: {available_markers}")
            
            # Calculate cell type score
            marker_expression = []
            for marker in available_markers:
                if hasattr(retnlb_cells.X, 'toarray'):
                    expr = retnlb_cells[:, marker].X.toarray().flatten()
                else:
                    expr = retnlb_cells[:, marker].X.flatten()
                marker_expression.append(expr)
            
            marker_matrix = np.column_stack(marker_expression)
            scaler = StandardScaler()
            marker_matrix_norm = scaler.fit_transform(marker_matrix)
            cell_type_score = np.mean(marker_matrix_norm, axis=1)
            
            cell_type_scores[cell_type] = cell_type_score
            
            # Add to metadata
            retnlb_cells.obs[f'{cell_type}_score'] = cell_type_score
            
            # Show statistics
            print(f"  Score range: {cell_type_score.min():.3f} to {cell_type_score.max():.3f}")
            print(f"  Mean score: {cell_type_score.mean():.3f}")
            
        else:
            print(f"\n{cell_type.upper()}: Only {len(available_markers)} markers available - skipping")
    
    # Determine dominant cell type for each cell
    if cell_type_scores:
        print(f"\n🏆 DETERMINING DOMINANT CELL TYPES...")
        
        # Create matrix of all scores
        score_matrix = np.column_stack([cell_type_scores[ct] for ct in cell_type_scores.keys()])
        cell_types_list = list(cell_type_scores.keys())
        
        # Find highest scoring cell type for each cell
        max_indices = np.argmax(score_matrix, axis=1)
        predicted_cell_types = [cell_types_list[i] for i in max_indices]
        
        # Add to metadata
        retnlb_cells.obs['predicted_cell_type'] = predicted_cell_types
        
        # Show distribution
        cell_type_distribution = pd.Series(predicted_cell_types).value_counts()
        print(f"Predicted cell type distribution:")
        for cell_type, count in cell_type_distribution.items():
            percentage = (count / len(predicted_cell_types)) * 100
            print(f"  {cell_type}: {count:,} cells ({percentage:.1f}%)")
        
        # Find cells with high confidence (score > threshold)
        print(f"\n🎯 HIGH CONFIDENCE CELL TYPE ASSIGNMENTS...")
        
        confidence_threshold = 1.0  # Above average + 1 std
        
        for i, cell_type in enumerate(cell_types_list):
            high_conf_mask = score_matrix[:, i] > confidence_threshold
            high_conf_count = np.sum(high_conf_mask)
            
            if high_conf_count > 0:
                percentage = (high_conf_count / len(retnlb_cells)) * 100
                print(f"  {cell_type}: {high_conf_count:,} high-confidence cells ({percentage:.1f}%)")
                
                # Extract high-confidence cells of this type
                if high_conf_count > 100:  # Only save if reasonable number
                    high_conf_cells = retnlb_cells[high_conf_mask].copy()
                    high_conf_cells.obs['confidence_level'] = 'high'
                    
                    # Save high-confidence cell type
                    output_file = f"/bigdata/godziklab/shared/Xinru/302006/data/raw/hct116_retnlb_{cell_type}_high_conf.h5ad"
                    high_conf_cells.write_h5ad(output_file)
                    print(f"    💾 Saved: {output_file}")
    
    # Save all RETNLB cells with annotations
    retnlb_cells.obs['dataset'] = 'HCT116'
    retnlb_cells.obs['perturbation'] = 'RETNLB'
    
    output_file = "/bigdata/godziklab/shared/Xinru/302006/data/raw/hct116_retnlb_all_annotated.h5ad"
    retnlb_cells.write_h5ad(output_file)
    print(f"\n💾 Saved all annotated RETNLB cells: {output_file}")
    
    # Summary report
    print(f"\n📋 SUMMARY REPORT:")
    print(f"  Total RETNLB cells: {retnlb_cells.n_obs:,}")
    print(f"  Genes analyzed: {retnlb_cells.n_vars:,}")
    print(f"  Cell types identified: {len(cell_type_scores)}")
    
    if cell_type_scores:
        dominant_type = cell_type_distribution.index[0]
        dominant_count = cell_type_distribution.iloc[0]
        dominant_pct = (dominant_count / len(retnlb_cells)) * 100
        print(f"  Dominant cell type: {dominant_type} ({dominant_count:,} cells, {dominant_pct:.1f}%)")

else:
    print("❌ RETNLB perturbations not found in any column")
    print("Available perturbation values in first column:")
    if perturbation_cols:
        print(adata.obs[perturbation_cols[0]].value_counts().head(10))