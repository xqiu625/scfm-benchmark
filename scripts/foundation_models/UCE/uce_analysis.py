#!/usr/bin/env python3
"""
UCE Embedding Results Analysis Script
Comprehensive analysis of UCE embeddings across all 6 datasets
Foundation Model Comparison Project - Phase 1 Analysis
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set up scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

class UCEAnalyzer:
    def __init__(self, base_path="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_path = Path(base_path)
        self.embeddings_path = self.base_path / "embeddings" / "uce"
        self.results_path = self.base_path / "results" / "uce_analysis"
        self.results_path.mkdir(parents=True, exist_ok=True)
        
        # Dataset information with actual file mapping and corrected cell counts
        self.datasets = {
            'hek293t_epithelial': {
                'file': 'hek293t_epithelial_uce_adata.h5ad', 
                'cells': 180768, 
                'type': 'epithelial'
            },
            'hct116_epithelial': {
                'file': 'hct116_epithelial_uce_adata.h5ad',
                'cells': 183797,  # Updated from inspection
                'type': 'epithelial'
            }, 
            'mouse_epithelial_retnla': {
                'file': 'mouse_epithelial_retnla_ko_wt_uce_adata.h5ad',
                'cells': 131,  # Much smaller than expected
                'type': 'epithelial'
            },
            'hek293t_platelets': {
                'file': 'hek293t_platelets_uce_adata.h5ad',
                'cells': 88434,  # Updated from inspection
                'type': 'platelets'
            },
            'hct116_platelets': {
                'file': 'hct116_platelets_uce_adata.h5ad',
                'cells': 89738,  # Updated from inspection
                'type': 'platelets'
            },
            'covid_platelets': {
                'file': 'human_platelet_covid_severity_uce_adata.h5ad',
                'cells': 47977,  # Updated from inspection
                'type': 'platelets'
            }
        }
        
    def load_uce_embeddings(self):
        """Load all UCE embedding files"""
        print("🔍 Loading UCE embeddings for all datasets...")
        self.adata_dict = {}
        
        for dataset_name, dataset_info in self.datasets.items():
            embedding_file = self.embeddings_path / dataset_info['file']
            
            if embedding_file.exists():
                print(f"  ✅ Loading {dataset_name}...")
                adata = sc.read_h5ad(embedding_file)
                print(f"     Shape: {adata.shape}")
                
                # Check for UCE embeddings in different possible locations
                uce_key = None
                for key in ['X_uce', 'X_UCE', 'uce', 'UCE']:
                    if key in adata.obsm.keys():
                        uce_key = key
                        break
                
                if uce_key:
                    print(f"     UCE embeddings: {adata.obsm[uce_key].shape} (key: {uce_key})")
                    # Standardize the key name
                    if uce_key != 'X_uce':
                        adata.obsm['X_uce'] = adata.obsm[uce_key].copy()
                else:
                    print(f"     ⚠️ No UCE embeddings found. Available keys: {list(adata.obsm.keys())}")
                    # Check if embeddings are in .X
                    if hasattr(adata.X, 'shape') and len(adata.X.shape) == 2:
                        print(f"     🔍 Using .X as embeddings: {adata.X.shape}")
                        adata.obsm['X_uce'] = adata.X.copy()
                
                self.adata_dict[dataset_name] = adata
            else:
                print(f"  ❌ Missing: {embedding_file}")
                
        print(f"\n📊 Successfully loaded {len(self.adata_dict)} datasets")
        return self.adata_dict
    
    def analyze_embedding_quality(self):
        """Analyze UCE embedding quality metrics"""
        print("\n🔬 Analyzing UCE embedding quality...")
        
        quality_results = {}
        
        for name, adata in self.adata_dict.items():
            print(f"\n📈 Analyzing {name}...")
            
            # Get UCE embeddings
            uce_embeddings = adata.obsm['X_uce']
            
            # Basic statistics
            stats = {
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'embedding_dim': uce_embeddings.shape[1],
                'embedding_mean': np.mean(uce_embeddings),
                'embedding_std': np.std(uce_embeddings),
                'embedding_range': (np.min(uce_embeddings), np.max(uce_embeddings))
            }
            
            # Check for any NaN or infinite values
            stats['has_nan'] = np.isnan(uce_embeddings).any()
            stats['has_inf'] = np.isinf(uce_embeddings).any()
            
            # Embedding variance analysis
            embedding_var = np.var(uce_embeddings, axis=0)
            stats['min_variance'] = np.min(embedding_var)
            stats['max_variance'] = np.max(embedding_var)
            stats['mean_variance'] = np.mean(embedding_var)
            
            quality_results[name] = stats
            
            print(f"  📏 Cells: {stats['n_cells']:,}")
            print(f"  🧬 Genes: {stats['n_genes']:,}")
            print(f"  🔢 Embedding dim: {stats['embedding_dim']}")
            print(f"  📊 Mean±SD: {stats['embedding_mean']:.3f}±{stats['embedding_std']:.3f}")
            print(f"  🎯 Range: ({stats['embedding_range'][0]:.3f}, {stats['embedding_range'][1]:.3f})")
            print(f"  ⚠️ Issues: NaN={stats['has_nan']}, Inf={stats['has_inf']}")
        
        return quality_results
    
    def perform_dimensionality_reduction(self):
        """Apply UMAP and PCA to UCE embeddings - FIXED CONNECTIVITY"""
        print("\n🗺️ Performing dimensionality reduction (UMAP & PCA)...")
        
        for name, adata in self.adata_dict.items():
            print(f"  🔄 Processing {name}...")
            
            try:
                # Get UCE embeddings directly
                uce_embeddings = adata.obsm['X_uce']
                print(f"    📊 UCE shape: {uce_embeddings.shape}")
                
                # Manual PCA using sklearn
                from sklearn.decomposition import PCA
                n_comps = min(50, uce_embeddings.shape[1] - 1, uce_embeddings.shape[0] - 1)
                pca = PCA(n_components=n_comps)
                pca_result = pca.fit_transform(uce_embeddings)
                adata.obsm['X_pca'] = pca_result
                adata.uns['pca'] = {'variance_ratio': pca.explained_variance_ratio_}
                
                # Manual UMAP using umap-learn
                import umap
                n_neighbors = min(15, uce_embeddings.shape[0] - 1)
                reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.5, spread=1.0, random_state=42)
                umap_result = reducer.fit_transform(uce_embeddings)
                adata.obsm['X_umap'] = umap_result
                
                # Create neighbors graph for clustering (simplified)
                from sklearn.neighbors import NearestNeighbors
                from scipy.sparse import csr_matrix
                
                nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
                nn.fit(uce_embeddings)
                distances, indices = nn.kneighbors(uce_embeddings)
                
                # Create proper connectivity matrix (cell x cell)
                n_cells = uce_embeddings.shape[0]
                row_ind = np.repeat(np.arange(n_cells), n_neighbors)
                col_ind = indices.ravel()
                data = np.ones(len(row_ind))
                connectivities = csr_matrix((data, (row_ind, col_ind)), shape=(n_cells, n_cells))
                connectivities = connectivities.maximum(connectivities.T)  # Make symmetric
                
                # Create distance matrix (sparse, only store k-nearest neighbors)
                distances_data = distances.ravel()
                distances_matrix = csr_matrix((distances_data, (row_ind, col_ind)), shape=(n_cells, n_cells))
                
                # Store in AnnData format
                adata.obsp['connectivities'] = connectivities
                adata.obsp['distances'] = distances_matrix
                
                # Store neighbors info in scanpy format
                adata.uns['neighbors'] = {
                    'connectivities_key': 'connectivities',
                    'distances_key': 'distances',
                    'params': {
                        'n_neighbors': n_neighbors,
                        'use_rep': 'X_uce',
                        'metric': 'euclidean'
                    }
                }
                
                print(f"    ✅ PCA: {adata.obsm['X_pca'].shape}")
                print(f"    ✅ UMAP: {adata.obsm['X_umap'].shape}")
                print(f"    ✅ Neighbors: {n_neighbors} neighbors, connectivity: {connectivities.shape}")
                
            except Exception as e:
                print(f"    ❌ Error processing {name}: {e}")
                # Minimal fallback - just PCA and UMAP without neighbors
                try:
                    uce_embeddings = adata.obsm['X_uce']
                    
                    # PCA
                    from sklearn.decomposition import PCA
                    n_comps = min(50, uce_embeddings.shape[1] - 1, uce_embeddings.shape[0] - 1)
                    pca = PCA(n_components=n_comps)
                    adata.obsm['X_pca'] = pca.fit_transform(uce_embeddings)
                    
                    # UMAP
                    import umap
                    n_neighbors = min(15, uce_embeddings.shape[0] - 1)
                    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.5, spread=1.0, random_state=42)
                    adata.obsm['X_umap'] = reducer.fit_transform(uce_embeddings)
                    
                    print(f"    ⚠️ Fallback successful: PCA {adata.obsm['X_pca'].shape}, UMAP {adata.obsm['X_umap'].shape}")
                    
                except Exception as e2:
                    print(f"    ❌ Fallback failed: {e2}")
                    continue
    
    def analyze_cell_type_clustering(self):
        """Analyze cell type separation in UCE embedding space - FIXED"""
        print("\n🎯 Analyzing cell type clustering...")
        
        clustering_results = {}
        
        for name, adata in self.adata_dict.items():
            print(f"  📊 Clustering analysis for {name}...")
            
            try:
                # Only proceed if neighbors graph exists
                if 'neighbors' in adata.uns:
                    # Perform Leiden clustering
                    sc.tl.leiden(adata, resolution=0.5, key_added='leiden_uce')
                    
                    # Calculate silhouette scores if cell_type annotation exists
                    if 'cell_type' in adata.obs.columns:
                        from sklearn.metrics import silhouette_score, adjusted_rand_score
                        
                        # Silhouette score using cell types
                        sil_score = silhouette_score(adata.obsm['X_uce'], adata.obs['cell_type'])
                        
                        # ARI between leiden clustering and cell types
                        ari_score = adjusted_rand_score(adata.obs['leiden_uce'], adata.obs['cell_type'])
                        
                        clustering_results[name] = {
                            'silhouette_score': sil_score,
                            'ari_score': ari_score,
                            'n_leiden_clusters': len(adata.obs['leiden_uce'].unique()),
                            'n_cell_types': len(adata.obs['cell_type'].unique())
                        }
                        
                        print(f"    🎯 Silhouette score: {sil_score:.3f}")
                        print(f"    🔗 ARI score: {ari_score:.3f}")
                        print(f"    📦 Leiden clusters: {clustering_results[name]['n_leiden_clusters']}")
                        print(f"    🧬 Cell types: {clustering_results[name]['n_cell_types']}")
                    else:
                        print(f"    ⚠️ No cell_type annotation found")
                        clustering_results[name] = {'note': 'No cell_type annotation'}
                else:
                    print(f"    ⚠️ No neighbors graph found, skipping clustering")
                    clustering_results[name] = {'note': 'No neighbors graph'}
            except Exception as e:
                print(f"    ❌ Error in clustering analysis for {name}: {e}")
                clustering_results[name] = {'error': str(e)}
        
        return clustering_results
    
    def generate_visualization_plots(self):
        """Generate comprehensive visualization plots - WITH ERROR HANDLING"""
        print("\n🎨 Generating visualization plots...")
        
        # Create figure directory
        fig_dir = self.results_path / "figures"
        fig_dir.mkdir(exist_ok=True)
        
        # Check which datasets have UMAP
        valid_datasets = []
        for name, adata in self.adata_dict.items():
            if 'X_umap' in adata.obsm:
                valid_datasets.append((name, adata))
            else:
                print(f"    ⚠️ Skipping {name} - no UMAP coordinates")
        
        if not valid_datasets:
            print("    ❌ No datasets with UMAP coordinates found")
            return
            
        # 1. UMAP plots for valid datasets
        n_datasets = len(valid_datasets)
        if n_datasets > 0:
            cols = 3
            rows = (n_datasets + cols - 1) // cols
            
            fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 4*rows))
            if rows == 1 and cols == 1:
                axes = [axes]
            elif rows == 1 or cols == 1:
                axes = axes.flatten()
            else:
                axes = axes.flatten()
            
            for i, (name, adata) in enumerate(valid_datasets):
                ax = axes[i] if i < len(axes) else None
                if ax is None:
                    continue
                    
                try:
                    umap_coords = adata.obsm['X_umap']
                    
                    # Color by cell type if available
                    if 'cell_type' in adata.obs.columns:
                        cell_types = adata.obs['cell_type']
                        unique_types = cell_types.unique()
                        
                        if len(unique_types) > 1:
                            # Multiple cell types - use colors
                            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_types)))
                            for j, cell_type in enumerate(unique_types):
                                if pd.isna(cell_type):
                                    continue
                                mask = cell_types == cell_type
                                ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1], 
                                          c=[colors[j]], label=str(cell_type), alpha=0.6, s=1)
                            if len(unique_types) <= 10:
                                ax.legend(markerscale=5, fontsize=8)
                        else:
                            # Single cell type
                            ax.scatter(umap_coords[:, 0], umap_coords[:, 1], alpha=0.6, s=1, c='blue')
                    else:
                        # No cell type info
                        ax.scatter(umap_coords[:, 0], umap_coords[:, 1], alpha=0.6, s=1, c='gray')
                    
                    ax.set_title(f"{name.replace('_', ' ').title()}\n({adata.n_obs:,} cells)", fontsize=10)
                    ax.set_xlabel("UMAP1")
                    ax.set_ylabel("UMAP2")
                    
                except Exception as e:
                    print(f"    ❌ Error plotting {name}: {e}")
                    ax.text(0.5, 0.5, f"Error plotting\n{name}", ha='center', va='center', 
                           transform=ax.transAxes)
            
            # Hide extra axes
            for i in range(len(valid_datasets), len(axes)):
                axes[i].set_visible(False)
            
            plt.tight_layout()
            plt.savefig(fig_dir / "uce_umap_overview.png", dpi=300, bbox_inches='tight')
            plt.close()
            print(f"    ✅ UMAP overview saved")
        
        # 2. Embedding distribution plots
        try:
            fig, axes = plt.subplots(2, 3, figsize=(18, 10))
            axes = axes.flatten()
            
            for i, (name, adata) in enumerate(self.adata_dict.items()):
                if i >= len(axes):
                    break
                    
                ax = axes[i]
                embeddings = adata.obsm['X_uce']
                
                # Plot first few dimensions
                for dim in range(min(5, embeddings.shape[1])):
                    ax.hist(embeddings[:, dim], alpha=0.6, bins=50, 
                           label=f'Dim {dim+1}')
                
                ax.set_title(f"{name.replace('_', ' ').title()}")
                ax.set_xlabel("Embedding Value")
                ax.set_ylabel("Frequency")
                if dim < 4:  # Only show legend if few dimensions
                    ax.legend()
            
            plt.tight_layout()
            plt.savefig(fig_dir / "uce_embedding_distributions.png", dpi=300, bbox_inches='tight')
            plt.close()
            print(f"    ✅ Embedding distributions saved")
        
        except Exception as e:
            print(f"    ❌ Error generating distribution plots: {e}")
        
        print(f"    📁 Plots saved to {fig_dir}")
    
    def covid_progression_analysis(self):
        """Enhanced analysis for COVID platelet progression"""
        print("\n🦠 COVID progression analysis...")
        
        if 'covid_platelets' in self.adata_dict:
            adata = self.adata_dict['covid_platelets']
            
            print(f"    📊 COVID dataset: {adata.shape}")
            
            # Analyze COVID severity progression
            if 'covid_severity' in adata.obs.columns:
                severity_counts = adata.obs['covid_severity'].value_counts()
                print(f"    🔍 COVID severity distribution:")
                for severity, count in severity_counts.items():
                    print(f"      {severity}: {count:,} cells")
                
                # Check for progression group
                if 'progression_group' in adata.obs.columns:
                    prog_counts = adata.obs['progression_group'].value_counts()
                    print(f"    📈 Progression groups:")
                    for group, count in prog_counts.items():
                        print(f"      {group}: {count:,} cells")
                
                # Analyze embedding space separation by severity
                from sklearn.metrics import silhouette_score
                
                # Remove NaN values for analysis
                valid_mask = ~adata.obs['covid_severity'].isna()
                if valid_mask.sum() > 0:
                    valid_embeddings = adata.obsm['X_uce'][valid_mask]
                    valid_severity = adata.obs['covid_severity'][valid_mask]
                    
                    if len(valid_severity.unique()) > 1:
                        sil_score = silhouette_score(valid_embeddings, valid_severity)
                        print(f"    🎯 Severity silhouette score: {sil_score:.3f}")
                        
                        # Create severity ordering for trajectory analysis
                        severity_order = ['healthy', 'mild', 'moderate', 'severe', 'fatal']
                        severity_in_data = [s for s in severity_order if s in valid_severity.unique()]
                        print(f"    📋 Severity progression available: {' → '.join(severity_in_data)}")
            
            # Analyze other COVID-related columns
            covid_cols = ['Disease', 'Disease_group', 'Death', 'Category', 'mild_vs_severe']
            for col in covid_cols:
                if col in adata.obs.columns:
                    unique_vals = adata.obs[col].dropna().unique()
                    print(f"    📊 {col}: {list(unique_vals)}")
                    
        else:
            print("    ❌ COVID platelets dataset not loaded")
    
    def mouse_genotype_analysis(self):
        """Analyze mouse KO/WT genotype differences"""
        print("\n🐭 Mouse genotype analysis...")
        
        if 'mouse_epithelial_retnla' in self.adata_dict:
            adata = self.adata_dict['mouse_epithelial_retnla']
            
            print(f"    📊 Mouse dataset: {adata.shape}")
            
            if 'genotype' in adata.obs.columns:
                geno_counts = adata.obs['genotype'].value_counts()
                print(f"    🧬 Genotype distribution:")
                for geno, count in geno_counts.items():
                    print(f"      {geno}: {count} cells")
                
                # Analyze condition combinations
                if 'condition' in adata.obs.columns:
                    condition_counts = adata.obs['condition'].value_counts()
                    print(f"    🔬 Experimental conditions:")
                    for cond, count in condition_counts.items():
                        print(f"      {cond}: {count} cells")
                
                # Embedding space separation by genotype
                from sklearn.metrics import silhouette_score
                embeddings = adata.obsm['X_uce']
                genotypes = adata.obs['genotype']
                
                if len(genotypes.unique()) > 1:
                    sil_score = silhouette_score(embeddings, genotypes)
                    print(f"    🎯 Genotype silhouette score: {sil_score:.3f}")
        else:
            print("    ❌ Mouse dataset not loaded")
    
    def generate_summary_report(self, quality_results, clustering_results):
        """Generate comprehensive summary report"""
        print("\n📝 Generating summary report...")
        
        report_file = self.results_path / "uce_analysis_report.md"
        
        with open(report_file, 'w') as f:
            f.write("# UCE Embedding Analysis Report\n\n")
            f.write("## Dataset Overview\n\n")
            
            # Quality results table
            f.write("### Embedding Quality Metrics\n\n")
            f.write("| Dataset | Cells | Genes | Embedding Dim | Mean±SD | Range |\n")
            f.write("|---------|-------|-------|---------------|---------|-------|\n")
            
            for name, stats in quality_results.items():
                f.write(f"| {name} | {stats['n_cells']:,} | {stats['n_genes']:,} | "
                       f"{stats['embedding_dim']} | {stats['embedding_mean']:.3f}±{stats['embedding_std']:.3f} | "
                       f"({stats['embedding_range'][0]:.3f}, {stats['embedding_range'][1]:.3f}) |\n")
            
            # Clustering results
            f.write("\n### Clustering Performance\n\n")
            f.write("| Dataset | Silhouette Score | ARI Score | Leiden Clusters | Cell Types |\n")
            f.write("|---------|------------------|-----------|-----------------|------------|\n")
            
            for name, results in clustering_results.items():
                if 'silhouette_score' in results:
                    f.write(f"| {name} | {results['silhouette_score']:.3f} | "
                           f"{results['ari_score']:.3f} | {results['n_leiden_clusters']} | "
                           f"{results['n_cell_types']} |\n")
                else:
                    f.write(f"| {name} | N/A | N/A | N/A | N/A |\n")
            
            f.write(f"\n## Analysis Completed: {pd.Timestamp.now()}\n")
            f.write(f"Results saved to: {self.results_path}\n")
        
        print(f"    ✅ Report saved: {report_file}")
    
    def run_complete_analysis(self):
        """Run the complete UCE analysis pipeline"""
        print("🚀 Starting comprehensive UCE embedding analysis...")
        print("="*60)
        
        # Load embeddings
        self.load_uce_embeddings()
        
        if not self.adata_dict:
            print("❌ No UCE embeddings found. Please check file locations.")
            return
        
        # Quality analysis
        quality_results = self.analyze_embedding_quality()
        
        # Dimensionality reduction
        self.perform_dimensionality_reduction()
        
        # Clustering analysis
        clustering_results = self.analyze_cell_type_clustering()
        
        # Generate visualizations
        self.generate_visualization_plots()
        
        # COVID specific analysis
        self.covid_progression_analysis()
        
        # Mouse genotype analysis
        self.mouse_genotype_analysis()
        
        # Generate report
        self.generate_summary_report(quality_results, clustering_results)
        
        print("\n" + "="*60)
        print("✅ UCE embedding analysis complete!")
        print(f"📁 Results saved to: {self.results_path}")
        print("🎯 Ready to proceed with scGPT processing completion")


if __name__ == "__main__":
    # Initialize analyzer
    analyzer = UCEAnalyzer()
    
    # Run complete analysis
    analyzer.run_complete_analysis()
    
    print("\n🎉 UCE Analysis Complete!")
    print("📋 Next steps:")
    print("   1. Review UCE results in /results/uce_analysis/")
    print("   2. Complete remaining scGPT chunks (26/31)")
    print("   3. Process remaining 5 datasets with scGPT")
    print("   4. Begin Geneformer integration")
    print("   5. Three-way foundation model comparison")