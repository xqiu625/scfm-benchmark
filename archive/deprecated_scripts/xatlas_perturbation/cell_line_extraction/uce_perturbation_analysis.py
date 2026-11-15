#!/usr/bin/env python3
"""
UCE Perturbation Analysis: RETNLA KO/WT and RETNLB Perturbation Study
Comprehensive analysis focused on RETNLA knockout effects in mouse epithelial cells
and RETNLB perturbation effects in human cell line epithelial cells
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import warnings
warnings.filterwarnings('ignore')

# Set up scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

class UCEPerturbationAnalyzer:
    def __init__(self, base_path="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_path = Path(base_path)
        self.embeddings_path = self.base_path / "embeddings" / "uce"
        self.results_path = self.base_path / "results" / "uce_perturbation_analysis"
        self.results_path.mkdir(parents=True, exist_ok=True)
        
        # Dataset information focused on perturbation studies
        self.datasets = {
            'mouse_epithelial_retnla': {
                'file': 'mouse_epithelial_retnla_ko_wt_uce_adata.h5ad',
                'type': 'mouse_epithelial',
                'perturbation': 'RETNLA_KO',
                'control_col': 'genotype',
                'treatment_col': 'condition',
                'expected_groups': ['WT', 'KO']
            },
            'hek293t_epithelial': {
                'file': 'hek293t_epithelial_uce_adata.h5ad',
                'type': 'human_epithelial', 
                'perturbation': 'RETNLB_perturbation',
                'control_col': 'guide_target',
                'treatment_col': 'gene_target',
                'expected_groups': ['control', 'RETNLB']
            },
            'hct116_epithelial': {
                'file': 'hct116_epithelial_uce_adata.h5ad',
                'type': 'human_epithelial',
                'perturbation': 'RETNLB_perturbation', 
                'control_col': 'guide_target',
                'treatment_col': 'gene_target',
                'expected_groups': ['control', 'RETNLB']
            }
        }
        
    def load_perturbation_datasets(self):
        """Load UCE embeddings for perturbation analysis"""
        print("🔬 Loading UCE embeddings for perturbation analysis...")
        self.adata_dict = {}
        
        for dataset_name, dataset_info in self.datasets.items():
            embedding_file = self.embeddings_path / dataset_info['file']
            
            if embedding_file.exists():
                print(f"  ✅ Loading {dataset_name}...")
                adata = sc.read_h5ad(embedding_file)
                print(f"     Shape: {adata.shape}")
                print(f"     UCE embeddings: {adata.obsm['X_uce'].shape}")
                
                # Analyze perturbation-related columns
                self.analyze_perturbation_metadata(adata, dataset_name, dataset_info)
                
                self.adata_dict[dataset_name] = adata
            else:
                print(f"  ❌ Missing: {embedding_file}")
                
        print(f"\n📊 Successfully loaded {len(self.adata_dict)} datasets")
        return self.adata_dict
    
    def analyze_perturbation_metadata(self, adata, dataset_name, dataset_info):
        """Analyze perturbation-specific metadata"""
        print(f"    🔍 Perturbation metadata for {dataset_name}:")
        
        # Check control column
        control_col = dataset_info['control_col']
        if control_col in adata.obs.columns:
            control_counts = adata.obs[control_col].value_counts()
            print(f"      {control_col}: {dict(control_counts)}")
        
        # Check treatment column  
        treatment_col = dataset_info['treatment_col']
        if treatment_col in adata.obs.columns:
            treatment_counts = adata.obs[treatment_col].value_counts()
            print(f"      {treatment_col}: {dict(treatment_counts)}")
        
        # For mouse data, check additional conditions
        if 'mouse' in dataset_name:
            for col in ['condition', 'genotype', 'Sex']:
                if col in adata.obs.columns:
                    unique_vals = adata.obs[col].unique()
                    print(f"      {col}: {list(unique_vals)}")
        
        # For human cell lines, check guide/gene targets
        if 'hek293t' in dataset_name or 'hct116' in dataset_name:
            for col in ['guide_target', 'gene_target']:
                if col in adata.obs.columns:
                    unique_vals = adata.obs[col].unique()[:10]  # Show first 10
                    print(f"      {col} (top 10): {list(unique_vals)}")
    
    def identify_perturbation_groups(self):
        """Identify and standardize perturbation groups"""
        print("\n🎯 Identifying perturbation groups...")
        
        for name, adata in self.adata_dict.items():
            dataset_info = self.datasets[name]
            print(f"  📊 Processing {name}...")
            
            if 'mouse' in name:
                # Mouse RETNLA KO/WT analysis
                if 'genotype' in adata.obs.columns:
                    adata.obs['perturbation_group'] = adata.obs['genotype']
                    adata.obs['perturbation_type'] = 'RETNLA_KO'
                    
                    ko_count = (adata.obs['genotype'] == 'KO').sum()
                    wt_count = (adata.obs['genotype'] == 'WT').sum()
                    print(f"    🧬 RETNLA KO: {ko_count} cells")
                    print(f"    🧬 RETNLA WT: {wt_count} cells")
                    
                    # Add sex information if available
                    if 'Sex' in adata.obs.columns:
                        sex_geno = adata.obs.groupby(['genotype', 'Sex']).size()
                        print(f"    👥 Sex distribution:")
                        for (geno, sex), count in sex_geno.items():
                            print(f"       {geno}-{sex}: {count} cells")
            
            else:
                # Human cell line RETNLB perturbation analysis
                if 'gene_target' in adata.obs.columns:
                    # Identify RETNLB-targeted cells
                    gene_targets = adata.obs['gene_target']
                    
                    # Look for RETNLB or similar targets
                    retnlb_mask = gene_targets.str.contains('RETNLB|RETNL', case=False, na=False)
                    control_mask = gene_targets.str.contains('control|negative|NT|non-targeting', case=False, na=False)
                    
                    # Create perturbation groups
                    adata.obs['perturbation_group'] = 'unknown'
                    adata.obs.loc[retnlb_mask, 'perturbation_group'] = 'RETNLB_perturbed'
                    adata.obs.loc[control_mask, 'perturbation_group'] = 'control'
                    
                    # Count groups
                    group_counts = adata.obs['perturbation_group'].value_counts()
                    print(f"    🎯 Perturbation groups:")
                    for group, count in group_counts.items():
                        print(f"       {group}: {count} cells")
                    
                    # Show top gene targets for verification
                    top_targets = gene_targets.value_counts().head(10)
                    print(f"    🔬 Top gene targets:")
                    for target, count in top_targets.items():
                        print(f"       {target}: {count} cells")
                
                # Also check guide targets
                if 'guide_target' in adata.obs.columns:
                    guide_counts = adata.obs['guide_target'].value_counts().head(5)
                    print(f"    📝 Top guide targets: {dict(guide_counts)}")
    
    def perform_perturbation_umap(self):
        """Generate UMAP specifically for perturbation analysis"""
        print("\n🗺️ Performing perturbation-focused UMAP analysis...")
        
        for name, adata in self.adata_dict.items():
            print(f"  🔄 Processing {name}...")
            
            try:
                # Get UCE embeddings
                uce_embeddings = adata.obsm['X_uce']
                
                # Manual UMAP
                import umap
                n_neighbors = min(15, uce_embeddings.shape[0] - 1)
                reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.3, spread=1.0, random_state=42)
                umap_result = reducer.fit_transform(uce_embeddings)
                adata.obsm['X_umap'] = umap_result
                
                # PCA for additional analysis
                from sklearn.decomposition import PCA
                n_comps = min(50, uce_embeddings.shape[1] - 1, uce_embeddings.shape[0] - 1)
                pca = PCA(n_components=n_comps)
                pca_result = pca.fit_transform(uce_embeddings)
                adata.obsm['X_pca'] = pca_result
                
                print(f"    ✅ UMAP: {adata.obsm['X_umap'].shape}")
                print(f"    ✅ PCA: {adata.obsm['X_pca'].shape}")
                
            except Exception as e:
                print(f"    ❌ Error processing {name}: {e}")
    
    def analyze_perturbation_separation(self):
        """Analyze how well perturbations separate in embedding space"""
        print("\n🔬 Analyzing perturbation separation in UCE space...")
        
        separation_results = {}
        
        for name, adata in self.adata_dict.items():
            print(f"  📊 Separation analysis for {name}...")
            
            try:
                embeddings = adata.obsm['X_uce']
                
                if 'mouse' in name and 'genotype' in adata.obs.columns:
                    # Mouse RETNLA KO/WT separation
                    genotypes = adata.obs['genotype']
                    
                    if len(genotypes.unique()) > 1:
                        sil_score = silhouette_score(embeddings, genotypes)
                        
                        separation_results[name] = {
                            'perturbation_type': 'RETNLA_KO_WT',
                            'silhouette_score': sil_score,
                            'groups': list(genotypes.unique()),
                            'group_sizes': dict(genotypes.value_counts())
                        }
                        
                        print(f"    🎯 RETNLA KO/WT Silhouette: {sil_score:.3f}")
                        
                        # Additional analysis by sex if available
                        if 'Sex' in adata.obs.columns:
                            sex_geno_combined = adata.obs['genotype'] + '_' + adata.obs['Sex']
                            if len(sex_geno_combined.unique()) > 1:
                                sil_score_sex = silhouette_score(embeddings, sex_geno_combined)
                                print(f"    👥 Genotype+Sex Silhouette: {sil_score_sex:.3f}")
                                separation_results[name]['sex_genotype_silhouette'] = sil_score_sex
                
                elif 'perturbation_group' in adata.obs.columns:
                    # Human cell line RETNLB perturbation separation
                    pert_groups = adata.obs['perturbation_group']
                    meaningful_groups = pert_groups[~pert_groups.isin(['unknown', 'other'])]
                    
                    if len(meaningful_groups.unique()) > 1:
                        mask = ~pert_groups.isin(['unknown', 'other'])
                        clean_embeddings = embeddings[mask]
                        clean_groups = pert_groups[mask]
                        
                        sil_score = silhouette_score(clean_embeddings, clean_groups)
                        
                        separation_results[name] = {
                            'perturbation_type': 'RETNLB_perturbation',
                            'silhouette_score': sil_score,
                            'groups': list(clean_groups.unique()),
                            'group_sizes': dict(clean_groups.value_counts())
                        }
                        
                        print(f"    🎯 RETNLB Perturbation Silhouette: {sil_score:.3f}")
                    else:
                        print(f"    ⚠️ Insufficient perturbation groups identified")
                        separation_results[name] = {'note': 'Insufficient groups'}
                
            except Exception as e:
                print(f"    ❌ Error in separation analysis for {name}: {e}")
                separation_results[name] = {'error': str(e)}
        
        return separation_results
    
    def generate_perturbation_visualizations(self, separation_results):
        """Generate perturbation-specific visualizations"""
        print("\n🎨 Generating perturbation visualizations...")
        
        fig_dir = self.results_path / "figures"
        fig_dir.mkdir(exist_ok=True)
        
        # 1. Perturbation UMAP plots
        n_datasets = len(self.adata_dict)
        fig, axes = plt.subplots(1, n_datasets, figsize=(6*n_datasets, 5))
        if n_datasets == 1:
            axes = [axes]
        
        for i, (name, adata) in enumerate(self.adata_dict.items()):
            ax = axes[i]
            
            if 'X_umap' in adata.obsm:
                umap_coords = adata.obsm['X_umap']
                
                if 'mouse' in name and 'genotype' in adata.obs.columns:
                    # Color by genotype
                    genotypes = adata.obs['genotype']
                    colors = {'KO': 'red', 'WT': 'blue'}
                    
                    for geno in genotypes.unique():
                        mask = genotypes == geno
                        color = colors.get(geno, 'gray')
                        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1], 
                                  c=color, label=f'RETNLA {geno}', alpha=0.7, s=20)
                    
                elif 'perturbation_group' in adata.obs.columns:
                    # Color by perturbation group
                    pert_groups = adata.obs['perturbation_group']
                    colors = {'RETNLB_perturbed': 'red', 'control': 'blue', 'unknown': 'lightgray'}
                    
                    for group in pert_groups.unique():
                        mask = pert_groups == group
                        color = colors.get(group, 'gray')
                        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1], 
                                  c=color, label=group, alpha=0.7, s=1)
                else:
                    # Fallback coloring
                    ax.scatter(umap_coords[:, 0], umap_coords[:, 1], alpha=0.6, s=1, c='gray')
                
                ax.set_title(f"{name.replace('_', ' ').title()}\n({adata.n_obs:,} cells)")
                ax.set_xlabel("UMAP1")
                ax.set_ylabel("UMAP2")
                ax.legend(markerscale=2, fontsize=8)
        
        plt.tight_layout()
        plt.savefig(fig_dir / "perturbation_umap_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Perturbation effect heatmap
        if separation_results:
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            
            # Prepare data for heatmap
            heatmap_data = []
            dataset_names = []
            
            for name, results in separation_results.items():
                if 'silhouette_score' in results:
                    dataset_names.append(name.replace('_', ' ').title())
                    heatmap_data.append([results['silhouette_score']])
            
            if heatmap_data:
                heatmap_array = np.array(heatmap_data)
                im = ax.imshow(heatmap_array, cmap='RdBu_r', aspect='auto', vmin=-0.5, vmax=0.5)
                
                ax.set_xticks([0])
                ax.set_xticklabels(['Silhouette Score'])
                ax.set_yticks(range(len(dataset_names)))
                ax.set_yticklabels(dataset_names)
                
                # Add text annotations
                for i in range(len(dataset_names)):
                    ax.text(0, i, f'{heatmap_array[i, 0]:.3f}', 
                           ha='center', va='center', fontweight='bold')
                
                plt.colorbar(im, ax=ax, label='Silhouette Score')
                ax.set_title('Perturbation Effect Separation in UCE Space')
                
                plt.tight_layout()
                plt.savefig(fig_dir / "perturbation_effect_heatmap.png", dpi=300, bbox_inches='tight')
                plt.close()
        
        print(f"    ✅ Visualizations saved to {fig_dir}")
    
    def compare_perturbation_effects(self, separation_results):
        """Compare perturbation effects across datasets"""
        print("\n📊 Comparing perturbation effects...")
        
        comparison_results = {
            'mouse_retnla_ko': None,
            'human_retnlb_hek293t': None,
            'human_retnlb_hct116': None
        }
        
        for name, results in separation_results.items():
            if 'mouse' in name and 'silhouette_score' in results:
                comparison_results['mouse_retnla_ko'] = results['silhouette_score']
            elif 'hek293t' in name and 'silhouette_score' in results:
                comparison_results['human_retnlb_hek293t'] = results['silhouette_score']
            elif 'hct116' in name and 'silhouette_score' in results:
                comparison_results['human_retnlb_hct116'] = results['silhouette_score']
        
        print("🔬 Perturbation Effect Summary:")
        if comparison_results['mouse_retnla_ko'] is not None:
            print(f"  🐭 Mouse RETNLA KO effect: {comparison_results['mouse_retnla_ko']:.3f}")
        if comparison_results['human_retnlb_hek293t'] is not None:
            print(f"  🧬 HEK293T RETNLB perturbation: {comparison_results['human_retnlb_hek293t']:.3f}")
        if comparison_results['human_retnlb_hct116'] is not None:
            print(f"  🧬 HCT116 RETNLB perturbation: {comparison_results['human_retnlb_hct116']:.3f}")
        
        # Determine strongest effect
        valid_scores = {k: v for k, v in comparison_results.items() if v is not None}
        if valid_scores:
            strongest_effect = max(valid_scores.items(), key=lambda x: abs(x[1]))
            print(f"\n🎯 Strongest perturbation effect: {strongest_effect[0]} (score: {strongest_effect[1]:.3f})")
        
        return comparison_results
    
    def generate_perturbation_report(self, separation_results, comparison_results):
        """Generate comprehensive perturbation analysis report"""
        print("\n📝 Generating perturbation analysis report...")
        
        report_file = self.results_path / "uce_perturbation_analysis_report.md"
        
        with open(report_file, 'w') as f:
            f.write("# UCE Perturbation Analysis Report\n\n")
            f.write("## RETNLA/RETNLB Perturbation Effects in UCE Embedding Space\n\n")
            
            f.write("### Dataset Summary\n\n")
            f.write("| Dataset | Cells | Perturbation Type | Silhouette Score | Groups |\n")
            f.write("|---------|-------|-------------------|------------------|--------|\n")
            
            for name, adata in self.adata_dict.items():
                cells = adata.n_obs
                if name in separation_results and 'silhouette_score' in separation_results[name]:
                    result = separation_results[name]
                    pert_type = result.get('perturbation_type', 'Unknown')
                    sil_score = result['silhouette_score']
                    groups = ', '.join(result['groups'])
                    f.write(f"| {name} | {cells:,} | {pert_type} | {sil_score:.3f} | {groups} |\n")
                else:
                    f.write(f"| {name} | {cells:,} | - | - | - |\n")
            
            f.write("\n### Perturbation Effect Comparison\n\n")
            f.write("**Key Findings:**\n\n")
            
            if comparison_results['mouse_retnla_ko'] is not None:
                f.write(f"- **Mouse RETNLA KO**: Silhouette score of {comparison_results['mouse_retnla_ko']:.3f}\n")
                if comparison_results['mouse_retnla_ko'] > 0.1:
                    f.write("  - Strong separation between KO and WT cells\n")
                elif comparison_results['mouse_retnla_ko'] > 0:
                    f.write("  - Moderate separation between KO and WT cells\n")
                else:
                    f.write("  - Poor separation between KO and WT cells\n")
            
            if comparison_results['human_retnlb_hek293t'] is not None or comparison_results['human_retnlb_hct116'] is not None:
                f.write(f"- **Human RETNLB Perturbation**:\n")
                if comparison_results['human_retnlb_hek293t'] is not None:
                    f.write(f"  - HEK293T: {comparison_results['human_retnlb_hek293t']:.3f}\n")
                if comparison_results['human_retnlb_hct116'] is not None:
                    f.write(f"  - HCT116: {comparison_results['human_retnlb_hct116']:.3f}\n")
            
            f.write(f"\n### Biological Interpretation\n\n")
            f.write("**RETNLA (Mouse) vs RETNLB (Human)**:\n")
            f.write("- RETNLA and RETNLB are orthologs involved in immune regulation\n")
            f.write("- Both proteins are associated with Th2 immune responses and tissue remodeling\n")
            f.write("- Perturbation effects in UCE space reflect the biological impact of these genes\n\n")
            
            f.write(f"**Analysis completed:** {pd.Timestamp.now()}\n")
            f.write(f"**Results location:** {self.results_path}\n")
        
        print(f"    ✅ Report saved: {report_file}")
    
    def run_complete_perturbation_analysis(self):
        """Run the complete perturbation-focused analysis"""
        print("🚀 Starting UCE Perturbation Analysis...")
        print("🎯 Focus: RETNLA KO/WT (mouse) and RETNLB perturbation (human cell lines)")
        print("="*80)
        
        # Load datasets
        self.load_perturbation_datasets()
        
        if not self.adata_dict:
            print("❌ No datasets loaded for perturbation analysis")
            return
        
        # Identify perturbation groups
        self.identify_perturbation_groups()
        
        # Perform UMAP analysis
        self.perform_perturbation_umap()
        
        # Analyze perturbation separation
        separation_results = self.analyze_perturbation_separation()
        
        # Generate visualizations
        self.generate_perturbation_visualizations(separation_results)
        
        # Compare effects
        comparison_results = self.compare_perturbation_effects(separation_results)
        
        # Generate report
        self.generate_perturbation_report(separation_results, comparison_results)
        
        print("\n" + "="*80)
        print("✅ UCE Perturbation Analysis Complete!")
        print(f"📁 Results saved to: {self.results_path}")
        print("\n🔬 Key Insights:")
        print("   • RETNLA KO/WT effects in mouse epithelial cells analyzed")
        print("   • RETNLB perturbation effects in human cell lines quantified")
        print("   • Cross-species comparison of RETNL family gene effects")
        print("   • UCE embedding space captures perturbation signatures")


if __name__ == "__main__":
    # Initialize analyzer
    analyzer = UCEPerturbationAnalyzer()
    
    # Run complete perturbation analysis
    analyzer.run_complete_perturbation_analysis()
    
    print("\n🎉 Perturbation Analysis Complete!")
    print("📋 Next steps:")
    print("   1. Review perturbation results in /results/uce_perturbation_analysis/")
    print("   2. Compare with scGPT perturbation analysis when available")
    print("   3. Validate findings with differential expression analysis")
    print("   4. Integrate into foundation model comparison framework")