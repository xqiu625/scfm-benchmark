#!/usr/bin/env python3
"""
Comprehensive validation of Geneformer embeddings
Checks quality, validity, and biological relevance
"""

import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from scipy.stats import pearsonr
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_geneformer_embeddings(base_path="/bigdata/godziklab/shared/Xinru/302006/results/geneformer/embeddings"):
    """Load and inspect Geneformer embedding files"""
    
    logger.info("Loading Geneformer embeddings...")
    
    # Load embeddings
    embeddings_file = f"{base_path}/retnla_ortholog_geneformer_embeddings.npy"
    embeddings = np.load(embeddings_file)
    logger.info(f"Embeddings shape: {embeddings.shape}")
    
    # Load metadata
    metadata_file = f"{base_path}/cell_metadata.csv"
    metadata = pd.read_csv(metadata_file)
    logger.info(f"Metadata shape: {metadata.shape}")
    
    # Load processing info
    processing_file = f"{base_path}/processing_metadata.pkl"
    with open(processing_file, 'rb') as f:
        processing_info = pickle.load(f)
    
    return embeddings, metadata, processing_info

def validate_embedding_quality(embeddings):
    """Comprehensive embedding quality assessment"""
    
    logger.info("\n=== EMBEDDING QUALITY VALIDATION ===")
    
    # Basic statistics
    logger.info(f"Embedding shape: {embeddings.shape}")
    logger.info(f"Data type: {embeddings.dtype}")
    logger.info(f"Memory usage: {embeddings.nbytes / 1024**2:.2f} MB")
    
    # Check for problematic values
    nan_count = np.isnan(embeddings).sum()
    inf_count = np.isinf(embeddings).sum()
    zero_rows = (embeddings == 0).all(axis=1).sum()
    
    logger.info(f"NaN values: {nan_count}")
    logger.info(f"Infinite values: {inf_count}")
    logger.info(f"All-zero rows: {zero_rows}")
    
    if zero_rows > 0:
        logger.warning(f"Found {zero_rows} cells with all-zero embeddings - these likely failed during processing")
    
    # Statistical properties
    mean_vals = np.mean(embeddings, axis=0)
    std_vals = np.std(embeddings, axis=0)
    
    logger.info(f"Mean embedding values: min={mean_vals.min():.4f}, max={mean_vals.max():.4f}, mean={mean_vals.mean():.4f}")
    logger.info(f"Std embedding values: min={std_vals.min():.4f}, max={std_vals.max():.4f}, mean={std_vals.mean():.4f}")
    
    # Check for constant features
    constant_features = (std_vals < 1e-8).sum()
    logger.info(f"Near-constant features: {constant_features}")
    
    # Distribution analysis
    logger.info(f"Overall embedding range: [{embeddings.min():.4f}, {embeddings.max():.4f}]")
    logger.info(f"Overall embedding mean: {embeddings.mean():.4f}")
    logger.info(f"Overall embedding std: {embeddings.std():.4f}")
    
    return {
        'valid_embeddings': embeddings.shape[0] - zero_rows,
        'failed_embeddings': zero_rows,
        'has_nan': nan_count > 0,
        'has_inf': inf_count > 0,
        'constant_features': constant_features,
        'mean_range': (mean_vals.min(), mean_vals.max()),
        'std_range': (std_vals.min(), std_vals.max())
    }

def analyze_biological_structure(embeddings, metadata):
    """Analyze whether embeddings capture biological structure"""
    
    logger.info("\n=== BIOLOGICAL STRUCTURE ANALYSIS ===")
    
    # Remove failed embeddings (all zeros) for analysis
    valid_mask = ~(embeddings == 0).all(axis=1)
    valid_embeddings = embeddings[valid_mask]
    valid_metadata = metadata[valid_mask].reset_index(drop=True)
    
    logger.info(f"Analyzing {valid_embeddings.shape[0]} valid embeddings")
    
    if valid_embeddings.shape[0] == 0:
        logger.error("No valid embeddings found for analysis")
        return {}
    
    results = {}
    
    # 1. Cross-species separation
    if 'species' in valid_metadata.columns:
        species_labels = valid_metadata['species']
        unique_species = species_labels.unique()
        logger.info(f"Species found: {unique_species}")
        
        if len(unique_species) > 1:
            # Calculate silhouette score for species separation
            species_encoded = pd.Categorical(species_labels).codes
            species_silhouette = silhouette_score(valid_embeddings, species_encoded)
            results['species_separation'] = species_silhouette
            logger.info(f"Species separation silhouette score: {species_silhouette:.4f}")
        
        # Count per species
        species_counts = species_labels.value_counts()
        logger.info(f"Cell counts per species:")
        for species, count in species_counts.items():
            logger.info(f"  {species}: {count} cells")
    
    # 2. Treatment status separation
    if 'treatment_status' in valid_metadata.columns:
        treatment_labels = valid_metadata['treatment_status']
        unique_treatments = treatment_labels.unique()
        logger.info(f"Treatment statuses: {unique_treatments}")
        
        if len(unique_treatments) > 1:
            treatment_encoded = pd.Categorical(treatment_labels).codes
            treatment_silhouette = silhouette_score(valid_embeddings, treatment_encoded)
            results['treatment_separation'] = treatment_silhouette
            logger.info(f"Treatment separation silhouette score: {treatment_silhouette:.4f}")
        
        # Count per treatment
        treatment_counts = treatment_labels.value_counts()
        logger.info(f"Cell counts per treatment:")
        for treatment, count in treatment_counts.items():
            logger.info(f"  {treatment}: {count} cells")
    
    # 3. Gene target separation
    if 'gene_target' in valid_metadata.columns:
        gene_labels = valid_metadata['gene_target']
        unique_genes = gene_labels.unique()
        logger.info(f"Gene targets: {unique_genes}")
        
        if len(unique_genes) > 1:
            gene_encoded = pd.Categorical(gene_labels).codes
            gene_silhouette = silhouette_score(valid_embeddings, gene_encoded)
            results['gene_separation'] = gene_silhouette
            logger.info(f"Gene target separation silhouette score: {gene_silhouette:.4f}")
        
        # Count per gene
        gene_counts = gene_labels.value_counts()
        logger.info(f"Cell counts per gene target:")
        for gene, count in gene_counts.items():
            logger.info(f"  {gene}: {count} cells")
    
    # 4. Dataset separation
    if 'dataset' in valid_metadata.columns:
        dataset_labels = valid_metadata['dataset']
        unique_datasets = dataset_labels.unique()
        logger.info(f"Datasets: {unique_datasets}")
        
        if len(unique_datasets) > 1:
            dataset_encoded = pd.Categorical(dataset_labels).codes
            dataset_silhouette = silhouette_score(valid_embeddings, dataset_encoded)
            results['dataset_separation'] = dataset_silhouette
            logger.info(f"Dataset separation silhouette score: {dataset_silhouette:.4f}")
    
    return results

def cross_species_ortholog_analysis(embeddings, metadata):
    """Specific analysis of RETNLA ortholog relationships"""
    
    logger.info("\n=== CROSS-SPECIES ORTHOLOG ANALYSIS ===")
    
    # Remove failed embeddings
    valid_mask = ~(embeddings == 0).all(axis=1)
    valid_embeddings = embeddings[valid_mask]
    valid_metadata = metadata[valid_mask].reset_index(drop=True)
    
    results = {}
    
    # Analyze mouse RETNLA vs human RETNLB/RETN
    mouse_mask = valid_metadata['species'] == 'mouse'
    human_mask = valid_metadata['species'] == 'human'
    
    if mouse_mask.sum() > 0 and human_mask.sum() > 0:
        mouse_embeddings = valid_embeddings[mouse_mask]
        human_embeddings = valid_embeddings[human_mask]
        
        logger.info(f"Mouse cells: {mouse_embeddings.shape[0]}")
        logger.info(f"Human cells: {human_embeddings.shape[0]}")
        
        # Calculate average embeddings per group
        mouse_avg = np.mean(mouse_embeddings, axis=0)
        human_avg = np.mean(human_embeddings, axis=0)
        
        # Cross-species correlation
        cross_species_corr, p_value = pearsonr(mouse_avg, human_avg)
        results['cross_species_correlation'] = cross_species_corr
        results['correlation_pvalue'] = p_value
        
        logger.info(f"Cross-species average embedding correlation: r={cross_species_corr:.4f}, p={p_value:.2e}")
        
        # Analyze perturbation effects within each species
        # Mouse: KO vs WT
        mouse_metadata = valid_metadata[mouse_mask]
        mouse_ko = mouse_metadata['treatment_status'].str.contains('KO', na=False)
        mouse_wt = mouse_metadata['treatment_status'].str.contains('WT', na=False)
        
        if mouse_ko.sum() > 0 and mouse_wt.sum() > 0:
            mouse_ko_avg = np.mean(mouse_embeddings[mouse_ko], axis=0)
            mouse_wt_avg = np.mean(mouse_embeddings[mouse_wt], axis=0)
            mouse_perturbation_corr, _ = pearsonr(mouse_ko_avg, mouse_wt_avg)
            results['mouse_perturbation_correlation'] = mouse_perturbation_corr
            
            logger.info(f"Mouse KO vs WT correlation: r={mouse_perturbation_corr:.4f}")
            logger.info(f"  KO cells: {mouse_ko.sum()}, WT cells: {mouse_wt.sum()}")
        
        # Human: targeted vs expressing
        human_metadata = valid_metadata[human_mask]
        human_targeted = human_metadata['treatment_status'].str.contains('targeted', na=False)
        human_expressing = human_metadata['treatment_status'].str.contains('expressing', na=False)
        
        if human_targeted.sum() > 0 and human_expressing.sum() > 0:
            human_targeted_avg = np.mean(human_embeddings[human_targeted], axis=0)
            human_expressing_avg = np.mean(human_embeddings[human_expressing], axis=0)
            human_perturbation_corr, _ = pearsonr(human_targeted_avg, human_expressing_avg)
            results['human_perturbation_correlation'] = human_perturbation_corr
            
            logger.info(f"Human targeted vs expressing correlation: r={human_perturbation_corr:.4f}")
            logger.info(f"  Targeted cells: {human_targeted.sum()}, Expressing cells: {human_expressing.sum()}")
    
    return results

def generate_visualization(embeddings, metadata, output_path="/bigdata/godziklab/shared/Xinru/302006/results/geneformer/embeddings"):
    """Generate visualization of embeddings"""
    
    logger.info("\n=== GENERATING VISUALIZATIONS ===")
    
    # Remove failed embeddings for visualization
    valid_mask = ~(embeddings == 0).all(axis=1)
    valid_embeddings = embeddings[valid_mask]
    valid_metadata = metadata[valid_mask].reset_index(drop=True)
    
    if valid_embeddings.shape[0] < 10:
        logger.warning("Too few valid embeddings for meaningful visualization")
        return
    
    # PCA
    logger.info("Computing PCA...")
    pca = PCA(n_components=2)
    pca_embeddings = pca.fit_transform(valid_embeddings)
    
    # Create visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # PCA by species
    if 'species' in valid_metadata.columns:
        species_colors = {'mouse': 'red', 'human': 'blue'}
        for species in valid_metadata['species'].unique():
            mask = valid_metadata['species'] == species
            ax1.scatter(pca_embeddings[mask, 0], pca_embeddings[mask, 1], 
                       c=species_colors.get(species, 'gray'), label=species, alpha=0.6)
        ax1.set_title('PCA by Species')
        ax1.legend()
    
    # PCA by treatment
    if 'treatment_status' in valid_metadata.columns:
        treatment_colors = plt.cm.Set1(np.linspace(0, 1, len(valid_metadata['treatment_status'].unique())))
        for i, treatment in enumerate(valid_metadata['treatment_status'].unique()):
            mask = valid_metadata['treatment_status'] == treatment
            ax2.scatter(pca_embeddings[mask, 0], pca_embeddings[mask, 1], 
                       c=[treatment_colors[i]], label=treatment, alpha=0.6)
        ax2.set_title('PCA by Treatment Status')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # PCA by gene target
    if 'gene_target' in valid_metadata.columns:
        gene_colors = plt.cm.Set2(np.linspace(0, 1, len(valid_metadata['gene_target'].unique())))
        for i, gene in enumerate(valid_metadata['gene_target'].unique()):
            mask = valid_metadata['gene_target'] == gene
            ax3.scatter(pca_embeddings[mask, 0], pca_embeddings[mask, 1], 
                       c=[gene_colors[i]], label=gene, alpha=0.6)
        ax3.set_title('PCA by Gene Target')
        ax3.legend()
    
    # PCA by dataset
    if 'dataset' in valid_metadata.columns:
        dataset_colors = plt.cm.Set3(np.linspace(0, 1, len(valid_metadata['dataset'].unique())))
        for i, dataset in enumerate(valid_metadata['dataset'].unique()):
            mask = valid_metadata['dataset'] == dataset
            ax4.scatter(pca_embeddings[mask, 0], pca_embeddings[mask, 1], 
                       c=[dataset_colors[i]], label=dataset, alpha=0.6)
        ax4.set_title('PCA by Dataset')
        ax4.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/geneformer_embedding_visualization.png", dpi=300, bbox_inches='tight')
    logger.info(f"Visualization saved to {output_path}/geneformer_embedding_visualization.png")

def main():
    """Main validation function"""
    
    logger.info("Starting Geneformer embedding validation...")
    
    try:
        # Load embeddings
        embeddings, metadata, processing_info = load_geneformer_embeddings()
        
        logger.info(f"\nProcessing info:")
        for key, value in processing_info.items():
            logger.info(f"  {key}: {value}")
        
        # Validate quality
        quality_results = validate_embedding_quality(embeddings)
        
        # Analyze biological structure
        bio_results = analyze_biological_structure(embeddings, metadata)
        
        # Cross-species analysis
        ortholog_results = cross_species_ortholog_analysis(embeddings, metadata)
        
        # Generate visualizations
        generate_visualization(embeddings, metadata)
        
        # Summary report
        logger.info("\n=== VALIDATION SUMMARY ===")
        logger.info(f"Total cells processed: {embeddings.shape[0]}")
        logger.info(f"Valid embeddings: {quality_results['valid_embeddings']}")
        logger.info(f"Failed embeddings: {quality_results['failed_embeddings']}")
        logger.info(f"Success rate: {quality_results['valid_embeddings']/embeddings.shape[0]*100:.1f}%")
        
        if bio_results:
            logger.info(f"\nBiological structure scores:")
            for key, value in bio_results.items():
                logger.info(f"  {key}: {value:.4f}")
        
        if ortholog_results:
            logger.info(f"\nCross-species ortholog analysis:")
            for key, value in ortholog_results.items():
                logger.info(f"  {key}: {value:.4f}")
        
        # Quality assessment
        if quality_results['failed_embeddings'] == 0:
            logger.info("\n✅ RESULT: All embeddings generated successfully")
        elif quality_results['failed_embeddings'] < embeddings.shape[0] * 0.1:
            logger.info(f"\n⚠️  RESULT: Mostly successful with {quality_results['failed_embeddings']} failed embeddings")
        else:
            logger.warning(f"\n❌ RESULT: High failure rate - {quality_results['failed_embeddings']} embeddings failed")
        
        if bio_results.get('species_separation', 0) > 0.1:
            logger.info("✅ Embeddings show good cross-species structure")
        elif bio_results.get('species_separation', 0) > 0:
            logger.info("⚠️  Embeddings show weak cross-species structure")
        
        logger.info("\nGeneformer embedding validation complete!")
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        raise

if __name__ == "__main__":
    main()