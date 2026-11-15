#!/usr/bin/env python3
"""
Perturbation-Severity Analysis

Map perturbation embeddings (HEK293T, HCT116) to COVID severity states.
Identify perturbations that shift transcriptional states toward recovery.

Goal: Find drug candidates for cardiovascular treatment post-infection.

Usage:
    python scripts/analysis/perturbation_severity_analysis.py --model UCE

Author: Project 302006
Date: December 2025
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import seaborn as sns
import umap
from datetime import datetime
import json
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300

# Colors
SEVERITY_COLORS = {
    'healthy': '#2ecc71',
    'mild': '#3498db',
    'moderate': '#f39c12',
    'severe': '#e74c3c',
    'fatal': '#8e44ad',
    'recovered': '#1abc9c'
}


class PerturbationSeverityAnalysis:
    """Analyze perturbation effects relative to COVID severity states"""

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006", model='UCE'):
        self.base_dir = Path(base_dir)
        self.embedding_dir = self.base_dir / "02_EMBEDDINGS"
        self.output_dir = self.base_dir / "05_RESULTS/perturbation_analysis"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.model = model

        # Model configurations
        # UCE 33-layer embeddings have HEK293T and HCT116 available
        # scGPT has combined file with dataset column
        # Geneformer has separate files (check availability)
        # TranscriptFormer only has COVID (HEK293T/HCT116 blocked due to gene ID issues)
        self.model_configs = {
            'STATE': {
                'covid_file': self.embedding_dir / "state/covid_state.h5ad",
                'hek293t_file': self.embedding_dir / "state/hek293t_state.h5ad",
                'hct116_file': self.embedding_dir / "state/hct116_state.h5ad",
                'obsm_keys': ['X_state', 'X_STATE', 'X_emb'],
            },
            'UCE': {
                'covid_file': self.embedding_dir / "uce_33layer/covid_uce_33layer.h5ad",
                'hek293t_file': self.embedding_dir / "uce_33layer/hek293t_uce_33layer.h5ad",
                'hct116_file': self.embedding_dir / "uce_33layer/hct116_uce_33layer.h5ad",
                'obsm_keys': ['X_uce', 'X_uce_33layer', 'X_uce_4layer'],
                # Fallback to 4-layer if 33-layer not found
                'fallback_covid': self.embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad",
                'fallback_hek293t': self.embedding_dir / "uce/human_platelet_hek293t_uce_adata.h5ad",
                'fallback_hct116': self.embedding_dir / "uce/human_platelet_hct116_uce_adata.h5ad",
            },
            'scGPT': {
                'covid_file': self.embedding_dir / "scgpt/platelet_scgpt_combined.h5ad",
                'hek293t_file': self.embedding_dir / "scgpt/platelet_scgpt_combined.h5ad",
                'hct116_file': self.embedding_dir / "scgpt/platelet_scgpt_combined.h5ad",
                'obsm_keys': ['X_scGPT', 'X_scgpt'],
            },
            'Geneformer': {
                'covid_file': self.embedding_dir / "geneformer/covid_geneformer.h5ad",
                'hek293t_file': self.embedding_dir / "geneformer/hek293t_geneformer.h5ad",
                'hct116_file': self.embedding_dir / "geneformer/hct116_geneformer.h5ad",
                'obsm_keys': ['X_geneformer', 'X_Geneformer'],
            },
            'TranscriptFormer': {
                'covid_file': self.embedding_dir / "transcriptformer/covid_transcriptformer.h5ad",
                'hek293t_file': None,  # Not available due to gene ID issues
                'hct116_file': None,   # Not available due to gene ID issues
                'obsm_keys': ['X_transcriptformer', 'X_TranscriptFormer'],
            }
        }

        self.data = {}

    def load_data(self):
        """Load COVID and perturbation embeddings"""
        config = self.model_configs[self.model]

        # Load COVID data
        logger.info(f"Loading {self.model} COVID embeddings...")
        covid_file = config['covid_file']

        # Try fallback for UCE if primary not found
        if not covid_file.exists() and 'fallback_covid' in config:
            logger.info(f"  Primary not found, trying fallback...")
            covid_file = config['fallback_covid']

        if covid_file.exists():
            adata = sc.read_h5ad(covid_file)

            # Filter for scGPT combined file
            if self.model == 'scGPT' and 'dataset' in adata.obs.columns:
                adata = adata[adata.obs['dataset'] == 'COVID'].copy()

            # Get embeddings
            obsm_key = self._get_obsm_key(adata, config['obsm_keys'])
            if obsm_key:
                embeddings = np.array(adata.obsm[obsm_key])

                # Get severity labels
                labels = self._get_severity_labels(adata)
                valid_mask = labels != 'unknown'

                self.data['covid'] = {
                    'embeddings': embeddings[valid_mask],
                    'labels': labels[valid_mask],
                    'adata': adata[valid_mask].copy()
                }
                logger.info(f"  COVID: {self.data['covid']['embeddings'].shape}")

        # Load HEK293T data
        logger.info(f"Loading {self.model} HEK293T embeddings...")
        hek_file = config['hek293t_file']

        # Skip if None (e.g., TranscriptFormer)
        if hek_file is None:
            logger.warning(f"  HEK293T embeddings not available for {self.model}")
        else:
            # Try fallback for UCE if primary not found
            if not hek_file.exists() and 'fallback_hek293t' in config:
                logger.info(f"  Primary not found, trying fallback...")
                hek_file = config['fallback_hek293t']

        if hek_file is not None and hek_file.exists():
            adata = sc.read_h5ad(hek_file)

            if self.model == 'scGPT' and 'dataset' in adata.obs.columns:
                adata = adata[adata.obs['dataset'] == 'HEK293T'].copy()

            obsm_key = self._get_obsm_key(adata, config['obsm_keys'])
            if obsm_key:
                embeddings = np.array(adata.obsm[obsm_key])
                perturbations = self._get_perturbation_labels(adata)

                self.data['hek293t'] = {
                    'embeddings': embeddings,
                    'perturbations': perturbations,
                    'adata': adata.copy()
                }
                logger.info(f"  HEK293T: {embeddings.shape}")
                logger.info(f"  Unique perturbations: {len(np.unique(perturbations))}")

        # Load HCT116 data
        logger.info(f"Loading {self.model} HCT116 embeddings...")
        hct_file = config['hct116_file']

        # Skip if None (e.g., TranscriptFormer)
        if hct_file is None:
            logger.warning(f"  HCT116 embeddings not available for {self.model}")
        else:
            # Try fallback for UCE if primary not found
            if not hct_file.exists() and 'fallback_hct116' in config:
                logger.info(f"  Primary not found, trying fallback...")
                hct_file = config['fallback_hct116']

        if hct_file is not None and hct_file.exists():
            adata = sc.read_h5ad(hct_file)

            if self.model == 'scGPT' and 'dataset' in adata.obs.columns:
                adata = adata[adata.obs['dataset'] == 'HCT116'].copy()

            obsm_key = self._get_obsm_key(adata, config['obsm_keys'])
            if obsm_key:
                embeddings = np.array(adata.obsm[obsm_key])
                perturbations = self._get_perturbation_labels(adata)

                self.data['hct116'] = {
                    'embeddings': embeddings,
                    'perturbations': perturbations,
                    'adata': adata.copy()
                }
                logger.info(f"  HCT116: {embeddings.shape}")
                logger.info(f"  Unique perturbations: {len(np.unique(perturbations))}")

        return len(self.data) > 0

    def _get_obsm_key(self, adata, obsm_keys):
        """Find the correct obsm key"""
        for key in obsm_keys:
            if key in adata.obsm:
                return key
        for k in adata.obsm.keys():
            if self.model.lower() in k.lower():
                return k
        return None

    def _get_severity_labels(self, adata):
        """Extract severity labels"""
        for col in ['covid_severity', 'severity', 'Category']:
            if col in adata.obs.columns:
                labels = np.array([str(x) if pd.notna(x) else 'unknown'
                                   for x in adata.obs[col]])
                return labels
        return np.array(['unknown'] * len(adata))

    def _get_perturbation_labels(self, adata):
        """Extract perturbation labels"""
        for col in ['guide_target', 'gene_target', 'perturbation', 'treatment']:
            if col in adata.obs.columns:
                labels = np.array([str(x) if pd.notna(x) else 'control'
                                   for x in adata.obs[col]])
                return labels
        return np.array(['unknown'] * len(adata))

    def compute_severity_centroids(self):
        """Compute centroid for each severity level"""
        if 'covid' not in self.data:
            logger.error("COVID data not loaded")
            return

        embeddings = self.data['covid']['embeddings']
        labels = self.data['covid']['labels']

        # Scale embeddings
        scaler = StandardScaler()
        embeddings_scaled = scaler.fit_transform(embeddings)
        self.scaler = scaler

        # Compute centroids
        self.severity_centroids = {}
        for severity in ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']:
            mask = labels == severity
            if mask.sum() > 0:
                self.severity_centroids[severity] = embeddings_scaled[mask].mean(axis=0)
                logger.info(f"  {severity}: {mask.sum()} cells")

        # Compute "recovery direction" - vector from severe to recovered
        if 'severe' in self.severity_centroids and 'recovered' in self.severity_centroids:
            self.recovery_direction = (self.severity_centroids['recovered'] -
                                       self.severity_centroids['severe'])
            self.recovery_direction = self.recovery_direction / np.linalg.norm(self.recovery_direction)
            logger.info("Computed recovery direction vector")

        # Compute "disease direction" - healthy to severe
        if 'healthy' in self.severity_centroids and 'severe' in self.severity_centroids:
            self.disease_direction = (self.severity_centroids['severe'] -
                                      self.severity_centroids['healthy'])
            self.disease_direction = self.disease_direction / np.linalg.norm(self.disease_direction)
            logger.info("Computed disease direction vector")

    def compute_perturbation_centroids(self, dataset='hek293t'):
        """Compute centroid for each perturbation"""
        if dataset not in self.data:
            logger.error(f"{dataset} data not loaded")
            return {}

        embeddings = self.data[dataset]['embeddings']
        perturbations = self.data[dataset]['perturbations']

        # Scale using COVID scaler
        embeddings_scaled = self.scaler.transform(embeddings)

        # Compute centroids per perturbation
        perturbation_centroids = {}
        unique_perts = np.unique(perturbations)

        for pert in unique_perts:
            mask = perturbations == pert
            if mask.sum() >= 5:  # Minimum cells
                perturbation_centroids[pert] = {
                    'centroid': embeddings_scaled[mask].mean(axis=0),
                    'n_cells': mask.sum()
                }

        logger.info(f"Computed centroids for {len(perturbation_centroids)} perturbations in {dataset}")
        return perturbation_centroids

    def score_perturbations(self, dataset='hek293t'):
        """Score perturbations by their similarity to severity states"""
        perturbation_centroids = self.compute_perturbation_centroids(dataset)

        results = []

        for pert, data in perturbation_centroids.items():
            centroid = data['centroid']
            n_cells = data['n_cells']

            # Distance to each severity centroid
            severity_distances = {}
            severity_similarities = {}
            for severity, sev_centroid in self.severity_centroids.items():
                dist = np.linalg.norm(centroid - sev_centroid)
                sim = 1 / (1 + dist)  # Convert to similarity
                severity_distances[severity] = dist
                severity_similarities[severity] = sim

            # Recovery score: projection onto recovery direction
            recovery_score = np.dot(centroid - self.severity_centroids['severe'],
                                   self.recovery_direction) if hasattr(self, 'recovery_direction') else 0

            # Disease score: projection onto disease direction
            disease_score = np.dot(centroid - self.severity_centroids['healthy'],
                                  self.disease_direction) if hasattr(self, 'disease_direction') else 0

            # Closest severity
            closest_severity = min(severity_distances, key=severity_distances.get)

            results.append({
                'perturbation': pert,
                'dataset': dataset,
                'n_cells': n_cells,
                'closest_severity': closest_severity,
                'recovery_score': recovery_score,
                'disease_score': disease_score,
                **{f'dist_{k}': v for k, v in severity_distances.items()},
                **{f'sim_{k}': v for k, v in severity_similarities.items()}
            })

        return pd.DataFrame(results)

    def find_therapeutic_candidates(self):
        """Find perturbations that promote recovery-like states"""
        logger.info("\n" + "="*60)
        logger.info("Finding therapeutic candidates...")
        logger.info("="*60)

        all_results = []

        for dataset in ['hek293t', 'hct116']:
            if dataset in self.data:
                df = self.score_perturbations(dataset)
                all_results.append(df)

        if not all_results:
            logger.error("No perturbation data available")
            return None

        results_df = pd.concat(all_results, ignore_index=True)

        # Sort by recovery score (higher = more recovery-like)
        results_df = results_df.sort_values('recovery_score', ascending=False)

        # Top candidates: high recovery score, close to recovered/healthy
        therapeutic_candidates = results_df[
            (results_df['recovery_score'] > 0) &
            (results_df['closest_severity'].isin(['healthy', 'recovered', 'mild']))
        ].head(50)

        # Harmful perturbations: low recovery score, close to severe/fatal
        harmful_perturbations = results_df[
            (results_df['recovery_score'] < 0) &
            (results_df['closest_severity'].isin(['severe', 'fatal']))
        ].head(50)

        logger.info(f"\nTop Therapeutic Candidates (promote recovery):")
        logger.info("-" * 50)
        for _, row in therapeutic_candidates.head(10).iterrows():
            logger.info(f"  {row['perturbation']}: recovery={row['recovery_score']:.3f}, "
                       f"closest={row['closest_severity']}, n={row['n_cells']}")

        logger.info(f"\nPotentially Harmful Perturbations (worsen severity):")
        logger.info("-" * 50)
        for _, row in harmful_perturbations.head(10).iterrows():
            logger.info(f"  {row['perturbation']}: recovery={row['recovery_score']:.3f}, "
                       f"closest={row['closest_severity']}, n={row['n_cells']}")

        return results_df, therapeutic_candidates, harmful_perturbations

    def visualize_perturbation_landscape(self, n_sample=5000):
        """Create UMAP visualization of perturbations in severity space"""
        if 'covid' not in self.data:
            return

        logger.info("Creating perturbation landscape visualization...")

        # Combine COVID and perturbation data
        all_embeddings = [self.data['covid']['embeddings']]
        all_labels = list(self.data['covid']['labels'])
        all_sources = ['COVID'] * len(self.data['covid']['embeddings'])

        for dataset in ['hek293t', 'hct116']:
            if dataset in self.data:
                all_embeddings.append(self.data[dataset]['embeddings'])
                all_labels.extend(list(self.data[dataset]['perturbations']))
                all_sources.extend([dataset.upper()] * len(self.data[dataset]['embeddings']))

        all_embeddings = np.vstack(all_embeddings)
        all_labels = np.array(all_labels)
        all_sources = np.array(all_sources)

        # Subsample if needed
        if len(all_embeddings) > n_sample:
            np.random.seed(42)
            idx = np.random.choice(len(all_embeddings), n_sample, replace=False)
            all_embeddings = all_embeddings[idx]
            all_labels = all_labels[idx]
            all_sources = all_sources[idx]

        # Scale and UMAP
        scaler = StandardScaler()
        embeddings_scaled = scaler.fit_transform(all_embeddings)

        reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=30, min_dist=0.3)
        coords = reducer.fit_transform(embeddings_scaled)

        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        # Plot 1: COVID severity
        ax = axes[0]
        covid_mask = all_sources == 'COVID'

        for severity in ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']:
            mask = covid_mask & (all_labels == severity)
            if mask.sum() > 0:
                ax.scatter(coords[mask, 0], coords[mask, 1],
                          c=SEVERITY_COLORS.get(severity, 'gray'),
                          label=severity, alpha=0.5, s=10)

        # Overlay perturbation data
        pert_mask = ~covid_mask
        ax.scatter(coords[pert_mask, 0], coords[pert_mask, 1],
                  c='black', marker='x', s=20, alpha=0.3, label='Perturbations')

        ax.set_title(f'{self.model}: COVID Severity + Perturbations', fontweight='bold')
        ax.legend(loc='best', markerscale=2)
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')

        # Plot 2: Data source
        ax = axes[1]
        source_colors = {'COVID': '#3498db', 'HEK293T': '#e74c3c', 'HCT116': '#2ecc71'}
        for source in ['COVID', 'HEK293T', 'HCT116']:
            mask = all_sources == source
            if mask.sum() > 0:
                ax.scatter(coords[mask, 0], coords[mask, 1],
                          c=source_colors.get(source, 'gray'),
                          label=f'{source} (n={mask.sum()})', alpha=0.5, s=10)

        ax.set_title(f'{self.model}: Data Sources', fontweight='bold')
        ax.legend(loc='best', markerscale=2)
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')

        plt.tight_layout()
        plt.savefig(self.output_dir / f'{self.model.lower()}_perturbation_landscape.png')
        plt.savefig(self.output_dir / f'{self.model.lower()}_perturbation_landscape.pdf')
        plt.close()
        logger.info(f"Saved: {self.model.lower()}_perturbation_landscape.png")

    def plot_therapeutic_rankings(self, results_df, top_n=30):
        """Plot ranking of therapeutic candidates"""
        if results_df is None or len(results_df) == 0:
            return

        fig, axes = plt.subplots(1, 2, figsize=(16, 10))

        # Top therapeutic (high recovery score)
        ax = axes[0]
        top_therapeutic = results_df.nlargest(top_n, 'recovery_score')
        colors = [SEVERITY_COLORS.get(s, 'gray') for s in top_therapeutic['closest_severity']]

        y_pos = np.arange(len(top_therapeutic))
        ax.barh(y_pos, top_therapeutic['recovery_score'], color=colors, alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_therapeutic['perturbation'], fontsize=8)
        ax.set_xlabel('Recovery Score', fontweight='bold')
        ax.set_title(f'Top {top_n} Therapeutic Candidates\n(Promote Recovery)', fontweight='bold')
        ax.invert_yaxis()

        # Bottom (harmful - low recovery score)
        ax = axes[1]
        bottom_harmful = results_df.nsmallest(top_n, 'recovery_score')
        colors = [SEVERITY_COLORS.get(s, 'gray') for s in bottom_harmful['closest_severity']]

        y_pos = np.arange(len(bottom_harmful))
        ax.barh(y_pos, bottom_harmful['recovery_score'], color=colors, alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(bottom_harmful['perturbation'], fontsize=8)
        ax.set_xlabel('Recovery Score', fontweight='bold')
        ax.set_title(f'Bottom {top_n} Perturbations\n(Worsen Severity)', fontweight='bold')
        ax.invert_yaxis()

        # Add legend
        legend_elements = [plt.Rectangle((0,0),1,1, facecolor=c, label=s)
                         for s, c in SEVERITY_COLORS.items()]
        fig.legend(handles=legend_elements, loc='upper center', ncol=6,
                  bbox_to_anchor=(0.5, 0.02))

        plt.tight_layout(rect=[0, 0.05, 1, 1])
        plt.savefig(self.output_dir / f'{self.model.lower()}_therapeutic_rankings.png')
        plt.savefig(self.output_dir / f'{self.model.lower()}_therapeutic_rankings.pdf')
        plt.close()
        logger.info(f"Saved: {self.model.lower()}_therapeutic_rankings.png")

    def save_results(self, results_df, therapeutic, harmful):
        """Save all results"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Full results
        results_df.to_csv(self.output_dir / f'{self.model.lower()}_perturbation_scores_{timestamp}.csv',
                         index=False)

        # Therapeutic candidates
        if therapeutic is not None:
            therapeutic.to_csv(self.output_dir / f'{self.model.lower()}_therapeutic_candidates_{timestamp}.csv',
                              index=False)

        # Harmful perturbations
        if harmful is not None:
            harmful.to_csv(self.output_dir / f'{self.model.lower()}_harmful_perturbations_{timestamp}.csv',
                          index=False)

        # Summary report
        with open(self.output_dir / f'{self.model.lower()}_analysis_report_{timestamp}.txt', 'w') as f:
            f.write(f"PERTURBATION-SEVERITY ANALYSIS REPORT\n")
            f.write(f"Model: {self.model}\n")
            f.write(f"Date: {timestamp}\n")
            f.write("="*60 + "\n\n")

            f.write("DATA SUMMARY:\n")
            for dataset, data in self.data.items():
                f.write(f"  {dataset}: {data['embeddings'].shape[0]} cells\n")

            f.write(f"\nTOTAL PERTURBATIONS ANALYZED: {len(results_df)}\n")
            f.write(f"THERAPEUTIC CANDIDATES: {len(therapeutic) if therapeutic is not None else 0}\n")
            f.write(f"HARMFUL PERTURBATIONS: {len(harmful) if harmful is not None else 0}\n")

            if therapeutic is not None and len(therapeutic) > 0:
                f.write("\n\nTOP 20 THERAPEUTIC CANDIDATES:\n")
                f.write("-"*60 + "\n")
                for _, row in therapeutic.head(20).iterrows():
                    f.write(f"{row['perturbation']}: recovery={row['recovery_score']:.3f}, "
                           f"closest={row['closest_severity']}\n")

        logger.info(f"Results saved to {self.output_dir}")

    def run(self):
        """Run full analysis pipeline"""
        logger.info("="*60)
        logger.info(f"PERTURBATION-SEVERITY ANALYSIS: {self.model}")
        logger.info("="*60)

        # Load data
        if not self.load_data():
            logger.error("Failed to load data")
            return

        # Compute severity centroids
        self.compute_severity_centroids()

        # Find therapeutic candidates
        results_df, therapeutic, harmful = self.find_therapeutic_candidates()

        if results_df is not None:
            # Visualizations
            self.visualize_perturbation_landscape()
            self.plot_therapeutic_rankings(results_df)

            # Save results
            self.save_results(results_df, therapeutic, harmful)

        logger.info("\n" + "="*60)
        logger.info("ANALYSIS COMPLETE")
        logger.info("="*60)


def main():
    parser = argparse.ArgumentParser(description='Perturbation-Severity Analysis')
    parser.add_argument('--model', type=str, default='UCE',
                        choices=['STATE', 'UCE', 'scGPT', 'Geneformer', 'TranscriptFormer'],
                        help='Foundation model to use')
    parser.add_argument('--base-dir', type=str,
                        default="/bigdata/godziklab/shared/Xinru/302006",
                        help='Base directory')
    args = parser.parse_args()

    analysis = PerturbationSeverityAnalysis(base_dir=args.base_dir, model=args.model)
    analysis.run()


if __name__ == "__main__":
    main()
