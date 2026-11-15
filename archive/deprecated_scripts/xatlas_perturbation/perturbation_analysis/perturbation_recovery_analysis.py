#!/usr/bin/env python3
"""
Perturbation Recovery Analysis

Analyze which gene perturbations shift cells toward a "recovery" direction
based on COVID severity embedding space.

Key Analysis:
1. Define recovery direction from COVID data (severe → recovered)
2. Compute perturbation effects as shifts in embedding space
3. Score each perturbation by its alignment with recovery direction
4. Rank therapeutic candidates

Usage:
    python perturbation_recovery_analysis.py --model STATE
    python perturbation_recovery_analysis.py --model UCE

Author: Project 302006
Date: December 2025
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from datetime import datetime
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import argparse
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "perturbation_recovery"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Severity colors
SEVERITY_COLORS = {
    'healthy': '#2ecc71', 'mild': '#3498db', 'moderate': '#f39c12',
    'severe': '#e74c3c', 'fatal': '#8e44ad', 'recovered': '#1abc9c'
}


class PerturbationRecoveryAnalysis:
    """Analyze perturbation effects in COVID severity space."""

    def __init__(self, model='STATE'):
        self.model = model
        self.data = {}
        self.scaler = None
        self.severity_centroids = {}
        self.recovery_direction = None
        self.disease_direction = None

        # Model-specific paths
        self.paths = self._get_model_paths()

    def _get_model_paths(self):
        """Get file paths for the model."""
        if self.model == 'STATE':
            return {
                'covid': EMB_DIR / "state/covid_state.h5ad",
                # X-Atlas embeddings (primary - full perturbation data)
                'hek293t_xatlas_full': EMB_DIR / "state/hek293t_xatlas_full_state.h5ad",
                'hct116_xatlas_full': EMB_DIR / "state/hct116_xatlas_full_state.h5ad",
                'hek293t_xatlas_500k': EMB_DIR / "state/hek293t_xatlas_500k_state.h5ad",
                'hct116_xatlas_500k': EMB_DIR / "state/hct116_xatlas_500k_state.h5ad",
                # Legacy full embeddings
                'hek293t': EMB_DIR / "state/hek293t_full_state.h5ad",
                'hct116': EMB_DIR / "state/hct116_full_state.h5ad",
                # Fallback to existing files
                'hek293t_alt': EMB_DIR / "state/hek293t_state.h5ad",
                'hct116_alt': EMB_DIR / "state/hct116_state.h5ad",
                'obsm_keys': ['X_state', 'X_STATE', 'X_emb'],
            }
        else:  # UCE
            return {
                'covid': EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
                # X-Atlas embeddings (primary - full perturbation data)
                'hek293t_xatlas_full': EMB_DIR / "uce/hek293t_xatlas_full_uce.h5ad",
                'hct116_xatlas_full': EMB_DIR / "uce/hct116_xatlas_full_uce.h5ad",
                'hek293t_xatlas_500k': EMB_DIR / "uce/hek293t_xatlas_500k_uce.h5ad",
                'hct116_xatlas_500k': EMB_DIR / "uce/hct116_xatlas_500k_uce.h5ad",
                # Legacy full embeddings
                'hek293t': EMB_DIR / "uce/hek293t_full_uce.h5ad",
                'hct116': EMB_DIR / "uce/hct116_full_uce.h5ad",
                # Fallback to existing files
                'hek293t_alt': EMB_DIR / "uce/human_platelet_hek293t_uce_adata.h5ad",
                'hct116_alt': EMB_DIR / "uce/human_platelet_hct116_uce_adata.h5ad",
                # 33-layer versions
                'covid_33': EMB_DIR / "uce_33layer/covid_uce_33layer.h5ad",
                'hek293t_33': EMB_DIR / "uce_33layer/hek293t_uce_33layer.h5ad",
                'hct116_33': EMB_DIR / "uce_33layer/hct116_uce_33layer.h5ad",
                'obsm_keys': ['X_uce', 'X_uce_33layer', 'X_uce_4layer'],
            }

    def _find_file(self, key):
        """Find an existing file for a given key."""
        # For perturbation datasets, prefer X-Atlas files
        if key in ['hek293t', 'hct116']:
            # Try X-Atlas versions first (500k subset, then full)
            for xatlas_suffix in ['xatlas_500k', 'xatlas_full']:
                xatlas_key = f'{key}_{xatlas_suffix}'
                xatlas_path = self.paths.get(xatlas_key)
                if xatlas_path and xatlas_path.exists():
                    logger.info(f"  Found X-Atlas embeddings: {xatlas_path.name}")
                    return xatlas_path

        # Try the direct key
        path = self.paths.get(key)
        if path and path.exists():
            return path

        # Try alternatives
        for alt_key in [f'{key}_alt', f'{key}_33']:
            alt_path = self.paths.get(alt_key)
            if alt_path and alt_path.exists():
                return alt_path
        return None

    def _get_embeddings(self, adata):
        """Extract embeddings from adata."""
        for key in self.paths.get('obsm_keys', []):
            if key in adata.obsm:
                return np.array(adata.obsm[key])
        # Try any key with model name
        for key in adata.obsm.keys():
            if self.model.lower() in key.lower():
                return np.array(adata.obsm[key])
        # Fallback to X if dense
        if hasattr(adata.X, 'toarray'):
            return adata.X.toarray()
        return np.array(adata.X)

    def _get_severity(self, adata):
        """Extract severity labels."""
        for col in ['covid_severity', 'severity', 'Severity']:
            if col in adata.obs.columns:
                return adata.obs[col].astype(str).str.lower().values
        return np.array(['unknown'] * len(adata))

    def _get_perturbations(self, adata):
        """Extract perturbation labels."""
        for col in ['guide_target', 'gene_target', 'perturbation', 'target_gene']:
            if col in adata.obs.columns:
                return adata.obs[col].astype(str).values
        return np.array(['unknown'] * len(adata))

    def load_data(self):
        """Load COVID and perturbation data."""
        logger.info(f"Loading {self.model} data...")

        # Load COVID data
        covid_path = self._find_file('covid')
        if covid_path:
            logger.info(f"  COVID: {covid_path}")
            adata = sc.read_h5ad(covid_path)
            embeddings = self._get_embeddings(adata)
            severity = self._get_severity(adata)

            # Filter valid labels
            valid_mask = ~np.isin(severity, ['unknown', 'nan', ''])
            self.data['covid'] = {
                'embeddings': embeddings[valid_mask],
                'severity': severity[valid_mask],
                'adata': adata[valid_mask].copy()
            }
            logger.info(f"    Cells: {self.data['covid']['embeddings'].shape[0]:,}")
            logger.info(f"    Severity distribution: {pd.Series(self.data['covid']['severity']).value_counts().to_dict()}")
        else:
            logger.error("COVID data not found!")
            return False

        # Load perturbation data
        for dataset in ['hek293t', 'hct116']:
            path = self._find_file(dataset)
            if path:
                logger.info(f"  {dataset.upper()}: {path}")
                adata = sc.read_h5ad(path)
                embeddings = self._get_embeddings(adata)
                perturbations = self._get_perturbations(adata)

                self.data[dataset] = {
                    'embeddings': embeddings,
                    'perturbations': perturbations,
                    'adata': adata.copy()
                }
                n_perts = len(np.unique(perturbations))
                logger.info(f"    Cells: {embeddings.shape[0]:,}")
                logger.info(f"    Unique perturbations: {n_perts:,}")
            else:
                logger.warning(f"  {dataset.upper()}: NOT FOUND")

        return 'covid' in self.data and len(self.data) > 1

    def compute_severity_space(self):
        """Compute severity centroids and direction vectors."""
        logger.info("Computing severity space...")

        embeddings = self.data['covid']['embeddings']
        severity = self.data['covid']['severity']

        # Scale embeddings
        self.scaler = StandardScaler()
        embeddings_scaled = self.scaler.fit_transform(embeddings)

        # Compute centroids
        for sev in ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']:
            mask = severity == sev
            if mask.sum() > 0:
                self.severity_centroids[sev] = embeddings_scaled[mask].mean(axis=0)
                logger.info(f"  {sev}: {mask.sum():,} cells")

        # Recovery direction: severe → recovered
        if 'severe' in self.severity_centroids and 'recovered' in self.severity_centroids:
            self.recovery_direction = (
                self.severity_centroids['recovered'] - self.severity_centroids['severe']
            )
            self.recovery_direction = self.recovery_direction / np.linalg.norm(self.recovery_direction)
            logger.info("  Recovery direction computed (severe → recovered)")

        # Disease direction: healthy → severe
        if 'healthy' in self.severity_centroids and 'severe' in self.severity_centroids:
            self.disease_direction = (
                self.severity_centroids['severe'] - self.severity_centroids['healthy']
            )
            self.disease_direction = self.disease_direction / np.linalg.norm(self.disease_direction)
            logger.info("  Disease direction computed (healthy → severe)")

    def compute_perturbation_effects(self, dataset='hek293t', min_cells=10):
        """Compute perturbation effects as shifts in embedding space."""
        if dataset not in self.data:
            return pd.DataFrame()

        logger.info(f"Computing perturbation effects for {dataset.upper()}...")

        embeddings = self.data[dataset]['embeddings']
        perturbations = self.data[dataset]['perturbations']

        # Scale using COVID scaler
        embeddings_scaled = self.scaler.transform(embeddings)

        # Find control cells
        control_names = ['non-targeting', 'control', 'nt', 'negative', 'negctrl']
        control_mask = np.array([
            any(c in str(p).lower() for c in control_names)
            for p in perturbations
        ])

        if control_mask.sum() == 0:
            logger.warning("  No control cells found!")
            control_centroid = embeddings_scaled.mean(axis=0)
        else:
            control_centroid = embeddings_scaled[control_mask].mean(axis=0)
            logger.info(f"  Control cells: {control_mask.sum():,}")

        # Compute effects for each perturbation
        results = []
        unique_perts = np.unique(perturbations)

        for pert in unique_perts:
            # Skip controls
            if any(c in str(pert).lower() for c in control_names):
                continue

            mask = perturbations == pert
            n_cells = mask.sum()

            if n_cells < min_cells:
                continue

            # Perturbation centroid
            pert_centroid = embeddings_scaled[mask].mean(axis=0)

            # Effect = shift from control
            effect = pert_centroid - control_centroid

            # Recovery score = projection onto recovery direction
            recovery_score = np.dot(effect, self.recovery_direction) if self.recovery_direction is not None else 0

            # Disease score = projection onto disease direction
            disease_score = np.dot(effect, self.disease_direction) if self.disease_direction is not None else 0

            # Distance to severity centroids
            dist_to_centroids = {}
            for sev, centroid in self.severity_centroids.items():
                dist_to_centroids[f'dist_{sev}'] = np.linalg.norm(pert_centroid - centroid)

            # Closest severity
            closest_sev = min(dist_to_centroids.items(), key=lambda x: x[1])[0].replace('dist_', '')

            results.append({
                'perturbation': pert,
                'dataset': dataset.upper(),
                'n_cells': n_cells,
                'recovery_score': recovery_score,
                'disease_score': disease_score,
                'closest_severity': closest_sev,
                'effect_magnitude': np.linalg.norm(effect),
                **dist_to_centroids
            })

        df = pd.DataFrame(results)
        logger.info(f"  Computed effects for {len(df):,} perturbations")
        return df

    def find_therapeutic_candidates(self):
        """Find top therapeutic candidates across all datasets."""
        logger.info("\n" + "=" * 60)
        logger.info("FINDING THERAPEUTIC CANDIDATES")
        logger.info("=" * 60)

        all_results = []
        for dataset in ['hek293t', 'hct116']:
            if dataset in self.data:
                df = self.compute_perturbation_effects(dataset)
                if len(df) > 0:
                    all_results.append(df)

        if not all_results:
            logger.error("No perturbation results!")
            return None, None, None

        results_df = pd.concat(all_results, ignore_index=True)

        # Sort by recovery score
        results_df = results_df.sort_values('recovery_score', ascending=False)

        # Top therapeutic candidates (positive recovery score, close to healthy/recovered)
        therapeutic = results_df[
            (results_df['recovery_score'] > 0) &
            (results_df['closest_severity'].isin(['healthy', 'recovered', 'mild']))
        ].head(100)

        # Harmful perturbations (negative recovery score, close to severe/fatal)
        harmful = results_df[
            (results_df['recovery_score'] < 0) &
            (results_df['closest_severity'].isin(['severe', 'fatal']))
        ].head(100)

        logger.info(f"\nTotal perturbations analyzed: {len(results_df):,}")
        logger.info(f"Therapeutic candidates: {len(therapeutic):,}")
        logger.info(f"Harmful perturbations: {len(harmful):,}")

        logger.info("\n--- TOP 20 THERAPEUTIC CANDIDATES ---")
        for i, row in therapeutic.head(20).iterrows():
            logger.info(f"  {row['perturbation']}: recovery={row['recovery_score']:.4f}, "
                       f"cells={row['n_cells']}, closest={row['closest_severity']}")

        logger.info("\n--- TOP 20 HARMFUL PERTURBATIONS ---")
        for i, row in harmful.head(20).iterrows():
            logger.info(f"  {row['perturbation']}: recovery={row['recovery_score']:.4f}, "
                       f"cells={row['n_cells']}, closest={row['closest_severity']}")

        return results_df, therapeutic, harmful

    def plot_results(self, results_df, therapeutic, harmful):
        """Create visualization plots."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Plot 1: Recovery score distribution
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        ax = axes[0]
        ax.hist(results_df['recovery_score'], bins=50, alpha=0.7, edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', label='Zero (neutral)')
        ax.set_xlabel('Recovery Score')
        ax.set_ylabel('Number of Perturbations')
        ax.set_title(f'{self.model}: Recovery Score Distribution')
        ax.legend()

        ax = axes[1]
        colors = [SEVERITY_COLORS.get(s, 'gray') for s in results_df['closest_severity']]
        ax.scatter(results_df['recovery_score'], results_df['disease_score'],
                  c=colors, alpha=0.5, s=10)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Recovery Score')
        ax.set_ylabel('Disease Score')
        ax.set_title(f'{self.model}: Perturbation Effect Space')

        plt.tight_layout()
        plt.savefig(RESULTS_DIR / f'{self.model.lower()}_recovery_distribution_{timestamp}.png', dpi=300)
        plt.close()

        # Plot 2: Top candidates bar chart
        fig, axes = plt.subplots(1, 2, figsize=(16, 10))

        # Top therapeutic
        ax = axes[0]
        top_ther = therapeutic.head(30)
        colors = [SEVERITY_COLORS.get(s, 'gray') for s in top_ther['closest_severity']]
        y_pos = np.arange(len(top_ther))
        ax.barh(y_pos, top_ther['recovery_score'], color=colors, alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_ther['perturbation'], fontsize=8)
        ax.set_xlabel('Recovery Score')
        ax.set_title(f'{self.model}: Top 30 Therapeutic Candidates')
        ax.invert_yaxis()

        # Top harmful
        ax = axes[1]
        top_harm = harmful.head(30)
        colors = [SEVERITY_COLORS.get(s, 'gray') for s in top_harm['closest_severity']]
        y_pos = np.arange(len(top_harm))
        ax.barh(y_pos, top_harm['recovery_score'], color=colors, alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_harm['perturbation'], fontsize=8)
        ax.set_xlabel('Recovery Score')
        ax.set_title(f'{self.model}: Top 30 Harmful Perturbations')
        ax.invert_yaxis()

        plt.tight_layout()
        plt.savefig(RESULTS_DIR / f'{self.model.lower()}_therapeutic_rankings_{timestamp}.png', dpi=300)
        plt.close()

        logger.info(f"Plots saved to {RESULTS_DIR}")

    def save_results(self, results_df, therapeutic, harmful):
        """Save all results."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Full results
        results_df.to_csv(RESULTS_DIR / f'{self.model.lower()}_all_perturbations_{timestamp}.csv', index=False)

        # Therapeutic candidates
        if therapeutic is not None:
            therapeutic.to_csv(RESULTS_DIR / f'{self.model.lower()}_therapeutic_candidates_{timestamp}.csv', index=False)

        # Harmful
        if harmful is not None:
            harmful.to_csv(RESULTS_DIR / f'{self.model.lower()}_harmful_perturbations_{timestamp}.csv', index=False)

        # Summary
        with open(RESULTS_DIR / f'{self.model.lower()}_summary_{timestamp}.txt', 'w') as f:
            f.write(f"PERTURBATION RECOVERY ANALYSIS\n")
            f.write(f"Model: {self.model}\n")
            f.write(f"Date: {timestamp}\n")
            f.write("=" * 60 + "\n\n")

            f.write("DATA SUMMARY:\n")
            for key, data in self.data.items():
                n_cells = data['embeddings'].shape[0]
                f.write(f"  {key}: {n_cells:,} cells\n")

            f.write(f"\nTOTAL PERTURBATIONS: {len(results_df):,}\n")
            f.write(f"THERAPEUTIC CANDIDATES: {len(therapeutic):,}\n")
            f.write(f"HARMFUL PERTURBATIONS: {len(harmful):,}\n")

            f.write("\n\nTOP 50 THERAPEUTIC CANDIDATES:\n")
            f.write("-" * 60 + "\n")
            for _, row in therapeutic.head(50).iterrows():
                f.write(f"{row['perturbation']}: recovery={row['recovery_score']:.4f}, "
                       f"cells={row['n_cells']}, closest={row['closest_severity']}\n")

        logger.info(f"Results saved to {RESULTS_DIR}")

    def run(self):
        """Run full analysis."""
        logger.info("=" * 80)
        logger.info(f"PERTURBATION RECOVERY ANALYSIS: {self.model}")
        logger.info("=" * 80)

        if not self.load_data():
            logger.error("Failed to load data")
            return

        self.compute_severity_space()
        results_df, therapeutic, harmful = self.find_therapeutic_candidates()

        if results_df is not None:
            self.plot_results(results_df, therapeutic, harmful)
            self.save_results(results_df, therapeutic, harmful)

        logger.info("\n" + "=" * 80)
        logger.info("ANALYSIS COMPLETE")
        logger.info("=" * 80)


def main():
    parser = argparse.ArgumentParser(description='Perturbation Recovery Analysis')
    parser.add_argument('--model', type=str, default='STATE',
                        choices=['STATE', 'UCE'],
                        help='Foundation model to analyze')
    args = parser.parse_args()

    analysis = PerturbationRecoveryAnalysis(model=args.model)
    analysis.run()


if __name__ == "__main__":
    main()
