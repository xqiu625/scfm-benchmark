#!/usr/bin/env python3
"""
Differential Diagnosis of Inflammatory Conditions (6-class)

Task 2 of the Two-Task Clinical Framework:
Can foundation models distinguish COVID-19 from other inflammatory conditions?

6-Class Problem:
- Healthy
- COVID-19
- Sepsis
- SLE (Systemic Lupus Erythematosus)
- Influenza (Flu)
- Pneumonia

Outputs:
- Per-model macro AUC
- Per-disease AUC (One-vs-Rest)
- 6x6 Confusion matrices
- Disease-specific precision/recall

Author: Project 302006
Date: January 2026
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score, confusion_matrix, classification_report,
    precision_recall_fscore_support
)
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import json
import logging
from joblib import Parallel, delayed

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DifferentialDiagnosisBenchmark:
    """6-class differential diagnosis benchmark for inflammatory conditions"""

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_dir = Path(base_dir)
        self.embedding_dir = self.base_dir / "02_EMBEDDINGS"
        self.output_dir = self.base_dir / "05_RESULTS/15_differential_diagnosis"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Also save to analyses folder structure
        self.analyses_dir = self.base_dir / "scripts/15_differential_diagnosis"
        self.analyses_results_dir = self.analyses_dir / "results"
        self.analyses_figures_dir = self.analyses_dir / "figures"
        self.analyses_results_dir.mkdir(parents=True, exist_ok=True)
        self.analyses_figures_dir.mkdir(parents=True, exist_ok=True)

        # Disease class mapping - standardize names
        # Based on raw data: COVID-19, Sepsis, HC, Lung_infect, Non_covid, SLE, Flu
        self.disease_classes = {
            # Healthy
            'healthy': 'Healthy',
            'control': 'Healthy',
            'hc': 'Healthy',  # HC in raw data
            # COVID-19
            'covid': 'COVID-19',
            'covid-19': 'COVID-19',
            'covid19': 'COVID-19',
            # Sepsis
            'sepsis': 'Sepsis',
            # SLE
            'sle': 'SLE',
            'lupus': 'SLE',
            # Influenza
            'flu': 'Influenza',
            'influenza': 'Influenza',
            # Pneumonia / Lung infection
            'pneumonia': 'Pneumonia',
            'lung_infect': 'Pneumonia',
            'lung_infection': 'Pneumonia',
            'bacterial_pneumonia': 'Pneumonia',
            # Non-COVID (exclude from analysis)
            'non_covid': 'Non_COVID',
        }

        # Target classes for differential diagnosis (5-class)
        # Raw data: COVID-19 (38673), Sepsis (3422), HC (3205), Lung_infect (1749), Non_covid (726), SLE (169), Flu (33)
        # Excluded: Non_COVID (catch-all), Influenza (n=33 too small)
        self.target_classes_full = ['Healthy', 'COVID-19', 'Sepsis', 'SLE', 'Pneumonia']
        self.target_classes = None  # Will be set based on available data
        self.min_samples_per_class = 10  # Minimum samples needed per class

        # Models to evaluate
        self.models = ['STATE', 'UCE', 'Geneformer', 'scGPT', 'TranscriptFormer', 'scPRINT2']

        # Embedding configurations - actual file names on HPCC
        self.embedding_configs = {
            'STATE': {
                'files': [
                    self.embedding_dir / "state/covid_state.h5ad",
                ],
                'obsm_keys': ['X_emb', 'X_state', 'X_STATE', 'state'],
            },
            'UCE': {
                'files': [
                    self.embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad",
                ],
                'obsm_keys': ['X_uce', 'X_UCE', 'X_uce_4layer'],
            },
            'scGPT': {
                'files': [
                    self.embedding_dir / "scgpt/covid_scgpt.h5ad",
                    self.base_dir / "04_ANALYSIS/scgpt/covid_scgpt.h5ad",
                ],
                'obsm_keys': ['X_scGPT', 'X_scgpt', 'scgpt'],
            },
            'TranscriptFormer': {
                'files': [
                    self.embedding_dir / "transcriptformer/covid_transcriptformer.h5ad",
                ],
                'obsm_keys': ['X_transcriptformer', 'X_TranscriptFormer', 'transcriptformer'],
            },
            'Geneformer': {
                'files': [
                    self.embedding_dir / "geneformer/covid_geneformer.h5ad",
                ],
                'obsm_keys': ['X_geneformer', 'X_Geneformer', 'geneformer'],
            },
            'scPRINT2': {
                'files': [
                    self.embedding_dir / "scprint2/covid_scprint2.h5ad",
                ],
                'obsm_keys': ['X_scprint2', 'scprint_emb'],
            }
        }

        self.results = {}

    def find_embedding_file(self, model_name):
        """Find the embedding file for a model"""
        config = self.embedding_configs.get(model_name)
        if not config:
            return None

        for fp in config['files']:
            if fp.exists():
                return fp
        return None

    def load_embeddings(self, model_name, adata):
        """Extract embeddings from adata for a specific model"""
        config = self.embedding_configs.get(model_name)
        if not config:
            return None

        for key in config['obsm_keys']:
            if key in adata.obsm:
                embeddings = adata.obsm[key]
                if hasattr(embeddings, 'toarray'):
                    embeddings = embeddings.toarray()
                return np.array(embeddings)

        # Try to find any matching key
        possible_keys = [k for k in adata.obsm.keys() if model_name.lower() in k.lower()]
        if possible_keys:
            embeddings = adata.obsm[possible_keys[0]]
            if hasattr(embeddings, 'toarray'):
                embeddings = embeddings.toarray()
            return np.array(embeddings)

        return None

    def standardize_disease_labels(self, labels):
        """Standardize disease labels to 6 target classes"""
        standardized = []
        for label in labels:
            label_lower = str(label).lower().strip()
            if label_lower in self.disease_classes:
                standardized.append(self.disease_classes[label_lower])
            else:
                standardized.append('Unknown')
        return np.array(standardized)

    def load_all_disease_data(self, model_name):
        """Load embeddings and disease labels for a model"""
        logger.info(f"Loading {model_name} embeddings...")

        file_path = self.find_embedding_file(model_name)
        if file_path is None:
            logger.warning(f"No embedding file found for {model_name}")
            return None, None, None

        logger.info(f"  File: {file_path}")
        adata = sc.read_h5ad(file_path)

        # Extract embeddings
        X = self.load_embeddings(model_name, adata)
        if X is None:
            logger.warning(f"  Could not extract embeddings for {model_name}")
            logger.info(f"  Available obsm keys: {list(adata.obsm.keys())}")
            return None, None, None

        # Find disease labels
        disease_cols = ['disease', 'disease_condition', 'condition', 'Disease', 'dataset']
        disease_col = None
        for col in disease_cols:
            if col in adata.obs.columns:
                disease_col = col
                break

        if disease_col is None:
            logger.warning(f"  No disease label column found for {model_name}")
            logger.info(f"  Available columns: {list(adata.obs.columns)}")
            return None, None, None

        # Standardize labels
        raw_labels = adata.obs[disease_col].values
        y = self.standardize_disease_labels(raw_labels)

        # Filter to target classes that have sufficient samples
        available_classes = []
        for cls in self.target_classes_full:
            count = np.sum(y == cls)
            if count >= self.min_samples_per_class:
                available_classes.append(cls)

        if len(available_classes) < 2:
            logger.warning(f"  Less than 2 classes with sufficient samples")
            return None, None, None

        # Update target classes for this run
        self.target_classes = available_classes

        mask = np.isin(y, available_classes)
        X_filtered = X[mask]
        y_filtered = y[mask]

        # Extract patient/sample-level groups for GroupKFold
        groups = None
        for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
            if col in adata.obs.columns:
                groups = adata.obs[col].values[mask]
                logger.info(f"  Using group column: {col} ({len(np.unique(groups))} groups)")
                break
        if groups is None:
            # External code-review finding (2026-05-24): silently falling back to
            # cell-index groups degrades the benchmark to cell-level CV without
            # disclosure. Fail loudly so the caller knows the input h5ad lacks
            # the required patient-level grouping column.
            raise RuntimeError(
                f"No patient-level group column found in input "
                f"(tried Sample_ID, sample_id, Data_NO, data_no, batch, donor). "
                f"Available obs columns: {list(adata.obs.columns)[:20]}. "
                f"Patient-level cross-validation is required for publication-grade "
                f"differential-diagnosis results; using cell index as a group label "
                f"would silently revert to cell-level CV and produce inflated AUCs."
            )

        logger.info(f"  Shape: {X_filtered.shape}")
        logger.info(f"  Available classes ({len(available_classes)}):")
        for cls in available_classes:
            count = np.sum(y_filtered == cls)
            logger.info(f"    {cls}: {count}")

        return X_filtered, y_filtered, adata.obs.index[mask].tolist(), groups

    def evaluate_model(self, X, y, n_splits=5, n_jobs=8, use_pca=None, groups=None):
        """
        Evaluate 6-class differential diagnosis for a single model

        Parameters
        ----------
        X : np.ndarray
            Embedding matrix (n_cells, n_dims)
        y : np.ndarray
            Disease labels
        n_splits : int
            Number of CV folds
        n_jobs : int
            Parallel jobs for CV
        use_pca : int or None
            If set, apply PCA reduction first (for scGPT denoising)
        groups : np.ndarray or None
            Patient/sample group labels for GroupKFold

        Returns
        -------
        dict with results
        """
        # Encode labels - use available classes
        le = LabelEncoder()
        le.fit(self.target_classes if self.target_classes else np.unique(y))
        y_encoded = le.transform(y)
        n_classes = len(le.classes_)

        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)

        all_y_true = []
        all_y_pred = []
        all_y_pred_proba = []
        fold_indices = []

        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y_encoded, groups=groups)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

            # Scale and train
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)

            # Optional PCA (fit on train fold only to prevent data leakage)
            if use_pca is not None:
                pca = PCA(n_components=min(use_pca, X_train_scaled.shape[1], X_train_scaled.shape[0]))
                X_train_scaled = pca.fit_transform(X_train_scaled)
                X_test_scaled = pca.transform(X_test_scaled)
                if fold_idx == 0:
                    logger.info(f"  Applied PCA: {pca.n_components_} components (fit per fold)")

            clf = LogisticRegression(
                solver='lbfgs',
                max_iter=1000,
                class_weight='balanced',
                multi_class='multinomial',
                random_state=42
            )
            clf.fit(X_train_scaled, y_train)

            y_pred = clf.predict(X_test_scaled)
            y_pred_proba = clf.predict_proba(X_test_scaled)

            all_y_true.extend(y_test)
            all_y_pred.extend(y_pred)
            all_y_pred_proba.append(y_pred_proba)
            fold_indices.extend([fold_idx] * len(y_test))

        all_y_true = np.array(all_y_true)
        all_y_pred = np.array(all_y_pred)
        all_y_pred_proba = np.vstack(all_y_pred_proba)

        # Calculate per-class AUC (One-vs-Rest)
        class_aucs = {}
        for i, class_name in enumerate(le.classes_):
            y_binary = (all_y_true == i).astype(int)
            if len(np.unique(y_binary)) < 2:
                class_aucs[class_name] = np.nan
            else:
                class_aucs[class_name] = roc_auc_score(y_binary, all_y_pred_proba[:, i])

        # Macro-averaged AUC
        try:
            macro_auc = roc_auc_score(all_y_true, all_y_pred_proba,
                                       multi_class='ovr', average='macro')
        except Exception:
            macro_auc = np.nan

        # Weighted AUC
        try:
            weighted_auc = roc_auc_score(all_y_true, all_y_pred_proba,
                                          multi_class='ovr', average='weighted')
        except Exception:
            weighted_auc = np.nan

        # Confusion matrix
        cm = confusion_matrix(all_y_true, all_y_pred)

        # Normalized confusion matrix (row-normalized = recall per class)
        cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

        # Per-class precision, recall, f1
        precision, recall, f1, support = precision_recall_fscore_support(
            all_y_true, all_y_pred, average=None
        )

        class_metrics = {}
        for i, class_name in enumerate(le.classes_):
            class_metrics[class_name] = {
                'auc': class_aucs[class_name],
                'precision': precision[i],
                'recall': recall[i],
                'f1': f1[i],
                'support': int(support[i])
            }

        return {
            'macro_auc': macro_auc,
            'weighted_auc': weighted_auc,
            'class_aucs': class_aucs,
            'class_metrics': class_metrics,
            'confusion_matrix': cm,
            'confusion_matrix_normalized': cm_normalized,
            'class_names': le.classes_.tolist(),
            'n_samples': len(all_y_true),
            'predictions': {
                'y_true': all_y_true.tolist(),
                'y_pred': all_y_pred.tolist(),
                'y_pred_proba': all_y_pred_proba.tolist()
            }
        }

    def run_benchmark(self, n_splits=5, n_jobs=8):
        """Run the full differential diagnosis benchmark"""
        logger.info("="*80)
        logger.info("DIFFERENTIAL DIAGNOSIS BENCHMARK (6-CLASS)")
        logger.info("Classes: Healthy, COVID-19, Sepsis, SLE, Influenza, Pneumonia")
        logger.info(f"CV Folds: {n_splits}")
        logger.info("="*80)

        results = {}

        for model_name in self.models:
            logger.info(f"\n{'='*60}")
            logger.info(f"Model: {model_name}")
            logger.info(f"{'='*60}")

            X, y, cell_ids, groups = self.load_all_disease_data(model_name)
            if X is None:
                logger.warning(f"Skipping {model_name} - data not available")
                continue

            # Check for minimum samples per class for CV
            min_samples = min([np.sum(y == cls) for cls in self.target_classes])
            if min_samples < n_splits:
                logger.warning(f"  Insufficient samples for {n_splits}-fold CV (min={min_samples})")
                logger.info(f"  Reducing to {min_samples}-fold CV")
                n_splits_model = max(2, min_samples)
            else:
                n_splits_model = n_splits

            # Run evaluation with patient-level grouping
            result = self.evaluate_model(X, y, n_splits=n_splits_model, n_jobs=n_jobs, groups=groups)
            result['n_classes'] = len(self.target_classes)
            result['classes_used'] = self.target_classes.copy()
            results[model_name] = result

            logger.info(f"\n  Macro AUC: {result['macro_auc']:.3f}")
            logger.info(f"  Per-disease AUC:")
            for disease, auc in result['class_aucs'].items():
                auc_str = f"{auc:.3f}" if not np.isnan(auc) else "N/A"
                logger.info(f"    {disease}: {auc_str}")

            # If scGPT, also run with PCA denoising
            if model_name == 'scGPT':
                logger.info(f"\n  Running with PCA-256 denoising...")
                result_pca = self.evaluate_model(X, y, n_splits=n_splits_model,
                                                  n_jobs=n_jobs, use_pca=256, groups=groups)
                result_pca['n_classes'] = len(self.target_classes)
                result_pca['classes_used'] = self.target_classes.copy()
                results[f'{model_name}_PCA256'] = result_pca
                logger.info(f"  scGPT+PCA256 Macro AUC: {result_pca['macro_auc']:.3f}")

        self.results = results
        return results

    def create_summary_table(self):
        """Create summary table of all results"""
        rows = []

        for model_name, result in self.results.items():
            row = {
                'Model': model_name,
                'Macro_AUC': result['macro_auc'],
                'Weighted_AUC': result['weighted_auc'],
                'N_Samples': result['n_samples']
            }

            # Add per-disease AUC
            for disease, auc in result['class_aucs'].items():
                row[f'AUC_{disease}'] = auc

            rows.append(row)

        df = pd.DataFrame(rows)
        df = df.sort_values('Macro_AUC', ascending=False)
        return df

    def plot_confusion_matrices(self, save=True):
        """Plot confusion matrices for all models"""
        n_models = len(self.results)
        if n_models == 0:
            return

        # Determine grid layout
        n_cols = min(3, n_models)
        n_rows = (n_models + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_models == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        for idx, (model_name, result) in enumerate(self.results.items()):
            ax = axes[idx]
            cm_norm = result['confusion_matrix_normalized']
            class_names = result['class_names']

            # Plot heatmap
            sns.heatmap(cm_norm, annot=True, fmt='.2f', cmap='Blues',
                       xticklabels=class_names, yticklabels=class_names,
                       ax=ax, vmin=0, vmax=1, cbar_kws={'label': 'Rate'})

            ax.set_title(f'{model_name}\nMacro AUC: {result["macro_auc"]:.3f}')
            ax.set_xlabel('Predicted')
            ax.set_ylabel('True')

            # Rotate labels for readability
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

        # Hide unused subplots
        for idx in range(len(self.results), len(axes)):
            axes[idx].set_visible(False)

        plt.tight_layout()

        if save:
            # Save to both output locations
            for output_dir in [self.output_dir, self.analyses_figures_dir]:
                plt.savefig(output_dir / 'fig_15_confusion_matrices_all.png',
                           dpi=150, bbox_inches='tight')
                plt.savefig(output_dir / 'fig_15_confusion_matrices_all.pdf',
                           bbox_inches='tight')

        plt.close()
        logger.info("Saved confusion matrix figures")

    def plot_best_model_confusion(self, save=True):
        """Plot confusion matrix for best model only (for main figure)"""
        if not self.results:
            return

        # Find best model
        best_model = max(self.results.items(),
                        key=lambda x: x[1]['macro_auc'] if not np.isnan(x[1]['macro_auc']) else 0)
        model_name, result = best_model

        fig, ax = plt.subplots(figsize=(8, 7))

        cm = result['confusion_matrix']
        cm_norm = result['confusion_matrix_normalized']
        class_names = result['class_names']

        # Create annotation with both count and percentage
        annot = []
        for i in range(len(class_names)):
            row = []
            for j in range(len(class_names)):
                row.append(f'{cm[i,j]}\n({cm_norm[i,j]:.1%})')
            annot.append(row)
        annot = np.array(annot)

        sns.heatmap(cm_norm, annot=annot, fmt='', cmap='Blues',
                   xticklabels=class_names, yticklabels=class_names,
                   ax=ax, vmin=0, vmax=1, cbar_kws={'label': 'Classification Rate'})

        ax.set_title(f'Differential Diagnosis Confusion Matrix\n{model_name} (Macro AUC: {result["macro_auc"]:.3f})',
                    fontsize=14)
        ax.set_xlabel('Predicted Diagnosis', fontsize=12)
        ax.set_ylabel('True Diagnosis', fontsize=12)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

        plt.tight_layout()

        if save:
            for output_dir in [self.output_dir, self.analyses_figures_dir]:
                plt.savefig(output_dir / f'fig_02_confusion_matrix_{model_name}.png',
                           dpi=300, bbox_inches='tight')
                plt.savefig(output_dir / f'fig_02_confusion_matrix_{model_name}.pdf',
                           bbox_inches='tight')

        plt.close()
        logger.info(f"Saved best model ({model_name}) confusion matrix")

    def plot_per_disease_auc(self, save=True):
        """Plot per-disease AUC comparison across models"""
        if not self.results:
            return

        # Prepare data
        models = []
        diseases = self.target_classes
        auc_data = {disease: [] for disease in diseases}

        for model_name, result in self.results.items():
            if '_PCA' in model_name:  # Skip PCA variants for main plot
                continue
            models.append(model_name)
            for disease in diseases:
                auc_data[disease].append(result['class_aucs'].get(disease, np.nan))

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=(12, 6))

        x = np.arange(len(models))
        width = 0.12
        n_diseases = len(diseases)

        colors = plt.cm.Set2(np.linspace(0, 1, n_diseases))

        for i, disease in enumerate(diseases):
            offset = (i - n_diseases/2 + 0.5) * width
            bars = ax.bar(x + offset, auc_data[disease], width,
                         label=disease, color=colors[i])

        ax.set_xlabel('Model', fontsize=12)
        ax.set_ylabel('AUC (One-vs-Rest)', fontsize=12)
        ax.set_title('Per-Disease Classification Performance', fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(models, rotation=45, ha='right')
        ax.legend(title='Disease', bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.set_ylim(0.5, 1.0)
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Random')
        ax.grid(axis='y', alpha=0.3)

        plt.tight_layout()

        if save:
            for output_dir in [self.output_dir, self.analyses_figures_dir]:
                plt.savefig(output_dir / 'fig_15_per_disease_auc.png',
                           dpi=150, bbox_inches='tight')
                plt.savefig(output_dir / 'fig_15_per_disease_auc.pdf',
                           bbox_inches='tight')

        plt.close()
        logger.info("Saved per-disease AUC figure")

    def plot_model_ranking(self, save=True):
        """Plot model ranking by macro AUC"""
        if not self.results:
            return

        # Get macro AUC for each model (exclude PCA variants)
        model_aucs = {name: res['macro_auc']
                     for name, res in self.results.items()
                     if '_PCA' not in name}

        # Sort by AUC
        sorted_models = sorted(model_aucs.items(), key=lambda x: x[1], reverse=True)
        models = [m[0] for m in sorted_models]
        aucs = [m[1] for m in sorted_models]

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = plt.cm.viridis(np.linspace(0.8, 0.2, len(models)))
        bars = ax.barh(models, aucs, color=colors)

        ax.set_xlabel('Macro AUC', fontsize=12)
        ax.set_title('Model Ranking: 6-Class Differential Diagnosis', fontsize=14)
        ax.set_xlim(0.5, 1.0)
        ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)

        # Add value labels
        for bar, auc in zip(bars, aucs):
            ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                   f'{auc:.3f}', va='center', fontsize=10)

        plt.tight_layout()

        if save:
            for output_dir in [self.output_dir, self.analyses_figures_dir]:
                plt.savefig(output_dir / 'fig_15_model_ranking.png',
                           dpi=150, bbox_inches='tight')
                plt.savefig(output_dir / 'fig_15_model_ranking.pdf',
                           bbox_inches='tight')

        plt.close()
        logger.info("Saved model ranking figure")

    def save_results(self):
        """Save all results to files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Save to both output directories
        for output_dir in [self.output_dir, self.analyses_results_dir]:
            # Save detailed JSON
            json_file = output_dir / f'differential_diagnosis_results_{timestamp}.json'

            # Prepare serializable results
            serializable_results = {}
            for model_name, result in self.results.items():
                serializable_results[model_name] = {
                    'macro_auc': float(result['macro_auc']) if not np.isnan(result['macro_auc']) else None,
                    'weighted_auc': float(result['weighted_auc']) if not np.isnan(result['weighted_auc']) else None,
                    'class_aucs': {k: float(v) if not np.isnan(v) else None
                                   for k, v in result['class_aucs'].items()},
                    'class_metrics': result['class_metrics'],
                    'confusion_matrix': result['confusion_matrix'].tolist(),
                    'confusion_matrix_normalized': result['confusion_matrix_normalized'].tolist(),
                    'class_names': result['class_names'],
                    'n_samples': result['n_samples']
                }

            with open(json_file, 'w') as f:
                json.dump(serializable_results, f, indent=2)

            # Save summary table
            df = self.create_summary_table()
            csv_file = output_dir / f'differential_diagnosis_summary_{timestamp}.csv'
            df.to_csv(csv_file, index=False)

            # Save confusion matrices as separate CSV files
            for model_name, result in self.results.items():
                cm_df = pd.DataFrame(
                    result['confusion_matrix'],
                    index=result['class_names'],
                    columns=result['class_names']
                )
                cm_file = output_dir / f'confusion_matrix_{model_name}.csv'
                cm_df.to_csv(cm_file)

            # Save best results text summary
            summary_file = output_dir / f'differential_diagnosis_summary_{timestamp}.txt'
            with open(summary_file, 'w') as f:
                f.write("DIFFERENTIAL DIAGNOSIS BENCHMARK RESULTS\n")
                f.write("="*60 + "\n")
                f.write(f"Task: 6-class classification\n")
                f.write(f"Classes: {', '.join(self.target_classes)}\n")
                f.write("="*60 + "\n\n")

                f.write("MODEL RANKING (by Macro AUC):\n")
                f.write("-"*40 + "\n")

                sorted_results = sorted(
                    [(name, res['macro_auc']) for name, res in self.results.items()
                     if '_PCA' not in name],
                    key=lambda x: x[1] if not np.isnan(x[1]) else 0,
                    reverse=True
                )

                for rank, (model, auc) in enumerate(sorted_results, 1):
                    f.write(f"  {rank}. {model}: {auc:.4f}\n")

                f.write("\n\nPER-DISEASE AUC:\n")
                f.write("-"*40 + "\n")

                for model_name, result in self.results.items():
                    if '_PCA' in model_name:
                        continue
                    f.write(f"\n{model_name}:\n")
                    for disease, auc in result['class_aucs'].items():
                        auc_str = f"{auc:.4f}" if not np.isnan(auc) else "N/A"
                        f.write(f"  {disease}: {auc_str}\n")

                f.write("\n\nKEY FINDINGS:\n")
                f.write("-"*40 + "\n")

                # Find best and worst performing disease
                best_model = sorted_results[0][0] if sorted_results else None
                if best_model and best_model in self.results:
                    class_aucs = self.results[best_model]['class_aucs']
                    best_disease = max(class_aucs.items(), key=lambda x: x[1] if not np.isnan(x[1]) else 0)
                    worst_disease = min(class_aucs.items(), key=lambda x: x[1] if not np.isnan(x[1]) else 1)
                    f.write(f"Best classified disease: {best_disease[0]} (AUC: {best_disease[1]:.3f})\n")
                    f.write(f"Most difficult disease: {worst_disease[0]} (AUC: {worst_disease[1]:.3f})\n")

        logger.info(f"Results saved to {self.output_dir} and {self.analyses_results_dir}")

    def run_all(self, n_splits=5, n_jobs=8):
        """Run complete analysis pipeline"""
        self.run_benchmark(n_splits=n_splits, n_jobs=n_jobs)

        if self.results:
            self.save_results()
            self.plot_confusion_matrices()
            self.plot_best_model_confusion()
            self.plot_per_disease_auc()
            self.plot_model_ranking()

            # Print summary
            df = self.create_summary_table()
            print("\n" + "="*80)
            print("DIFFERENTIAL DIAGNOSIS RESULTS SUMMARY")
            print("="*80)
            print(df.to_string(index=False))
            print("="*80)


def main():
    """Main execution"""
    import argparse
    parser = argparse.ArgumentParser(
        description='Differential Diagnosis Benchmark (6-class)'
    )
    parser.add_argument('--base-dir', type=str,
                       default="/bigdata/godziklab/shared/Xinru/302006",
                       help='Base directory for the project')
    parser.add_argument('--cv-folds', type=int, default=5,
                       help='Number of cross-validation folds')
    parser.add_argument('--n-jobs', type=int, default=8,
                       help='Number of parallel jobs')
    args = parser.parse_args()

    # Initialize and run benchmark
    benchmark = DifferentialDiagnosisBenchmark(base_dir=args.base_dir)
    benchmark.run_all(n_splits=args.cv_folds, n_jobs=args.n_jobs)

    return 0


if __name__ == "__main__":
    sys.exit(main())
