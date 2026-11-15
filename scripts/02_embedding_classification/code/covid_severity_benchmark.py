#!/usr/bin/env python3
"""
COVID Severity Prediction Benchmark

Compare 4 foundation model embeddings (UCE, scGPT, TranscriptFormer, Geneformer)
for predicting COVID-19 disease severity.

Classification tasks:
1. 6-class: healthy, mild, moderate, severe, fatal, recovered
2. 3-class: control, mild, severe
3. Binary: severe (severe+fatal) vs non-severe

Metrics: Accuracy, Balanced Accuracy, AUC (macro/weighted), F1-score

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

from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    roc_auc_score, classification_report, confusion_matrix
)
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import json
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class COVIDSeverityBenchmark:
    """Benchmark foundation model embeddings for COVID severity prediction"""

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_dir = Path(base_dir)
        self.embedding_dir = self.base_dir / "02_EMBEDDINGS"
        self.output_dir = self.base_dir / "05_RESULTS/covid_severity_benchmark"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Embedding configurations - try multiple possible file names
        self.embedding_configs = {
            'UCE': {
                'files': [
                    self.embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad",
                    self.embedding_dir / "uce/covid_uce.h5ad",
                ],
                'obsm_keys': ['X_uce', 'X_uce_4layer'],
                'dims': 1280
            },
            'scGPT': {
                'files': [
                    self.embedding_dir / "scgpt/covid_scgpt.h5ad",
                    self.embedding_dir / "scgpt/platelet_scgpt_combined.h5ad",
                ],
                'obsm_keys': ['X_scGPT', 'X_scgpt'],
                'dims': 512
            },
            'TranscriptFormer': {
                'files': [
                    self.embedding_dir / "transcriptformer/covid_transcriptformer.h5ad",
                ],
                'obsm_keys': ['X_transcriptformer', 'X_TranscriptFormer'],
                'dims': 2048
            },
            'Geneformer': {
                'files': [
                    self.embedding_dir / "geneformer/covid_geneformer.h5ad",
                ],
                'obsm_keys': ['X_geneformer', 'X_Geneformer'],
                'dims': 1152
            },
            'scPRINT2': {
                'files': [
                    self.embedding_dir / "scprint2/covid_scprint2.h5ad",
                ],
                'obsm_keys': ['X_scprint2', 'scprint_emb'],
                'dims': 256
            }
        }

        # Classification tasks
        self.tasks = {
            '6-class': None,  # Original 6 severity levels
            '3-class': {      # Simplified 3 groups
                'control': ['healthy'],
                'mild': ['mild', 'moderate', 'recovered'],
                'severe': ['severe', 'fatal']
            },
            'binary': {       # Binary severe vs non-severe
                'non_severe': ['healthy', 'mild', 'moderate', 'recovered'],
                'severe': ['severe', 'fatal']
            }
        }

        # Classifiers to evaluate
        self.classifiers = {
            'LogisticRegression': LogisticRegression(
                max_iter=2000, random_state=42, n_jobs=-1, class_weight='balanced',
                solver='saga'  # saga converges faster for large datasets
            ),
            'RandomForest': RandomForestClassifier(
                n_estimators=100, random_state=42, n_jobs=-1, class_weight='balanced'
            ),
            'SVM': SVC(
                kernel='rbf', probability=True, random_state=42, class_weight='balanced'
            ),
            'MLP': MLPClassifier(
                hidden_layer_sizes=(256, 128), max_iter=500, random_state=42,
                early_stopping=True, validation_fraction=0.1
            ),
            'GradientBoosting': GradientBoostingClassifier(
                n_estimators=100, random_state=42
            )
        }

        self.results = {}

    def load_embeddings(self):
        """Load all embedding files and extract COVID severity labels"""
        logger.info("Loading embeddings from all models...")
        self.embeddings = {}
        self.model_labels = {}  # Store labels per model
        self.model_cell_ids = {}  # Store cell IDs per model
        self.model_groups = {}  # Store patient/sample groups per model

        for model_name, config in self.embedding_configs.items():
            # Try multiple possible file paths
            file_path = None
            for fp in config['files']:
                if fp.exists():
                    file_path = fp
                    break

            if file_path is None:
                logger.warning(f"No file found for {model_name}")
                logger.warning(f"  Tried: {[str(f) for f in config['files']]}")
                continue

            logger.info(f"Loading {model_name} from {file_path}")
            adata = sc.read_h5ad(file_path)

            # For scGPT combined file, filter to COVID only (cells with valid severity labels)
            if model_name == 'scGPT' and 'platelet_scgpt_combined' in str(file_path):
                # Try filtering by dataset column first
                if 'dataset' in adata.obs.columns:
                    covid_mask = adata.obs['dataset'] == 'covid'
                    if covid_mask.sum() > 0:
                        adata = adata[covid_mask].copy()
                        logger.info(f"  Filtered by dataset==covid: {adata.n_obs} cells")

                # Also filter by valid severity labels
                if 'covid_severity' in adata.obs.columns:
                    valid_severities = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']
                    severity_labels = np.array([str(x) if pd.notna(x) else 'unknown' for x in adata.obs['covid_severity'].values])
                    valid_mask = np.isin(severity_labels, valid_severities)
                    if valid_mask.sum() < len(adata):
                        adata = adata[valid_mask].copy()
                        logger.info(f"  Filtered by valid severity: {adata.n_obs} cells")

            # Get embeddings - try multiple possible keys
            obsm_key = None
            for key in config['obsm_keys']:
                if key in adata.obsm:
                    obsm_key = key
                    break

            if obsm_key is None:
                # Try to find any key containing model name
                possible_keys = [k for k in adata.obsm.keys() if model_name.lower() in k.lower()]
                if possible_keys:
                    obsm_key = possible_keys[0]
                    logger.info(f"  Using discovered key: {obsm_key}")
                else:
                    logger.warning(f"  Embedding key not found for {model_name}")
                    logger.info(f"  Available keys: {list(adata.obsm.keys())}")
                    continue

            embeddings = adata.obsm[obsm_key]
            if hasattr(embeddings, 'toarray'):
                embeddings = embeddings.toarray()
            embeddings = np.array(embeddings)

            # Get labels for this model
            labels = None
            label_cols = ['covid_severity', 'severity', 'Category', 'three_group']
            for col in label_cols:
                if col in adata.obs.columns:
                    raw_labels = adata.obs[col].values
                    # Convert to string array and handle NaN/None
                    labels = np.array([str(x) if pd.notna(x) else 'unknown' for x in raw_labels])
                    logger.info(f"  Using label column: {col}")
                    break

            if labels is None:
                logger.warning(f"  No severity label column found for {model_name}")
                logger.info(f"  Available columns: {list(adata.obs.columns)}")
                continue

            self.embeddings[model_name] = embeddings
            self.model_labels[model_name] = labels
            self.model_cell_ids[model_name] = adata.obs.index.tolist()

            # Extract patient/sample-level groups for GroupKFold
            groups = None
            for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
                if col in adata.obs.columns:
                    groups = adata.obs[col].values
                    logger.info(f"  Using group column: {col} ({len(np.unique(groups))} groups)")
                    break
            if groups is None:
                logger.warning("  No group column found — using cell index as group")
                groups = np.arange(adata.n_obs)
            self.model_groups[model_name] = groups

            logger.info(f"  Shape: {embeddings.shape}")
            logger.info(f"  Obsm key: {obsm_key}")
            logger.info(f"  Label distribution: {pd.Series(labels).value_counts().to_dict()}")

        logger.info(f"\nLoaded {len(self.embeddings)} embedding models")

        return len(self.embeddings) > 0

    def prepare_labels(self, task_name, labels):
        """Prepare labels for a specific classification task"""
        # Valid severity labels (exclude unknown/nan)
        valid_severities = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']

        if task_name == '6-class':
            # For 6-class, return labels and mask for valid severity labels only
            valid_mask = np.isin(labels, valid_severities)
            return labels, valid_mask

        task_mapping = self.tasks[task_name]
        new_labels = np.array(['unknown'] * len(labels), dtype=object)

        for new_label, old_labels in task_mapping.items():
            for old_label in old_labels:
                mask = labels == old_label
                new_labels[mask] = new_label

        # Remove unknown labels
        valid_mask = new_labels != 'unknown'
        return new_labels, valid_mask

    def evaluate_classifier(self, X, y, classifier, n_splits=5, groups=None):
        """Evaluate a classifier using cross-validation"""
        # Encode labels
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        n_classes = len(le.classes_)

        # Create pipeline with scaling
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', classifier)
        ])

        # Cross-validation with patient-level grouping
        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)

        # Get predictions
        y_pred = cross_val_predict(pipeline, X, y_encoded, cv=cv, groups=groups)

        # Get probabilities for AUC (if supported)
        try:
            y_prob = cross_val_predict(pipeline, X, y_encoded, cv=cv, groups=groups, method='predict_proba')

            # Calculate AUC
            if n_classes == 2:
                auc_macro = roc_auc_score(y_encoded, y_prob[:, 1])
                auc_weighted = auc_macro
            else:
                auc_macro = roc_auc_score(y_encoded, y_prob, multi_class='ovr', average='macro')
                auc_weighted = roc_auc_score(y_encoded, y_prob, multi_class='ovr', average='weighted')
        except Exception:
            auc_macro = None
            auc_weighted = None

        # Calculate metrics
        metrics = {
            'accuracy': accuracy_score(y_encoded, y_pred),
            'balanced_accuracy': balanced_accuracy_score(y_encoded, y_pred),
            'f1_macro': f1_score(y_encoded, y_pred, average='macro'),
            'f1_weighted': f1_score(y_encoded, y_pred, average='weighted'),
            'auc_macro': auc_macro,
            'auc_weighted': auc_weighted,
            'n_samples': len(y),
            'n_classes': n_classes,
            'classes': le.classes_.tolist()
        }

        # Confusion matrix
        metrics['confusion_matrix'] = confusion_matrix(y_encoded, y_pred).tolist()

        return metrics

    def run_benchmark(self, tasks=None, classifiers=None, n_splits=5):
        """Run the full benchmark"""
        if tasks is None:
            tasks = list(self.tasks.keys())
        if classifiers is None:
            classifiers = list(self.classifiers.keys())

        self.n_splits = n_splits  # Store for use in evaluate_classifier

        logger.info("="*80)
        logger.info("COVID SEVERITY PREDICTION BENCHMARK")
        logger.info(f"CV Folds: {n_splits}")
        logger.info("="*80)

        results = {}

        for task_name in tasks:
            logger.info(f"\n{'='*60}")
            logger.info(f"Task: {task_name}")
            logger.info(f"{'='*60}")

            results[task_name] = {}

            for model_name, X in self.embeddings.items():
                logger.info(f"\n--- {model_name} ---")

                # Get labels for this model
                model_labels = self.model_labels[model_name]

                # Prepare labels for this task
                y, valid_mask = self.prepare_labels(task_name, model_labels)

                X_valid = X[valid_mask]
                y_valid = y[valid_mask]

                # Get patient/sample groups for GroupKFold
                model_groups = self.model_groups[model_name]
                groups_valid = model_groups[valid_mask]

                logger.info(f"Samples: {len(y_valid)}, Classes: {np.unique(y_valid)}")

                results[task_name][model_name] = {}

                for clf_name in classifiers:
                    clf = self.classifiers[clf_name]

                    try:
                        # Clone classifier for fresh evaluation
                        from sklearn.base import clone
                        clf_clone = clone(clf)

                        metrics = self.evaluate_classifier(X_valid, y_valid, clf_clone, n_splits=self.n_splits, groups=groups_valid)
                        results[task_name][model_name][clf_name] = metrics

                        auc_str = f"{metrics['auc_macro']:.3f}" if metrics['auc_macro'] is not None else 'N/A'
                        logger.info(f"  {clf_name}: Acc={metrics['accuracy']:.3f}, "
                                   f"BalAcc={metrics['balanced_accuracy']:.3f}, "
                                   f"AUC={auc_str}")

                    except Exception as e:
                        logger.error(f"  {clf_name}: Failed - {str(e)}")
                        results[task_name][model_name][clf_name] = {'error': str(e)}

        self.results = results
        return results

    def create_summary_table(self):
        """Create summary table of results"""
        rows = []

        for task_name, task_results in self.results.items():
            for model_name, model_results in task_results.items():
                for clf_name, metrics in model_results.items():
                    if 'error' not in metrics:
                        rows.append({
                            'Task': task_name,
                            'Embedding': model_name,
                            'Classifier': clf_name,
                            'Accuracy': metrics['accuracy'],
                            'Balanced_Accuracy': metrics['balanced_accuracy'],
                            'F1_Macro': metrics['f1_macro'],
                            'F1_Weighted': metrics['f1_weighted'],
                            'AUC_Macro': metrics['auc_macro'],
                            'AUC_Weighted': metrics['auc_weighted'],
                            'N_Samples': metrics['n_samples'],
                            'N_Classes': metrics['n_classes']
                        })

        df = pd.DataFrame(rows)
        return df

    def plot_results(self, save=True):
        """Generate visualization plots"""
        df = self.create_summary_table()

        if df.empty:
            logger.warning("No results to plot")
            return

        # 1. Comparison by embedding model (best classifier for each)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        for idx, task in enumerate(['binary', '3-class', '6-class']):
            task_df = df[df['Task'] == task]
            if task_df.empty:
                continue

            # Get best classifier for each embedding
            best_df = task_df.loc[task_df.groupby('Embedding')['Balanced_Accuracy'].idxmax()]

            ax = axes[idx]
            x = np.arange(len(best_df))
            width = 0.35

            bars1 = ax.bar(x - width/2, best_df['Balanced_Accuracy'], width, label='Balanced Acc', color='steelblue')
            bars2 = ax.bar(x + width/2, best_df['AUC_Macro'].fillna(0), width, label='AUC Macro', color='coral')

            ax.set_xlabel('Embedding Model')
            ax.set_ylabel('Score')
            ax.set_title(f'{task} Classification')
            ax.set_xticks(x)
            ax.set_xticklabels(best_df['Embedding'], rotation=45, ha='right')
            ax.legend()
            ax.set_ylim(0, 1)

            # Add value labels
            for bar in bars1:
                height = bar.get_height()
                ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width()/2, height),
                           xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8)

        plt.tight_layout()
        if save:
            plt.savefig(self.output_dir / 'embedding_comparison.png', dpi=150, bbox_inches='tight')
            plt.savefig(self.output_dir / 'embedding_comparison.pdf', bbox_inches='tight')
        plt.close()

        # 2. Heatmap of all results
        for task in df['Task'].unique():
            task_df = df[df['Task'] == task]

            # Pivot for heatmap
            pivot_df = task_df.pivot_table(
                index='Embedding', columns='Classifier',
                values='Balanced_Accuracy', aggfunc='mean'
            )

            fig, ax = plt.subplots(figsize=(10, 6))
            sns.heatmap(pivot_df, annot=True, fmt='.3f', cmap='RdYlGn',
                       vmin=0, vmax=1, ax=ax)
            ax.set_title(f'{task} Classification - Balanced Accuracy')
            plt.tight_layout()

            if save:
                plt.savefig(self.output_dir / f'heatmap_{task}.png', dpi=150, bbox_inches='tight')
            plt.close()

        # 3. Bar plot comparing all classifiers
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        for idx, task in enumerate(['binary', '3-class', '6-class']):
            task_df = df[df['Task'] == task]
            if task_df.empty:
                continue

            ax = axes[idx]

            # Group by embedding and classifier
            pivot = task_df.pivot(index='Classifier', columns='Embedding', values='Balanced_Accuracy')
            pivot.plot(kind='bar', ax=ax, rot=45)

            ax.set_title(f'{task} Classification')
            ax.set_ylabel('Balanced Accuracy')
            ax.set_ylim(0, 1)
            ax.legend(title='Embedding', bbox_to_anchor=(1.02, 1), loc='upper left')

        plt.tight_layout()
        if save:
            plt.savefig(self.output_dir / 'classifier_comparison.png', dpi=150, bbox_inches='tight')
        plt.close()

        logger.info(f"Plots saved to {self.output_dir}")

    def save_results(self):
        """Save all results to files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Save detailed JSON results
        json_file = self.output_dir / f'benchmark_results_{timestamp}.json'
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)
        logger.info(f"Saved JSON results: {json_file}")

        # Save summary table
        df = self.create_summary_table()
        csv_file = self.output_dir / f'benchmark_summary_{timestamp}.csv'
        df.to_csv(csv_file, index=False)
        logger.info(f"Saved CSV summary: {csv_file}")

        # Save best results summary
        summary_file = self.output_dir / f'best_results_{timestamp}.txt'
        with open(summary_file, 'w') as f:
            f.write("COVID SEVERITY PREDICTION BENCHMARK - BEST RESULTS\n")
            f.write("="*60 + "\n\n")

            for task in df['Task'].unique():
                task_df = df[df['Task'] == task]
                best_idx = task_df['Balanced_Accuracy'].idxmax()
                best = task_df.loc[best_idx]

                f.write(f"\n{task.upper()} CLASSIFICATION:\n")
                f.write(f"  Best Embedding: {best['Embedding']}\n")
                f.write(f"  Best Classifier: {best['Classifier']}\n")
                f.write(f"  Balanced Accuracy: {best['Balanced_Accuracy']:.4f}\n")
                f.write(f"  AUC Macro: {best['AUC_Macro']:.4f}\n")
                f.write(f"  F1 Macro: {best['F1_Macro']:.4f}\n")

            # Embedding ranking
            f.write("\n\nEMBEDDING RANKING (by average Balanced Accuracy):\n")
            f.write("-"*40 + "\n")
            ranking = df.groupby('Embedding')['Balanced_Accuracy'].mean().sort_values(ascending=False)
            for i, (emb, score) in enumerate(ranking.items(), 1):
                f.write(f"  {i}. {emb}: {score:.4f}\n")

        logger.info(f"Saved best results: {summary_file}")

        return csv_file


def main():
    """Main execution"""
    import argparse
    parser = argparse.ArgumentParser(description='COVID Severity Prediction Benchmark')
    parser.add_argument('--base-dir', type=str,
                        default="/bigdata/godziklab/shared/Xinru/302006",
                        help='Base directory for the project')
    parser.add_argument('--tasks', nargs='+', default=None,
                        choices=['6-class', '3-class', 'binary'],
                        help='Classification tasks to run')
    parser.add_argument('--classifiers', nargs='+', default=None,
                        help='Classifiers to evaluate')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: only LogisticRegression and RandomForest')
    parser.add_argument('--fastest', action='store_true',
                        help='Fastest mode: binary task only, LogisticRegression, 3-fold CV')
    parser.add_argument('--cv-folds', type=int, default=5,
                        help='Number of cross-validation folds (default: 5, use 3 for faster)')
    args = parser.parse_args()

    # Initialize benchmark
    benchmark = COVIDSeverityBenchmark(base_dir=args.base_dir)

    # Load embeddings
    if not benchmark.load_embeddings():
        logger.error("Failed to load embeddings")
        return 1

    # Set classifiers and tasks based on mode
    classifiers = args.classifiers
    tasks = args.tasks
    n_splits = args.cv_folds

    if args.fastest:
        # Fastest mode: binary only, LogisticRegression only, 3-fold CV
        classifiers = ['LogisticRegression']
        tasks = ['binary']
        n_splits = 3
        logger.info("FASTEST MODE: binary task, LogisticRegression, 3-fold CV")
    elif args.quick:
        classifiers = ['LogisticRegression', 'RandomForest']

    # Run benchmark
    benchmark.run_benchmark(tasks=tasks, classifiers=classifiers, n_splits=n_splits)

    # Generate outputs
    benchmark.save_results()
    benchmark.plot_results()

    # Print summary
    df = benchmark.create_summary_table()
    print("\n" + "="*80)
    print("SUMMARY: Best Results by Task")
    print("="*80)

    for task in df['Task'].unique():
        task_df = df[df['Task'] == task]
        best_idx = task_df['Balanced_Accuracy'].idxmax()
        best = task_df.loc[best_idx]
        print(f"\n{task}: {best['Embedding']} + {best['Classifier']}")
        print(f"  Balanced Acc: {best['Balanced_Accuracy']:.4f}, AUC: {best['AUC_Macro']:.4f}")

    print("\n" + "="*80)
    print("EMBEDDING RANKING (avg Balanced Accuracy across all tasks):")
    print("="*80)
    ranking = df.groupby('Embedding')['Balanced_Accuracy'].mean().sort_values(ascending=False)
    for i, (emb, score) in enumerate(ranking.items(), 1):
        print(f"  {i}. {emb}: {score:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
