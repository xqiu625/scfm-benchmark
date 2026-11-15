#!/usr/bin/env python3
"""
Cross-Disease Generalization Analysis

Train on COVID, test on Sepsis (and vice versa) to test true out-of-distribution generalization.
This tests whether foundation model embeddings capture generalizable disease patterns.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, classification_report
from pathlib import Path
from datetime import datetime
import logging
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "cross_disease_generalization"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Models
MODELS = {
    'STATE': EMB_DIR / "state/covid_state.h5ad",
    'UCE': EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
    'TranscriptFormer': EMB_DIR / "transcriptformer/covid_transcriptformer.h5ad",
    'scGPT': EMB_DIR / "scgpt/covid_scgpt.h5ad",
    'Geneformer': EMB_DIR / "geneformer/covid_geneformer.h5ad",
}

OBSM_KEYS = {
    'STATE': ['X_emb', 'X_state', 'X_STATE'],
    'UCE': ['X_uce', 'X_uce_4layer'],
    'scGPT': ['X_scGPT', 'X_scgpt'],
    'TranscriptFormer': ['X_transcriptformer', 'X_TranscriptFormer'],
    'Geneformer': ['X_geneformer', 'X_Geneformer'],
}


def load_embeddings(model_name):
    """Load embeddings and metadata."""
    path = MODELS[model_name]
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    logger.info(f"Loading {model_name} from {path}")
    adata = sc.read_h5ad(path)

    # Get embeddings — try all known keys for this model
    X = None
    keys_to_try = OBSM_KEYS.get(model_name, [])
    for key in keys_to_try:
        if key in adata.obsm:
            X = adata.obsm[key]
            logger.info(f"  Using embedding key: {key}")
            break

    if X is None:
        raise ValueError(
            f"No embedding key found for {model_name} in {path}. "
            f"Tried: {keys_to_try}. Available obsm keys: {list(adata.obsm.keys())}. "
            f"Refusing to fall back to raw expression (adata.X)."
        )

    if hasattr(X, 'toarray'):
        X = X.toarray()

    return np.array(X), adata

def get_disease_and_severity(adata):
    """Extract disease type and severity labels."""
    obs = adata.obs

    # Get disease type
    if 'Disease' in obs.columns:
        disease = obs['Disease'].values
    elif 'disease' in obs.columns:
        disease = obs['disease'].values
    else:
        # Try to infer from Data_NO or other columns
        disease = np.array(['Unknown'] * len(obs))

    # Get severity
    severity = None
    for col in ['covid_severity', 'severity', 'Severity', 'covid_severity_label']:
        if col in obs.columns:
            severity = obs[col].astype(str).str.lower().values
            logger.info(f"  Using severity column: {col}")
            break
    if severity is None:
        raise ValueError(
            f"No severity column found. Available: {list(obs.columns)[:20]}"
        )

    return disease, severity

def create_binary_labels(severity):
    """Create binary severe/non-severe labels."""
    severe_classes = ['severe', 'fatal']
    binary = np.array([1 if s in severe_classes else 0 for s in severity])
    return binary

def split_by_disease(X, disease, severity):
    """Split data by disease type."""
    # Find COVID and Sepsis samples
    covid_mask = np.isin(disease, ['COVID', 'COVID-19', 'covid', 'COVID19'])
    sepsis_mask = np.isin(disease, ['Sepsis', 'sepsis', 'SEPSIS'])

    if covid_mask.sum() == 0 and sepsis_mask.sum() == 0:
        # No disease labels found - cannot perform cross-disease analysis
        raise ValueError(
            "No disease labels found in data. Cross-disease generalization analysis "
            "requires explicit COVID and Sepsis labels. Please ensure the data contains "
            "'Disease' column with values like 'COVID', 'COVID-19', 'Sepsis', etc."
        )

    logger.info(f"  COVID samples: {covid_mask.sum()}")
    logger.info(f"  Sepsis samples: {sepsis_mask.sum()}")

    return covid_mask, sepsis_mask

def train_test_cross_disease(X, severity, disease, train_disease='COVID-19', test_disease='Sepsis'):
    """Train on one disease, test on another."""

    # Create masks — match known aliases
    DISEASE_ALIASES = {
        'COVID-19': ['COVID-19', 'covid-19', 'COVID', 'covid', 'COVID19'],
        'Sepsis': ['Sepsis', 'sepsis', 'SEPSIS'],
    }
    train_variants = DISEASE_ALIASES.get(train_disease, [train_disease])
    test_variants = DISEASE_ALIASES.get(test_disease, [test_disease])
    train_mask = np.isin(disease, train_variants)
    test_mask = np.isin(disease, test_variants)

    if train_mask.sum() == 0 or test_mask.sum() == 0:
        return None

    X_train = X[train_mask]
    X_test = X[test_mask]
    y_train = create_binary_labels(severity[train_mask])
    y_test = create_binary_labels(severity[test_mask])

    # Check class balance
    logger.info(f"  Train ({train_disease}): {len(y_train)} samples, {y_train.sum()} severe")
    logger.info(f"  Test ({test_disease}): {len(y_test)} samples, {y_test.sum()} severe")

    if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        logger.warning("  Not enough class diversity, skipping...")
        return None

    # Scale
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    results = {}

    # Logistic Regression
    clf_lr = LogisticRegression(max_iter=1000, random_state=42)
    clf_lr.fit(X_train_scaled, y_train)
    y_prob_lr = clf_lr.predict_proba(X_test_scaled)[:, 1]
    y_pred_lr = clf_lr.predict(X_test_scaled)

    results['LogReg'] = {
        'auc': roc_auc_score(y_test, y_prob_lr),
        'balanced_accuracy': balanced_accuracy_score(y_test, y_pred_lr)
    }

    # Random Forest
    clf_rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
    clf_rf.fit(X_train_scaled, y_train)
    y_prob_rf = clf_rf.predict_proba(X_test_scaled)[:, 1]
    y_pred_rf = clf_rf.predict(X_test_scaled)

    results['RandomForest'] = {
        'auc': roc_auc_score(y_test, y_prob_rf),
        'balanced_accuracy': balanced_accuracy_score(y_test, y_pred_rf)
    }

    return results

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 80)
    logger.info("CROSS-DISEASE GENERALIZATION ANALYSIS")
    logger.info("Train on COVID → Test on Sepsis (and vice versa)")
    logger.info("=" * 80)

    all_results = []

    for model_name, path in MODELS.items():
        if not path.exists():
            logger.warning(f"Skipping {model_name}: file not found")
            continue

        logger.info(f"\n{'='*60}")
        logger.info(f"Model: {model_name}")
        logger.info(f"{'='*60}")

        try:
            X, adata = load_embeddings(model_name)
            disease, severity = get_disease_and_severity(adata)

            unique_diseases = np.unique(disease)
            logger.info(f"  Unique diseases: {unique_diseases}")

            # COVID → Sepsis
            logger.info("\n  Direction: COVID → Sepsis")
            results_covid_sepsis = train_test_cross_disease(X, severity, disease, 'COVID-19', 'Sepsis')

            if results_covid_sepsis:
                for clf_name, metrics in results_covid_sepsis.items():
                    all_results.append({
                        'model': model_name,
                        'train': 'COVID',
                        'test': 'Sepsis',
                        'classifier': clf_name,
                        'auc': metrics['auc'],
                        'balanced_accuracy': metrics['balanced_accuracy']
                    })
                    logger.info(f"    {clf_name}: AUC={metrics['auc']:.4f}, BalAcc={metrics['balanced_accuracy']:.4f}")

            # Sepsis → COVID
            logger.info("\n  Direction: Sepsis → COVID")
            results_sepsis_covid = train_test_cross_disease(X, severity, disease, 'Sepsis', 'COVID-19')

            if results_sepsis_covid:
                for clf_name, metrics in results_sepsis_covid.items():
                    all_results.append({
                        'model': model_name,
                        'train': 'Sepsis',
                        'test': 'COVID',
                        'classifier': clf_name,
                        'auc': metrics['auc'],
                        'balanced_accuracy': metrics['balanced_accuracy']
                    })
                    logger.info(f"    {clf_name}: AUC={metrics['auc']:.4f}, BalAcc={metrics['balanced_accuracy']:.4f}")

        except Exception as e:
            logger.error(f"Error processing {model_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Save results
    if all_results:
        results_df = pd.DataFrame(all_results)
        output_file = RESULTS_DIR / f"cross_disease_generalization_{timestamp}.csv"
        results_df.to_csv(output_file, index=False)
        logger.info(f"\nResults saved to {output_file}")

        # Summary
        logger.info("\n" + "=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)

        for direction in [('COVID', 'Sepsis'), ('Sepsis', 'COVID')]:
            logger.info(f"\n{direction[0]} → {direction[1]}:")
            subset = results_df[(results_df['train'] == direction[0]) & (results_df['test'] == direction[1])]
            if len(subset) > 0:
                for model in subset['model'].unique():
                    model_data = subset[subset['model'] == model]
                    best = model_data.loc[model_data['auc'].idxmax()]
                    logger.info(f"  {model}: AUC={best['auc']:.4f} ({best['classifier']})")
    else:
        logger.warning("No results generated - check if disease labels exist in data")

if __name__ == "__main__":
    main()
