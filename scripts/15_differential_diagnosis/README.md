# 15 - Differential Diagnosis of Inflammatory Conditions

**Accession:** 15
**Status:** ✓ DONE — results in `results/differential_diagnosis_5class_results.csv`
**Priority:** HIGH - Key differentiator for paper

**Headline numbers (v5, patient-level CV):** UCE Macro AUC 0.938 (best),
STATE 0.936, TF 0.871, scGPT 0.842, GF 0.795. UCE > STATE here is a ranking
reversal vs severity prediction.

**Note:** to be re-run on Option C canonical-input embeddings (TODO §9.5)
once those land.

---

## Purpose

Assess whether foundation models can distinguish COVID-19 from other inflammatory conditions that present with similar symptoms. This is the **clinically critical task** that elevates this benchmark from academic exercise to clinical utility demonstration.

---

## Clinical Motivation

A physician seeing a patient with fever and elevated inflammatory markers needs to distinguish between:
- COVID-19 (requires isolation, specific antivirals)
- Sepsis (requires broad-spectrum antibiotics, source control)
- SLE flare (requires immunosuppression)
- Flu (supportive care, antivirals if early)
- Bacterial pneumonia (antibiotics)

**A model that can make this distinction from platelet signatures is clinically transformative.**

---

## Task Definition

### 6-Class Classification Problem

| Class | Disease | Expected Samples | Clinical Action |
|-------|---------|------------------|-----------------|
| 0 | Healthy | ~X,XXX | Reassurance |
| 1 | COVID-19 | ~43,000 (90%) | Isolation, antivirals |
| 2 | Sepsis | ~X,XXX | Antibiotics, ICU |
| 3 | SLE | ~X,XXX | Immunosuppression |
| 4 | Flu | ~X,XXX | Antivirals if early |
| 5 | Pneumonia | ~X,XXX | Antibiotics |

### Dataset
- **Full platelet dataset** (all diseases included)
- Stratified sampling to handle class imbalance
- 5-fold cross-validation

---

## Planned Analysis

### 1. Multi-class Classification
```python
# 6-class: Healthy, COVID, Sepsis, SLE, Flu, Pneumonia
y = adata.obs['disease_condition']
# One-vs-Rest AUC for each class
```

### 2. Confusion Matrix Analysis
- Which diseases are most confused?
- Do models confuse COVID with Flu? (both respiratory)
- Do models confuse Sepsis with severe COVID? (both critical)

### 3. Disease-Specific Performance
| Model | COVID AUC | Sepsis AUC | SLE AUC | Flu AUC | Pneumonia AUC |
|-------|-----------|------------|---------|---------|---------------|
| STATE | ? | ? | ? | ? | ? |
| UCE | ? | ? | ? | ? | ? |
| ... | | | | | |

### 4. Hierarchical Analysis
```
Level 1: Healthy vs Any Disease
Level 2: Infectious vs Autoimmune (COVID/Sepsis/Flu/Pneumonia vs SLE)
Level 3: Viral vs Bacterial (COVID/Flu vs Sepsis/Pneumonia)
Level 4: Specific disease identification
```

---

## Expected Output

### Results
| File | Description |
|------|-------------|
| `results/differential_diagnosis_6class.csv` | Per-model, per-class AUC |
| `results/confusion_matrix_STATE.csv` | 6x6 confusion matrix |
| `results/confusion_matrix_UCE.csv` | 6x6 confusion matrix |
| `results/disease_specific_metrics.csv` | Precision/Recall per disease |
| `results/hierarchical_classification.csv` | Results at each hierarchy level |

### Figures
| File | Description |
|------|-------------|
| `figures/fig_15_confusion_matrix_heatmap.png` | 6x6 confusion matrices for all models |
| `figures/fig_15_disease_auc_comparison.png` | Per-disease AUC bar chart |
| `figures/fig_15_umap_6disease.png` | UMAP colored by disease |
| `figures/fig_15_hierarchical_performance.png` | Performance at each hierarchy level |

---

## Key Questions to Answer

1. **Can models distinguish COVID from other respiratory infections (Flu, Pneumonia)?**
   - This tests specificity beyond "inflammation detection"

2. **Which disease pairs are most confusable?**
   - Clinically important: COVID vs Flu confusion is common

3. **Do models learn disease-specific signatures or generic inflammation?**
   - If all diseases cluster together → models see generic inflammation
   - If diseases separate → models learn specific signatures

4. **Which model is best for differential diagnosis?**
   - May differ from COVID severity ranking!

---

## Hypothesis

Based on prior results:
- **STATE**: Likely best (best OOD generalization suggests disease-specific learning)
- **UCE**: Good robustness may help with smaller disease classes
- **scGPT**: Likely fails (failed on COVID→Sepsis transfer)
- **Geneformer**: Distributed information may help multi-class

---

## Code Template

```python
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.preprocessing import LabelEncoder

def evaluate_differential_diagnosis(X, y_disease, n_splits=5):
    """
    Evaluate 6-class differential diagnosis performance.

    Parameters
    ----------
    X : np.ndarray
        Embedding matrix (n_cells, n_dims)
    y_disease : np.ndarray
        Disease labels (Healthy, COVID, Sepsis, SLE, Flu, Pneumonia)

    Returns
    -------
    dict with AUC per class, confusion matrix, etc.
    """
    le = LabelEncoder()
    y = le.fit_transform(y_disease)
    n_classes = len(le.classes_)

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    all_y_true = []
    all_y_pred_proba = []

    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        clf = LogisticRegression(
            solver='lbfgs',
            max_iter=1000,
            class_weight='balanced',
            multi_class='multinomial'
        )
        clf.fit(X_train, y_train)

        y_pred_proba = clf.predict_proba(X_test)

        all_y_true.extend(y_test)
        all_y_pred_proba.append(y_pred_proba)

    all_y_true = np.array(all_y_true)
    all_y_pred_proba = np.vstack(all_y_pred_proba)

    # Per-class AUC (One-vs-Rest)
    class_aucs = {}
    for i, class_name in enumerate(le.classes_):
        y_binary = (all_y_true == i).astype(int)
        class_aucs[class_name] = roc_auc_score(y_binary, all_y_pred_proba[:, i])

    # Macro-averaged AUC
    macro_auc = roc_auc_score(all_y_true, all_y_pred_proba, multi_class='ovr')

    # Confusion matrix
    y_pred = np.argmax(all_y_pred_proba, axis=1)
    cm = confusion_matrix(all_y_true, y_pred)

    return {
        'class_aucs': class_aucs,
        'macro_auc': macro_auc,
        'confusion_matrix': cm,
        'class_names': le.classes_
    }

# Run for each model
models = ['STATE', 'UCE', 'Geneformer', 'scGPT', 'TranscriptFormer']
results = {}

for model in models:
    X = load_embeddings(model)
    y_disease = adata.obs['disease']
    results[model] = evaluate_differential_diagnosis(X, y_disease)

    print(f"\n{model}:")
    print(f"  Macro AUC: {results[model]['macro_auc']:.3f}")
    for disease, auc in results[model]['class_aucs'].items():
        print(f"  {disease}: {auc:.3f}")
```

---

## Manuscript Integration

This becomes **Section 3.2** of Results:

> **3.2 Differential Diagnosis of Inflammatory Conditions**
>
> Beyond severity prediction within COVID-19, we assessed whether foundation models can distinguish COVID-19 from other inflammatory conditions—a clinically critical capability. We framed this as a 6-class classification problem: Healthy, COVID-19, Sepsis, SLE, Influenza, and Pneumonia.
>
> [Results Table]
>
> STATE achieved the highest macro-averaged AUC (X.XX), demonstrating superior ability to learn disease-specific signatures rather than generic inflammatory patterns. Notably, the confusion matrix revealed that [key finding about which diseases are confused].

---

## Dependencies

- `scikit-learn>=1.3`
- `matplotlib`
- `seaborn`
- Pre-computed embeddings from `01_embedding_integrity`
- Full dataset with disease labels
