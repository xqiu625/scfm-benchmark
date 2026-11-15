#!/usr/bin/env python3
"""
Shared cross-disease transfer utilities -- code review Fix-2/Fix-3
(2026-09-04, reports/manuscript/revision/code_review_2026-09-04.md), used by
both cross_disease_fatal_vs_survived_patient_v2.py (Job H2, legacy 12-patient
sepsis arm) and jobS_cross_disease.py (Chongqing 91-patient sepsis arm).

Fix-2: the AUC point estimate was cell-weighted (bootstrap resampled patients,
but the score itself was computed on raw cells) -- inconsistent with P0's
estimand (everything scored patient-level). PRIMARY here is patient-level:
each patient's cell predictions are averaged first, then AUC + the 200-iter
bootstrap are both computed on that small patient-level array. Cell-weighted
AUC is kept as a secondary column (`auc_cell_weighted`), computed once, no
bootstrap (matches how it's used elsewhere in this revision: a representation-
level sanity reference, not a headline).

Fix-3: P0 has decided severe=survived (severe-INCLUDED) is the primary COVID
protocol; severe-excluded is a sensitivity variant. Both are computed and
reported side by side (column prefixes primary_ / sens_severe_excluded_),
no adjudication in code -- that's a manuscript decision, not a script's job.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

N_BOOTSTRAP = 200

COVID_FATAL = {'fatal'}
COVID_SURVIVED_PRIMARY = {'mild', 'moderate', 'recovered', 'severe'}  # Fix-3: severe=survived
COVID_SURVIVED_SENS = {'mild', 'moderate', 'recovered'}               # severe-excluded (sensitivity)
SEPSIS_FATAL = {'fatal'}
SEPSIS_SURVIVED = {'severe'}


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return float('nan')
    return float(roc_auc_score(y, p))


def aggregate_patient(y_cell, p_cell, patient_ids):
    """Fix-2: mean cell-level score per patient -> one row per patient."""
    df = pd.DataFrame({'patient_id': patient_ids, 'y': y_cell, 'p': p_cell})
    g = df.groupby('patient_id').agg(y=('y', 'max'), p=('p', 'mean'))
    return g['y'].to_numpy(), g['p'].to_numpy()


def patient_bootstrap_ci(y_pat, p_pat, n_iter=N_BOOTSTRAP, seed=42):
    """Bootstrap resampling PATIENTS on the already patient-aggregated scores."""
    rng = np.random.default_rng(seed)
    n = len(y_pat)
    aucs = []
    for _ in range(n_iter):
        idx = rng.integers(0, n, n)
        a = safe_auc(y_pat[idx], p_pat[idx])
        if not np.isnan(a):
            aucs.append(a)
    if not aucs:
        return float('nan'), float('nan'), float('nan')
    arr = np.asarray(aucs)
    return float(arr.mean()), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def transfer_one(Xtr, ytr, Xte, yte, gte):
    """Fit LR+RF on one train/test split. Returns dict per classifier:
    primary AUC = patient-level (Fix-2), plus boot CI on that, plus a
    secondary `auc_cell_weighted` column."""
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        return {}
    scaler = StandardScaler()
    Xtr_s = scaler.fit_transform(Xtr)
    Xte_s = scaler.transform(Xte)

    out = {'n_train': int(len(ytr)), 'n_test': int(len(yte)),
           'fatal_train': int(ytr.sum()), 'fatal_test': int(yte.sum()),
           'n_test_patients': int(len(np.unique(gte))),
           'n_fatal_test_patients': int(pd.DataFrame({'p': gte, 'y': yte})
                                        .groupby('p')['y'].max().sum())}

    for clf_name, clf in [
        ('LogReg', LogisticRegression(max_iter=2000, random_state=42, n_jobs=-1,
                                       class_weight='balanced')),
        ('RandomForest', RandomForestClassifier(n_estimators=200, random_state=42,
                                                 n_jobs=-1, class_weight='balanced')),
    ]:
        clf.fit(Xtr_s, ytr)
        p_cell = clf.predict_proba(Xte_s)[:, 1]
        auc_cell = safe_auc(yte, p_cell)
        bal_acc = float(balanced_accuracy_score(yte, clf.predict(Xte_s)))

        y_pat, p_pat = aggregate_patient(yte, p_cell, gte)
        auc_patient = safe_auc(y_pat, p_pat)
        boot_mean, boot_lo, boot_hi = patient_bootstrap_ci(y_pat, p_pat)

        out[clf_name] = {
            'auc': auc_patient,  # PRIMARY (Fix-2: patient-level)
            'boot_auc_mean': boot_mean, 'boot_ci_lo': boot_lo, 'boot_ci_hi': boot_hi,
            'balanced_accuracy': bal_acc,
            'auc_cell_weighted': auc_cell,  # secondary
        }
    return out
