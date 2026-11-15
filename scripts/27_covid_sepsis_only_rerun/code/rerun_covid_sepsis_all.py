#!/usr/bin/env python3
"""
Comprehensive COVID+Sepsis-only re-run: regenerates all cohort-affected tables.

Outputs (all under 05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/):
  - mortality_lr_rf.csv          -- Table 1 (already produced separately; recomputed for OOF)
  - oof_predictions/<FM>.csv     -- per-FM out-of-fold predictions (for downstream Table 3 + 4)
  - calibration.csv              -- Table 4 (Brier + ECE; from OOF predictions)
  - delong_pairwise.csv          -- Table 3 (pairwise DeLong)
  - pca_baseline.csv             -- Table 2 (PCA-50/100/200/500 + LR on raw counts)
  - lobo.csv                     -- Table 6 (leave-one-cohort-out)
  - per_cohort_oof.csv           -- Table 8 (per-cohort breakdown of OOF predictions)

Filter: Disease in {COVID-19, Sepsis} AND covid_severity in {fatal, mild, moderate, recovered, severe}
After filter: ~42,095 cells / 322 donors / 8,850 fatal / 33,245 survived (per-FM cell counts vary).

Memory-safe: h5py-only access to FM obsm; sparse X loaded once for PCA via h5py CSR reader.
"""
import h5py
import numpy as np
import pandas as pd
import os
import sys
import argparse

from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import TruncatedSVD
from sklearn.model_selection import StratifiedGroupKFold, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, brier_score_loss

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark/covid_sepsis_only'
PART_DIR = os.path.join(OUT_DIR, '_partials')

# When set (array per-FM mode), the per-FM stages (1 + 4) process only this model
# and write per-FM partials to PART_DIR instead of the aggregated CSVs. The gather
# step concatenates the partials back in EMBEDDINGS order — byte-identical to the
# monolithic run because all seeds are fixed (RANDOM_STATE) and CSV rows are
# formatted per-row independently of the others.
_ONLY_FM = None


def _selected_embeddings():
    if _ONLY_FM is None:
        return EMBEDDINGS
    return [e for e in EMBEDDINGS if e[0] == _ONLY_FM]


EMBEDDINGS = [
    ('STATE',            'state_covid_canonical.h5ad',            'X_state'),
    ('UCE',              'uce_covid_canonical.h5ad',              'X_uce'),
    ('TranscriptFormer', 'transcriptformer_covid_canonical.h5ad', 'X_transcriptformer'),
    ('scPRINT-2',        'scprint2_covid_canonical.h5ad',         'X_scprint2'),
    ('STACK',            'stack_covid_canonical.h5ad',            'X_stack'),
    ('Geneformer',       'geneformer_covid_canonical.h5ad',       'X_geneformer'),
    ('scGPT',            'scgpt_covid_canonical.h5ad',             'X_scgpt'),
]

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']
RANDOM_STATE = 42

# ============================================================================
# Helpers
# ============================================================================

def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def read_obs(h, col):
    if col in h['obs']:
        return read_categorical(h[f'obs/{col}'])
    return None


def load_filtered(path, obsm_key):
    """Load filtered FM embedding + obs columns (no X)."""
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_obs(h, 'Disease')
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
        cell_id = read_obs(h, 'Cell_ID')
        cohort = read_obs(h, 'Data_NO')

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    return {
        'X': X[keep],
        'severity': severity[keep],
        'sample_id': sample_id[keep],
        'cell_id': (cell_id[keep] if cell_id is not None else None),
        'cohort': (cohort[keep] if cohort is not None else None),
        'disease': disease[keep],
    }


def load_sparse_X(path):
    """Load X (sparse CSR) from an h5ad's X group, applying same filter."""
    with h5py.File(path, 'r') as h:
        # X is stored as a sparse CSR group: data, indices, indptr
        if 'X' in h and isinstance(h['X'], h5py.Group):
            data = h['X/data'][:]
            indices = h['X/indices'][:]
            indptr = h['X/indptr'][:]
            shape = tuple(h['X'].attrs.get('shape', [len(indptr) - 1, indices.max() + 1]))
            X = sparse.csr_matrix((data, indices, indptr), shape=shape)
        else:
            X = sparse.csr_matrix(h['X'][:])
        disease = read_obs(h, 'Disease')
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    return X[keep], severity[keep], sample_id[keep]


def make_skf():
    return StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)


def lr_classifier():
    return LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)


def rf_classifier():
    return RandomForestClassifier(n_estimators=100, class_weight='balanced',
                                  n_jobs=-1, random_state=RANDOM_STATE)


def expected_calibration_error(y_true, y_prob, n_bins=10):
    """ECE with equal-width bins."""
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    n = len(y_true)
    for i in range(n_bins):
        m = (y_prob >= bins[i]) & (y_prob < bins[i + 1])
        if i == n_bins - 1:
            m = (y_prob >= bins[i]) & (y_prob <= bins[i + 1])
        if m.sum() == 0:
            continue
        bin_acc = y_true[m].mean()
        bin_conf = y_prob[m].mean()
        ece += (m.sum() / n) * abs(bin_conf - bin_acc)
    return ece


def delong_var(y_true, y_pred):
    """DeLong covariance computation. Returns AUC + variance for one model."""
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    n_pos, n_neg = len(pos_idx), len(neg_idx)
    p = y_pred[pos_idx]
    n = y_pred[neg_idx]
    # Compute structural components V10 (pos) and V01 (neg)
    # V10[i] = (1/n_neg) * sum_j ind(p[i] > n[j]) + 0.5*ind(p[i]==n[j])
    V10 = np.array([
        (np.sum(p[i] > n) + 0.5 * np.sum(p[i] == n)) / n_neg
        for i in range(n_pos)
    ])
    V01 = np.array([
        (np.sum(p > n[j]) + 0.5 * np.sum(p == n[j])) / n_pos
        for j in range(n_neg)
    ])
    auc = V10.mean()
    var = V10.var(ddof=1) / n_pos + V01.var(ddof=1) / n_neg
    return auc, var, V10, V01


def delong_pvalue(y_true, y_a, y_b):
    """Two-sample DeLong test for paired ROC AUCs."""
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    n_pos, n_neg = len(pos_idx), len(neg_idx)
    pa, na = y_a[pos_idx], y_a[neg_idx]
    pb, nb = y_b[pos_idx], y_b[neg_idx]

    V10a = np.array([(np.sum(pa[i] > na) + 0.5 * np.sum(pa[i] == na)) / n_neg for i in range(n_pos)])
    V01a = np.array([(np.sum(pa > na[j]) + 0.5 * np.sum(pa == na[j])) / n_pos for j in range(n_neg)])
    V10b = np.array([(np.sum(pb[i] > nb) + 0.5 * np.sum(pb[i] == nb)) / n_neg for i in range(n_pos)])
    V01b = np.array([(np.sum(pb > nb[j]) + 0.5 * np.sum(pb == nb[j])) / n_pos for j in range(n_neg)])

    auc_a = V10a.mean()
    auc_b = V10b.mean()

    cov_pos = np.cov(np.vstack([V10a, V10b]), ddof=1)
    cov_neg = np.cov(np.vstack([V01a, V01b]), ddof=1)
    var_diff = (cov_pos[0, 0] + cov_pos[1, 1] - 2 * cov_pos[0, 1]) / n_pos + \
               (cov_neg[0, 0] + cov_neg[1, 1] - 2 * cov_neg[0, 1]) / n_neg
    if var_diff <= 0:
        return auc_a, auc_b, np.nan
    z = (auc_a - auc_b) / np.sqrt(var_diff)
    from scipy.stats import norm
    p = 2 * (1 - norm.cdf(abs(z)))
    return auc_a, auc_b, p


# ============================================================================
# Stage 1: Primary benchmark + OOF predictions + Calibration
# ============================================================================

def stage1_primary_calibration():
    print('\n' + '=' * 70, flush=True)
    print('STAGE 1: Primary mortality + OOF predictions + Calibration', flush=True)
    print('=' * 70, flush=True)

    os.makedirs(os.path.join(OUT_DIR, 'oof_predictions'), exist_ok=True)
    primary_rows = []
    cal_rows = []

    for fm, fname, key in _selected_embeddings():
        path = os.path.join(ROOT, fname)
        if not os.path.exists(path):
            print(f'SKIP: {path}', flush=True); continue
        print(f'\n--- {fm} ---', flush=True)
        D = load_filtered(path, key)
        X, sev, sid = D['X'], D['severity'], D['sample_id']
        cid = D['cell_id']
        y = (sev == 'fatal').astype(int)
        n_donors = len(np.unique(sid))
        print(f'  {X.shape[0]} cells, {n_donors} donors, fatal={y.sum()}', flush=True)

        skf = make_skf()
        oof_pred_lr = np.full(len(y), np.nan)
        oof_pred_rf = np.full(len(y), np.nan)
        fold_aucs_lr, fold_aucs_rf = [], []

        for fold, (tr, te) in enumerate(skf.split(X, y, groups=sid)):
            # LR
            clf_lr = lr_classifier()
            clf_lr.fit(X[tr], y[tr])
            p_lr = clf_lr.predict_proba(X[te])[:, 1]
            oof_pred_lr[te] = p_lr
            fold_aucs_lr.append(roc_auc_score(y[te], p_lr))
            # RF
            clf_rf = rf_classifier()
            clf_rf.fit(X[tr], y[tr])
            p_rf = clf_rf.predict_proba(X[te])[:, 1]
            oof_pred_rf[te] = p_rf
            fold_aucs_rf.append(roc_auc_score(y[te], p_rf))

        # Save OOF predictions
        oof_df = pd.DataFrame({
            'cell_id': cid if cid is not None else np.arange(len(y)),
            'sample_id': sid,
            'severity': sev,
            'y': y,
            'p_lr': oof_pred_lr,
            'p_rf': oof_pred_rf,
        })
        oof_df.to_csv(os.path.join(OUT_DIR, 'oof_predictions', f'{fm}.csv'), index=False)

        # Primary metrics
        primary_rows.append({
            'fm': fm, 'classifier': 'LR',
            'n_cells': int(X.shape[0]), 'n_donors': n_donors,
            'auc_mean': float(np.mean(fold_aucs_lr)), 'auc_std': float(np.std(fold_aucs_lr)),
            'fold_aucs': ';'.join(f'{a:.3f}' for a in fold_aucs_lr),
        })
        primary_rows.append({
            'fm': fm, 'classifier': 'RF',
            'n_cells': int(X.shape[0]), 'n_donors': n_donors,
            'auc_mean': float(np.mean(fold_aucs_rf)), 'auc_std': float(np.std(fold_aucs_rf)),
            'fold_aucs': ';'.join(f'{a:.3f}' for a in fold_aucs_rf),
        })

        # Calibration on OOF predictions
        for clf_name, p_oof in [('LR', oof_pred_lr), ('RF', oof_pred_rf)]:
            mask_valid = ~np.isnan(p_oof)
            yv = y[mask_valid]; pv = p_oof[mask_valid]
            auc = roc_auc_score(yv, pv)
            brier = brier_score_loss(yv, pv)
            ece = expected_calibration_error(yv, pv, n_bins=10)
            cal_rows.append({
                'fm': fm, 'classifier': clf_name,
                'n_cells_valid': int(mask_valid.sum()),
                'auc': auc, 'brier': brier, 'ece': ece,
            })
            print(f'  {clf_name}: AUC={auc:.3f}  Brier={brier:.3f}  ECE={ece:.3f}', flush=True)

    if _ONLY_FM is not None:
        os.makedirs(PART_DIR, exist_ok=True)
        pd.DataFrame(primary_rows).to_csv(os.path.join(PART_DIR, f'mortality_{_ONLY_FM}.csv'), index=False)
        pd.DataFrame(cal_rows).to_csv(os.path.join(PART_DIR, f'calibration_{_ONLY_FM}.csv'), index=False)
        print(f'\nSaved partials for {_ONLY_FM}: mortality, calibration, oof_predictions/{_ONLY_FM}.csv', flush=True)
    else:
        pd.DataFrame(primary_rows).to_csv(os.path.join(OUT_DIR, 'mortality_lr_rf.csv'), index=False)
        pd.DataFrame(cal_rows).to_csv(os.path.join(OUT_DIR, 'calibration.csv'), index=False)
        print(f'\nSaved: mortality_lr_rf.csv, calibration.csv, oof_predictions/<FM>.csv', flush=True)


# ============================================================================
# Stage 2: DeLong pairwise on aligned cell set
# ============================================================================

def stage2_delong():
    print('\n' + '=' * 70, flush=True)
    print('STAGE 2: DeLong pairwise tests on Cell_ID-aligned subset', flush=True)
    print('=' * 70, flush=True)

    oof_dfs = {}
    for fm, _, _ in EMBEDDINGS:
        p = os.path.join(OUT_DIR, 'oof_predictions', f'{fm}.csv')
        if not os.path.exists(p):
            continue
        oof_dfs[fm] = pd.read_csv(p)

    # Find Cell_ID intersection across all FMs
    common_cids = None
    for fm, df in oof_dfs.items():
        cids = set(df['cell_id'].astype(str))
        common_cids = cids if common_cids is None else (common_cids & cids)
    print(f'  Cell_ID-aligned subset: {len(common_cids)} cells', flush=True)

    # Reduce each to the aligned set, sort by cell_id for consistent ordering
    aligned = {}
    y_aligned = None
    for fm, df in oof_dfs.items():
        df_a = df[df['cell_id'].astype(str).isin(common_cids)].sort_values('cell_id').reset_index(drop=True)
        aligned[fm] = df_a['p_lr'].values
        if y_aligned is None:
            y_aligned = df_a['y'].values

    fms = list(aligned.keys())
    rows = []
    for i in range(len(fms)):
        for j in range(i + 1, len(fms)):
            a_fm, b_fm = fms[i], fms[j]
            auc_a, auc_b, p = delong_pvalue(y_aligned, aligned[a_fm], aligned[b_fm])
            rows.append({
                'model_a': a_fm, 'model_b': b_fm,
                'auc_a': auc_a, 'auc_b': auc_b,
                'delta_auc': auc_a - auc_b, 'p_value': p,
            })
            print(f'  {a_fm} vs {b_fm}: ΔAUC={auc_a-auc_b:+.4f}, p={p:.2e}', flush=True)

    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'delong_pairwise.csv'), index=False)
    print(f'\nSaved: delong_pairwise.csv', flush=True)


# ============================================================================
# Stage 3: PCA baseline (raw counts via X from STATE h5ad)
# ============================================================================

def stage3_pca_baseline():
    print('\n' + '=' * 70, flush=True)
    print('STAGE 3: PCA baseline on raw counts', flush=True)
    print('=' * 70, flush=True)

    # Use STATE's h5ad as the source of raw X (canonical 47977-cell input)
    state_path = os.path.join(ROOT, 'state_covid_canonical.h5ad')
    print(f'  Loading sparse X from {os.path.basename(state_path)}', flush=True)
    X, sev, sid = load_sparse_X(state_path)
    y = (sev == 'fatal').astype(int)
    print(f'  Filtered X: {X.shape}, fatal={y.sum()}, donors={len(np.unique(sid))}', flush=True)

    skf = make_skf()
    rows = []

    for n_comp in [50, 100, 200, 500]:
        print(f'  -- PCA-{n_comp} + LR --', flush=True)
        fold_aucs = []
        for fold, (tr, te) in enumerate(skf.split(X, y, groups=sid)):
            svd = TruncatedSVD(n_components=n_comp, random_state=RANDOM_STATE)
            Z_tr = svd.fit_transform(X[tr])
            Z_te = svd.transform(X[te])
            clf = lr_classifier()
            clf.fit(Z_tr, y[tr])
            p = clf.predict_proba(Z_te)[:, 1]
            auc = roc_auc_score(y[te], p)
            fold_aucs.append(auc)
            print(f'    fold {fold}: AUC={auc:.3f}', flush=True)
        m, s = float(np.mean(fold_aucs)), float(np.std(fold_aucs))
        rows.append({
            'method': f'PCA-{n_comp}+LR',
            'n_components': n_comp,
            'auc_mean': m, 'auc_std': s,
            'fold_aucs': ';'.join(f'{a:.3f}' for a in fold_aucs),
        })
        print(f'    mean: {m:.3f} +/- {s:.3f}', flush=True)

    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'pca_baseline.csv'), index=False)
    print(f'\nSaved: pca_baseline.csv', flush=True)


# ============================================================================
# Stage 4: LOBO (leave-one-cohort-out)
# ============================================================================

def stage4_lobo():
    print('\n' + '=' * 70, flush=True)
    print('STAGE 4: LOBO (leave-one-cohort-out)', flush=True)
    print('=' * 70, flush=True)

    rows = []
    for fm, fname, key in _selected_embeddings():
        path = os.path.join(ROOT, fname)
        if not os.path.exists(path):
            continue
        print(f'\n--- {fm} ---', flush=True)
        D = load_filtered(path, key)
        if D['cohort'] is None:
            print('  SKIP: no Data_NO column', flush=True); continue
        X, sev, sid, coh = D['X'], D['severity'], D['sample_id'], D['cohort']
        y = (sev == 'fatal').astype(int)

        unique_cohorts = sorted(set(coh))
        for c in unique_cohorts:
            test_mask = coh == c
            train_mask = ~test_mask
            n_test = test_mask.sum()
            n_train = train_mask.sum()
            n_test_fatal = int(y[test_mask].sum())
            n_test_surv = n_test - n_test_fatal
            if n_test_fatal == 0 or n_test_surv == 0 or n_train == 0:
                print(f'  cohort={c}: skip (insufficient class balance)', flush=True)
                continue
            clf = lr_classifier()
            clf.fit(X[train_mask], y[train_mask])
            p = clf.predict_proba(X[test_mask])[:, 1]
            yh = clf.predict(X[test_mask])
            auc = roc_auc_score(y[test_mask], p)
            bacc = balanced_accuracy_score(y[test_mask], yh)
            rows.append({
                'fm': fm, 'held_cohort': c,
                'n_train': int(n_train), 'n_test': int(n_test),
                'n_test_fatal': n_test_fatal, 'n_test_survived': n_test_surv,
                'auc': auc, 'balanced_accuracy': bacc,
            })
            print(f'  cohort={c}: AUC={auc:.3f}, BAcc={bacc:.3f} (n_test={n_test})', flush=True)

    if _ONLY_FM is not None:
        os.makedirs(PART_DIR, exist_ok=True)
        pd.DataFrame(rows).to_csv(os.path.join(PART_DIR, f'lobo_{_ONLY_FM}.csv'), index=False)
        print(f'\nSaved partial lobo for {_ONLY_FM}', flush=True)
    else:
        pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'lobo.csv'), index=False)
        print(f'\nSaved: lobo.csv', flush=True)


# ============================================================================
# Stage 5: Per-cohort OOF breakdown (uses Stage 1 OOF predictions)
# ============================================================================

def stage5_per_cohort_oof():
    print('\n' + '=' * 70, flush=True)
    print('STAGE 5: Per-cohort OOF breakdown', flush=True)
    print('=' * 70, flush=True)

    # Load cohort labels per cell (from any FM h5ad — same source)
    cohort_path = os.path.join(ROOT, 'state_covid_canonical.h5ad')
    with h5py.File(cohort_path, 'r') as h:
        disease = read_obs(h, 'Disease')
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
        cell_id = read_obs(h, 'Cell_ID')
        cohort = read_obs(h, 'Data_NO')
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    cohort_lookup = pd.DataFrame({
        'cell_id': cell_id[keep] if cell_id is not None else np.arange(int(keep.sum())),
        'cohort': cohort[keep],
    }).astype(str)

    rows = []
    for fm, _, _ in EMBEDDINGS:
        oof_p = os.path.join(OUT_DIR, 'oof_predictions', f'{fm}.csv')
        if not os.path.exists(oof_p):
            continue
        oof = pd.read_csv(oof_p, dtype={'cell_id': str})
        merged = oof.merge(cohort_lookup, on='cell_id', how='left')
        for c in sorted(merged['cohort'].dropna().unique()):
            sub = merged[merged['cohort'] == c]
            yv = sub['y'].values
            pv = sub['p_lr'].values
            mask_valid = ~np.isnan(pv)
            if mask_valid.sum() < 5 or yv[mask_valid].sum() == 0 or (1 - yv[mask_valid]).sum() == 0:
                continue
            auc = roc_auc_score(yv[mask_valid], pv[mask_valid])
            rows.append({
                'fm': fm, 'cohort': c,
                'n_cells': int(mask_valid.sum()),
                'auc': auc,
            })
            print(f'  {fm} | cohort={c}: AUC={auc:.3f} (n={mask_valid.sum()})', flush=True)

    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'per_cohort_oof.csv'), index=False)
    print(f'\nSaved: per_cohort_oof.csv', flush=True)


# ============================================================================
# Main
# ============================================================================

def _concat_partials(prefix, out_name):
    """Concatenate per-FM partials in EMBEDDINGS order into the aggregated CSV.

    Byte-identical to the monolithic write: each partial holds that FM's rows
    formatted by pandas exactly as they appear in the full file (CSV formats each
    value independently of the others), and we re-emit them in EMBEDDINGS order.
    """
    frames = []
    for fm, _, _ in EMBEDDINGS:
        p = os.path.join(PART_DIR, f'{prefix}_{fm}.csv')
        if os.path.exists(p):
            frames.append(pd.read_csv(p))
    if not frames:
        print(f'  WARN: no partials found for {prefix}', flush=True)
        return
    pd.concat(frames, ignore_index=True).to_csv(os.path.join(OUT_DIR, out_name), index=False)
    print(f'  gathered {out_name} from {len(frames)} per-FM partials', flush=True)


def run_gather():
    print('\n' + '=' * 70, flush=True)
    print('GATHER: concat per-FM partials + DeLong + per-cohort', flush=True)
    print('=' * 70, flush=True)
    _concat_partials('mortality', 'mortality_lr_rf.csv')
    _concat_partials('calibration', 'calibration.csv')
    _concat_partials('lobo', 'lobo.csv')
    stage2_delong()       # reads all oof_predictions/<FM>.csv (cross-FM)
    stage5_per_cohort_oof()
    print('\nGATHER DONE.', flush=True)


def main():
    global _ONLY_FM
    ap = argparse.ArgumentParser(description='COVID+Sepsis-only re-run (monolithic or array+gather)')
    ap.add_argument('--mode', choices=['all', 'per_fm', 'pca', 'gather'], default='all',
                    help="all=original single run; per_fm=one FM's stage1+stage4 partials; "
                         "pca=stage3 only; gather=concat partials + DeLong + per-cohort")
    ap.add_argument('--index', type=int, default=None,
                    help='per_fm: index into EMBEDDINGS; falls back to $SLURM_ARRAY_TASK_ID')
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    print(f'COVID + Sepsis-only re-run | mode={args.mode}', flush=True)
    print(f'  ROOT      = {ROOT}', flush=True)
    print(f'  OUT_DIR   = {OUT_DIR}', flush=True)

    if args.mode == 'all':
        stage1_primary_calibration()
        stage2_delong()
        stage3_pca_baseline()
        stage4_lobo()
        stage5_per_cohort_oof()
        print('\nALL DONE.', flush=True)
    elif args.mode == 'per_fm':
        idx = args.index
        if idx is None:
            idx = int(os.environ.get('SLURM_ARRAY_TASK_ID', '0'))
        _ONLY_FM = EMBEDDINGS[idx][0]
        print(f'  per_fm: index={idx} fm={_ONLY_FM}', flush=True)
        stage1_primary_calibration()
        stage4_lobo()
        print(f'\nPER-FM DONE ({_ONLY_FM}).', flush=True)
    elif args.mode == 'pca':
        stage3_pca_baseline()
        print('\nPCA DONE.', flush=True)
    elif args.mode == 'gather':
        run_gather()


if __name__ == '__main__':
    main()
