#!/usr/bin/env python3
"""
Shared repeated-CV utilities -- CV Protocol v2 (2026-09-04 user decision, see
HPCC_RUNBOOK.md "CV Protocol v2" section and reports/manuscript/revision/
P0_onepager.md Sec.2).

Why: a fold-composition check found StratifiedGroupKFold(shuffle=True) reshuffles
~80% of patients across folds after a 2-3 label correction -- a single train/test
split's AUC ranking is dominated by which patients happened to land in the test
fold, not by real model differences. This is now the PRIMARY protocol for every
CV-dependent step in the 302006 revision: Job H2 (within-COVID CV, PCA-200
baseline, null model), Job S (within-sepsis CV, transfer, adaptation curve),
and Job F's own evaluation.

Protocol: 10 seeds x 5-fold patient-level StratifiedGroupKFold = 50 fits/model.
Reporting spec (apply exactly, in this order):
  1. Per-repeat pooled-OOF patient-level AUC -> mean +/- sd across the 10
     repeats. This is the headline number.
  2. Pairwise dAUC, SAME seed paired across models -> mean delta + 95% CI
     across the 10 seeds (t-based), Wilcoxon signed-rank as secondary.
  3. The seed-42 single split is kept only as a labeled secondary view, never
     as the headline.

Usage pattern for a CV-fitting script (e.g. rerun_covid_sepsis_mortality_patient_v2.py):
    from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold
    def fit_predict(X_tr, y_tr, X_te):
        lr = LogisticRegression(...).fit(X_tr, y_tr)
        rf = RandomForestClassifier(...).fit(X_tr, y_tr)
        return {'p_lr': lr.predict_proba(X_te)[:, 1], 'p_rf': rf.predict_proba(X_te)[:, 1]}
    oof_by_seed = run_repeated_kfold(X, y, groups=patient_id, fit_predict_fn=fit_predict,
                                      seeds=seeds_arg)
    # then write one oof_<name>_seed<NN>.csv per seed (schema unchanged per-seed,
    # just one file per seed instead of one file total)

Usage pattern for an aggregation/reporting script (e.g. donor_level_reanalysis.py):
    from repeated_cv import repeated_patient_auc, summarize_repeated, paired_delta_auc
    aucs_a, seeds = repeated_patient_auc(oof_by_seed_a, patient_id, y, score_key='p_lr')
    aucs_b, _     = repeated_patient_auc(oof_by_seed_b, patient_id, y, score_key='p_lr')
    headline = summarize_repeated(aucs_a)               # {'mean':..., 'sd':..., 'n':10}
    pair     = paired_delta_auc(aucs_a, aucs_b)          # same-seed-paired dAUC + CI + Wilcoxon

---

LOCO protocol (benchmark_design_decisions.md, 2026-09-15 rebuild design) -- the PRIMARY
evaluation for the rebuilt benchmark, alongside (not replacing) the StratifiedGroupKFold
protocol above, which becomes the SECONDARY evaluation once grouped on `cv_group` instead
of `patient_id` (cv_group = connected components of the "shares a GEM well" graph; a plain
patient_id group-key would let COMBAT's single 124-donor block split across folds, which is
exactly the leak this redesign exists to prevent -- see project_302006_rebuild_architecture
memory). `run_repeated_kfold`'s `groups` argument is already caller-supplied and needs no
change for this: once `canonical_donor_labels.csv` carries `cv_group`, callers just pass
`groups=cv_group` instead of `groups=patient_id`.

LOCO itself has no seed/shuffle (the partition is fully determined by cohort membership) and
no meaningful stratification choice -- each held-out fold's event rate is whatever that
cohort's real rate is. Every cohort is held out exactly once by construction, so no group can
ever split across a LOCO train/test boundary regardless of well structure.

Usage pattern:
    from repeated_cv import run_loco, summarize_loco
    oof, held_out_cohorts = run_loco(X, y, cohort, fit_predict_fn)
    report = summarize_loco(oof, patient_id, y, cohort, score_key='p_lr')
    # -> DataFrame: one row per cohort (patient-level AUC + n + n_events),
    #    plus a final 'POOLED' row (pooled patient-level AUC, not an average
    #    of the per-cohort rows)
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, LeaveOneGroupOut

DEFAULT_SEEDS = list(range(42, 52))  # 42..51 inclusive, 10 seeds (42 = the historical single-split seed)
N_SPLITS = 5


def safe_auc(y, p):
    y = np.asarray(y)
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, p)


def run_repeated_kfold(X, y, groups, fit_predict_fn, seeds=None, n_splits=N_SPLITS):
    """Repeated StratifiedGroupKFold OOF predictions.

    fit_predict_fn(X_train, y_train, X_test) -> dict {name: proba_test_array}
    (e.g. {'p_lr': ..., 'p_rf': ...}); every seed must return the same set of
    names.

    Returns: dict seed -> dict name -> full-length OOF array (float), row order
    matching the input X/y/groups. With group-k-fold every row is held out
    exactly once per seed, so there should be no NaNs left; callers should
    assert this (a NaN indicates a group had too few members for n_splits).
    """
    if seeds is None:
        seeds = DEFAULT_SEEDS
    y = np.asarray(y)
    groups = np.asarray(groups)
    n = len(y)
    results = {}
    for seed in seeds:
        skf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        preds = None
        for tr, te in skf.split(X, y, groups=groups):
            out = fit_predict_fn(X[tr] if not hasattr(X, 'iloc') else X.iloc[tr],
                                  y[tr], X[te] if not hasattr(X, 'iloc') else X.iloc[te])
            if preds is None:
                preds = {k: np.full(n, np.nan) for k in out}
            for k, v in out.items():
                preds[k][te] = v
        n_missing = {k: int(np.isnan(v).sum()) for k, v in preds.items()}
        if any(n_missing.values()):
            raise AssertionError(f'seed {seed}: rows never held out {n_missing}')
        results[seed] = preds
    return results


def run_loco(X, y, cohort, fit_predict_fn):
    """Leave-one-cohort-out OOF predictions -- see module docstring "LOCO protocol".

    fit_predict_fn(X_train, y_train, X_test) -> dict {name: proba_test_array},
    same contract as run_repeated_kfold. X/y/cohort are at whatever row
    granularity the caller works in (per-cell, matching run_repeated_kfold's
    usual usage) -- this function does not compute AUC itself, deliberately:
    every other headline number in this module is *patient*-level (module
    docstring "Reporting spec"), which requires the separate
    aggregate_patient_level step: cell-level AUC here would silently be on
    the wrong granularity, one cell-count away from every other number this
    pipeline reports. Use summarize_loco for the patient-level per-cohort
    table.

    Returns:
      oof: dict name -> full-length OOF array (float), row order matching
        the input X/y/cohort. Every row held out exactly once, so there
        should be no NaNs left (asserted, same as run_repeated_kfold).
      held_out_cohorts: list of the cohort labels in the order LeaveOneGroupOut
        produced them (informational only -- summarize_loco re-derives its own
        grouping from patient_id/cohort, it does not depend on this order).
    """
    y = np.asarray(y)
    cohort = np.asarray(cohort)
    n = len(y)
    logo = LeaveOneGroupOut()
    oof = None
    held_out_cohorts = []
    for tr, te in logo.split(X, y, groups=cohort):
        held_out_cohorts.append(cohort[te][0])
        out = fit_predict_fn(X[tr] if not hasattr(X, 'iloc') else X.iloc[tr],
                              y[tr], X[te] if not hasattr(X, 'iloc') else X.iloc[te])
        if oof is None:
            oof = {k: np.full(n, np.nan) for k in out}
        for k, v in out.items():
            oof[k][te] = v
    n_missing = {k: int(np.isnan(v).sum()) for k, v in oof.items()}
    if any(n_missing.values()):
        raise AssertionError(f'LOCO: rows never held out {n_missing}')
    return oof, held_out_cohorts


def summarize_loco(oof, patient_id, y, cohort, score_key='p_lr', agg='mean'):
    """Patient-level per-cohort AUC table (DataFrame, sorted by n_events
    descending) plus the pooled patient-level AUC across all cohorts combined
    reported as its own final row ('cohort'=='POOLED'). Pooled is NOT an
    average of the per-cohort AUCs -- that would misweight cohorts this
    different in size (COMBAT n=124/16 events vs. Qiu n=7/2), so it is
    computed from the full pooled patient-level OOF instead, same as the rest
    of this module's pooled-OOF convention.

    patient_id/y/cohort must be the same row-length/order as the oof arrays
    from run_loco (i.e. the original per-cell X/y/cohort/patient_id, not
    already patient-aggregated) -- aggregation to patient-level happens here,
    once, per cohort.
    """
    patient_id = np.asarray(patient_id)
    cohort = np.asarray(cohort)
    p = oof[score_key]
    rows = []
    for c in sorted(set(cohort)):
        mask = cohort == c
        pat = aggregate_patient_level(patient_id[mask], np.asarray(y)[mask], p[mask], agg=agg)
        rows.append({'cohort': str(c), 'n': len(pat), 'n_events': int(pat['y'].sum()),
                     'auc': safe_auc(pat['y'].to_numpy(), pat['score'].to_numpy())})
    pat_pooled = aggregate_patient_level(patient_id, y, p, agg=agg)
    rows.append({'cohort': 'POOLED', 'n': len(pat_pooled), 'n_events': int(pat_pooled['y'].sum()),
                 'auc': safe_auc(pat_pooled['y'].to_numpy(), pat_pooled['score'].to_numpy())})
    df = pd.DataFrame(rows)
    pooled_row = df[df['cohort'] == 'POOLED']
    df = pd.concat([df[df['cohort'] != 'POOLED'].sort_values('n_events', ascending=False),
                    pooled_row]).reset_index(drop=True)
    return df


def aggregate_patient_level(patient_id, y, score, agg='mean'):
    """One row per patient_id: aggregated score + patient label (max y over draws)."""
    df = pd.DataFrame({'patient_id': patient_id, 'y': y, 'score': score})
    g = df.groupby('patient_id')
    if agg == 'mean':
        s = g['score'].mean()
    elif agg == 'median':
        s = g['score'].median()
    else:
        raise ValueError(f'unknown agg {agg}')
    return pd.DataFrame({'y': g['y'].max(), 'score': s})


def repeated_patient_auc(oof_by_seed, patient_id, y, score_key='p_lr', agg='mean'):
    """Per-seed patient-level pooled-OOF AUC.

    oof_by_seed: dict seed -> {score_key: full-length array}, same row order
    as patient_id/y.
    Returns (aucs: np.ndarray shape (n_seeds,), seeds: sorted list) -- aucs[i]
    corresponds to seeds[i], so two calls on different models with the SAME
    oof_by_seed key set can be paired index-for-index for paired_delta_auc.
    """
    seeds = sorted(oof_by_seed.keys())
    aucs = []
    for s in seeds:
        p = oof_by_seed[s][score_key]
        pat = aggregate_patient_level(patient_id, y, p, agg=agg)
        aucs.append(safe_auc(pat['y'].to_numpy(), pat['score'].to_numpy()))
    return np.array(aucs), seeds


def summarize_repeated(aucs):
    """Headline reporting stat #1: repeat-mean +/- sd."""
    aucs = np.asarray(aucs, dtype=float)
    valid = aucs[~np.isnan(aucs)]
    return {'mean': float(np.mean(valid)) if len(valid) else float('nan'),
            'sd': float(np.std(valid, ddof=1)) if len(valid) > 1 else float('nan'),
            'n': int(len(valid))}


def paired_delta_auc(aucs_a, aucs_b):
    """Reporting stat #2: same-seed-paired dAUC, mean + 95% CI (t-based) + Wilcoxon.

    aucs_a, aucs_b: arrays of per-seed AUC for two models, SAME seed order
    (i.e. both produced by repeated_patient_auc with the same seed list).
    """
    aucs_a = np.asarray(aucs_a, dtype=float)
    aucs_b = np.asarray(aucs_b, dtype=float)
    mask = ~(np.isnan(aucs_a) | np.isnan(aucs_b))
    delta = aucs_a[mask] - aucs_b[mask]
    n = len(delta)
    if n < 2:
        return {'mean_delta': float('nan'), 'ci_lo': float('nan'), 'ci_hi': float('nan'),
                'wilcoxon_p': float('nan'), 'n_seeds': n}
    mean_d = float(delta.mean())
    sd_d = float(delta.std(ddof=1))
    se = sd_d / np.sqrt(n)
    tcrit = stats.t.ppf(0.975, df=n - 1)
    ci_lo, ci_hi = mean_d - tcrit * se, mean_d + tcrit * se
    try:
        _, w_p = stats.wilcoxon(delta)
    except ValueError:
        w_p = float('nan')
    return {'mean_delta': mean_d, 'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
            'wilcoxon_p': float(w_p), 'n_seeds': n}
