"""
Wires `cv_group` into the SECONDARY evaluation protocol (grouped 10-seed x
5-fold StratifiedGroupKFold on the pooled set) for the rebuilt benchmark --
the leftover TODO flagged since Phase 2 (benchmark_design_decisions.md
"Leakage design"; project_302006_rebuild_architecture memory). LOCO (the
PRIMARY evaluation, see pca200_loco_before_after_decontx.py) needs no such
wiring -- cohort-level splits can never break a well/pool apart. The
secondary pooled evaluation CAN: a plain patient_id group-key would let
COMBAT's single 124-donor block (all sharing one multiplexed pool) split
across folds, which is exactly the ambient-RNA/doublet leakage this
`cv_group` column exists to prevent (cv_group = connected components of the
"shares a GEM well" graph, see canonical_donor_labels.csv).

Smoke-tested here (not a benchmark run) on the Platelet-only subset of the
final 8-cohort assembled dataset (665 donors, 503 cv_groups, 80 events),
using the same TruncatedSVD-200 + LogisticRegression baseline validated in
pca200_loco_before_after_decontx.py -- Platelet-only because it is the
project's actual thesis subset and small enough (112,967 cells) to smoke-
test quickly; swapping in real FM embeddings later needs no change to this
script's harness logic, only to what produces X.

repeated_cv.py's run_repeated_kfold needs no code change for this -- its
`groups` argument is already caller-supplied; this script is the first
caller in the rebuild to actually pass `groups=cv_group` instead of
`groups=patient_id`.

Usage:
    python pca200_pooled_secondary_cv.py --which raw|decontx --lineage Platelet|Monocyte_Myeloid|all
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code')
from repeated_cv import run_repeated_kfold, repeated_patient_auc, summarize_repeated
sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code')
from assemble_platelet_monocyte_benchmark_input import read_mtx_columns

DATA_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/pooled_secondary_cv'
N_COMPONENTS = 200


def fit_predict(X_tr, y_tr, X_te):
    scaler = StandardScaler(with_mean=False).fit(X_tr)
    Xs_tr, Xs_te = scaler.transform(X_tr), scaler.transform(X_te)
    svd = TruncatedSVD(n_components=min(N_COMPONENTS, min(Xs_tr.shape) - 1), random_state=42)
    Z_tr, Z_te = svd.fit_transform(Xs_tr), svd.transform(Xs_te)
    clf = LogisticRegression(C=1.0, solver='lbfgs', max_iter=2000, class_weight='balanced', n_jobs=4)
    clf.fit(Z_tr, y_tr)
    return {'p_lr': clf.predict_proba(Z_te)[:, 1]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--which', required=True, choices=['raw', 'decontx'])
    ap.add_argument('--lineage', default='Platelet', choices=['Platelet', 'Monocyte_Myeloid', 'all'])
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f'Loading {args.which} cells.csv...', flush=True)
    cells_full = pd.read_csv(f'{DATA_DIR}/{args.which}/cells.csv', low_memory=False)
    assert 'cv_group' in cells_full.columns, 'cv_group missing -- run the cv_group backfill first'

    # event is NaN for donors with a genuinely undefined outcome (e.g.
    # Stephenson's "unknown"-outcome patients) -- real bug caught here
    # 2026-09-19: casting NaN to int silently produces garbage labels
    # instead of raising, so this must be filtered BEFORE the cast, not
    # assumed absent. 50,588 / 1,236,881 cells (53 donors across 5 cohorts)
    # have no defined outcome and cannot contribute to this classification
    # task at all.
    keep = cells_full['event'].notna()
    if args.lineage != 'all':
        keep &= cells_full['harmonized_lineage'] == args.lineage
    n_undefined = (cells_full['event'].isna() & (
        cells_full['harmonized_lineage'] == args.lineage if args.lineage != 'all'
        else pd.Series(True, index=cells_full.index))).sum()
    idx = np.where(keep.to_numpy())[0]
    cells = cells_full.iloc[idx].reset_index(drop=True)
    print(f'  excluded {n_undefined} cells with undefined event (unknown outcome)', flush=True)

    # Stream the (22-50G) combined counts.mtx, keeping only the wanted
    # lineage's columns while reading -- same reasoning as
    # assemble_platelet_monocyte_benchmark_input.py's read_mtx_columns: a
    # plain mmread-then-subset on a file this size risks the same memory/
    # walltime blowup that hit Chongqing's DecontX write step.
    print(f'Streaming {args.which} counts.mtx, filtering to {len(idx)} {args.lineage} cells...',
          flush=True)
    mat = read_mtx_columns(f'{DATA_DIR}/{args.which}/counts.mtx', idx)  # genes x cells
    print(f'  {mat.shape[1]} cells x {mat.shape[0]} genes, {cells.donor_id.nunique()} donors, '
          f'{cells.cv_group.nunique()} cv_groups', flush=True)

    X = mat.T.tocsr()
    y = cells['event'].to_numpy().astype(int)
    patient_id = cells['donor_id'].to_numpy()
    cv_group = cells['cv_group'].to_numpy()

    # Direct, independent proof that no cv_group is ever split across a fold
    # (the actual property this script exists to guarantee) -- checked for
    # every seed, not assumed from StratifiedGroupKFold's contract alone.
    print('Verifying cv_group purity across all folds/seeds...', flush=True)
    for seed in [42, 47, 51]:
        skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        group_to_folds = {}
        for fold_idx, (_, te) in enumerate(skf.split(X, y, groups=cv_group)):
            for g in np.unique(cv_group[te]):
                group_to_folds.setdefault(g, set()).add(fold_idx)
        bad = {g: f for g, f in group_to_folds.items() if len(f) > 1}
        assert not bad, f'seed {seed}: {len(bad)} cv_groups split across folds: {list(bad)[:5]}'
    print(f'  Verified: no cv_group split across folds for seeds 42/47/51 '
          f'({cells.cv_group.nunique()} groups checked each time)', flush=True)

    print('Running grouped 10-seed x 5-fold CV (groups=cv_group)...', flush=True)
    oof_by_seed = run_repeated_kfold(X, y, groups=cv_group, fit_predict_fn=fit_predict)

    aucs, seeds = repeated_patient_auc(oof_by_seed, patient_id, y, score_key='p_lr')
    headline = summarize_repeated(aucs)
    print(f'\n=== Pooled secondary CV (grouped by cv_group), lineage={args.lineage}, {args.which} ===')
    print(f'Patient-level AUC across {len(seeds)} seeds: mean={headline["mean"]:.4f} '
          f'sd={headline["sd"]:.4f} (n={headline["n"]})')
    print(f'Per-seed AUCs: {dict(zip(seeds, np.round(aucs, 4)))}')

    report = pd.DataFrame({'seed': seeds, 'auc': aucs})
    out_path = f'{OUT_DIR}/pooled_secondary_cv_{args.which}_{args.lineage}.csv'
    report.to_csv(out_path, index=False)
    print(f'\nSAVED -> {out_path}')


if __name__ == '__main__':
    main()
