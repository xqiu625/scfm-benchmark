#!/usr/bin/env python3
"""Primary benchmark suite on Option C 5-FM set (in-hospital mortality binary).

Five analyses run in sequence per FM:
  1. Mortality benchmark      → 5-fold StratifiedGroupKFold(Sample_ID), LR + RF
                                Outputs per-fold + mean AUC, BalAcc
  2. Clinical metrics         → AUC-ROC, AUC-PR, Sens@90%Spec, Cohen's kappa
                                Computed from out-of-fold predictions
  3. DeLong pairwise tests    → Across all 5 FMs on pooled OOF predictions
  4. Sample efficiency        → train on 5/10/20/50/100% of training set
                                (subsample at the patient level to preserve groups)
  5. LOBO cross-batch         → leave-one-Data_NO-out generalization

Skipped: gene dropout (requires re-embedding through each FM — separate effort).

Outputs (in --output-dir):
  mortality_per_fold.csv          # FM × classifier × fold → AUC, BalAcc
  mortality_summary.csv           # FM × classifier → mean AUC, std, BalAcc
  clinical_metrics.csv            # FM → AUC, AUC-PR, Sens@90Spec, Kappa
  delong_pairwise.csv             # FM_a, FM_b → AUC_diff, p-value, sig
  sample_efficiency.csv           # FM × frac → AUC, BalAcc
  lobo.csv                        # FM × cohort_held_out → AUC
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import norm
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    cohen_kappa_score,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

# Embedding sanity check (catches degenerate FM outputs like the May 2026 TF-Metazoa
# tokenization-mismatch bug). Path-relative import so this works regardless of cwd.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
from embedding_sanity import check_embedding_or_raise, EmbeddingSanityError  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("primary_benchmark")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer", "embeddings"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}

N_FOLDS = 5
SAMPLE_FRACS = [0.05, 0.10, 0.20, 0.50, 1.00]


def find_embedding(adata, model_key: str):
    """Find an FM embedding by known obsm keys. Refuses the previous unsafe
    fallback (first 2D obsm) — that could silently grab UMAP coords or other
    auxiliary arrays. Raises a clear error if no known key is present."""
    keys = OBSM_KEYS.get(model_key, [])
    for k in keys:
        if k in adata.obsm:
            return np.asarray(adata.obsm[k]), k
    raise RuntimeError(
        f"No known embedding key for model_key={model_key!r}. "
        f"Tried {keys}; available obsm: {list(adata.obsm.keys())}. "
        f"If a new FM is added, register its obsm key(s) in OBSM_KEYS."
    )


def build_mortality(severity: np.ndarray):
    sev = np.asarray([str(s).lower() for s in severity])
    fatal = sev == "fatal"
    survived = np.isin(sev, ["mild", "moderate", "recovered", "severe"])
    valid = fatal | survived
    y = np.zeros(len(sev), dtype=np.int8)
    y[fatal] = 1
    return y, valid


def fit_classifier(X_tr, y_tr, X_te, kind="lr"):
    """Return (proba_pos, pred). Standardize inside the fold (no leakage)."""
    if kind == "lr":
        clf = LogisticRegression(
            penalty="l2", C=1.0, solver="lbfgs",
            max_iter=2000, class_weight="balanced", n_jobs=4,
        )
        scaler = StandardScaler().fit(X_tr)
        clf.fit(scaler.transform(X_tr), y_tr)
        proba = clf.predict_proba(scaler.transform(X_te))[:, 1]
        pred = clf.predict(scaler.transform(X_te))
    elif kind == "rf":
        clf = RandomForestClassifier(
            n_estimators=100, n_jobs=8, class_weight="balanced", random_state=42,
        )
        clf.fit(X_tr, y_tr)
        proba = clf.predict_proba(X_te)[:, 1]
        pred = clf.predict(X_te)
    else:
        raise ValueError(kind)
    return proba, pred


def sens_at_spec(y_true, y_score, target_spec=0.90):
    """Sensitivity at fixed specificity."""
    fpr, tpr, _ = roc_curve(y_true, y_score)
    spec = 1 - fpr
    valid = spec >= target_spec
    if not valid.any():
        return float("nan")
    return float(tpr[valid].max())


def delong_test(y_true, score_a, score_b):
    """DeLong test for two correlated ROC curves (Hanley-McNeil variance).
    Returns (auc_a, auc_b, z, p_two_sided)."""
    y = np.asarray(y_true).astype(int)
    pos = y == 1
    neg = y == 0
    n_pos = int(pos.sum())
    n_neg = int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")

    def _placement(scores, y_eq=1):
        # Midrank version
        from scipy.stats import rankdata
        s_pos = scores[pos]
        s_neg = scores[neg]
        # V10[i] = P(scores_a < neg_score_i) for each pos i (and so on)
        # Use Fast DeLong via rank-based estimator
        return s_pos, s_neg

    auc_a = roc_auc_score(y, score_a)
    auc_b = roc_auc_score(y, score_b)

    # Compute V components via the unbiased estimator (Sun & Xu 2014 fast DeLong)
    # Implementation derived from the standard reference.
    def fast_delong(scores):
        s_pos = scores[pos]
        s_neg = scores[neg]
        # tx, ty are placement values
        from scipy.stats import rankdata
        all_scores = np.concatenate([s_pos, s_neg])
        ranks = rankdata(all_scores, method="average")
        r_pos = ranks[:n_pos]
        r_neg = ranks[n_pos:]
        # V10[i] (pos): ( rank_i - i_within_pos ) / n_neg
        # Using the formula: V10 = (rank_pos_in_all - rank_pos_within_pos) / n_neg
        rank_pos_within = rankdata(s_pos, method="average")
        rank_neg_within = rankdata(s_neg, method="average")
        v10 = (r_pos - rank_pos_within) / n_neg
        v01 = 1 - (r_neg - rank_neg_within) / n_pos
        return v10, v01

    v10_a, v01_a = fast_delong(score_a)
    v10_b, v01_b = fast_delong(score_b)

    s10 = np.cov(np.vstack([v10_a, v10_b]), ddof=1) if n_pos > 1 else np.zeros((2, 2))
    s01 = np.cov(np.vstack([v01_a, v01_b]), ddof=1) if n_neg > 1 else np.zeros((2, 2))
    var = (s10[0, 0] + s10[1, 1] - 2 * s10[0, 1]) / n_pos + \
          (s01[0, 0] + s01[1, 1] - 2 * s01[0, 1]) / n_neg
    if var <= 0:
        return auc_a, auc_b, float("nan"), float("nan")
    z = (auc_a - auc_b) / np.sqrt(var)
    p = 2 * (1 - norm.cdf(abs(z)))
    return auc_a, auc_b, z, p


# ============================================================================
# 1. Mortality benchmark (LR + RF, 5-fold)
# ============================================================================
def benchmark_mortality(model_key: str, X, y, groups):
    """Returns:
      per_fold_rows: list[dict] one per (classifier, fold)
      oof_predictions: dict[classifier -> (oof_proba, oof_pred)]
    """
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    per_fold = []
    oof = {kind: (np.full(len(y), np.nan), np.full(len(y), np.nan)) for kind in ["lr", "rf"]}

    for fold, (tr, te) in enumerate(skf.split(X, y, groups)):
        for kind in ["lr", "rf"]:
            proba, pred = fit_classifier(X[tr], y[tr], X[te], kind=kind)
            auc = roc_auc_score(y[te], proba)
            bacc = balanced_accuracy_score(y[te], pred)
            per_fold.append(dict(model=model_key, classifier=kind, fold=fold,
                                 auc=auc, balanced_accuracy=bacc,
                                 n_train=len(tr), n_test=len(te),
                                 fatal_train=int(y[tr].sum()), fatal_test=int(y[te].sum())))
            oof[kind][0][te] = proba
            oof[kind][1][te] = pred

    return per_fold, oof


# ============================================================================
# 2. Clinical metrics (computed from OOF preds)
# ============================================================================
def clinical_metrics(model_key, classifier, y, oof_proba, oof_pred):
    valid = ~np.isnan(oof_proba)
    y_v = y[valid]
    p_v = oof_proba[valid]
    pred_v = oof_pred[valid]
    return dict(
        model=model_key,
        classifier=classifier,
        auc_roc=roc_auc_score(y_v, p_v),
        auc_pr=average_precision_score(y_v, p_v),
        sens_at_90_spec=sens_at_spec(y_v, p_v, 0.90),
        cohen_kappa=cohen_kappa_score(y_v, pred_v),
    )


# ============================================================================
# 3. DeLong pairwise (LR OOF preds across FMs)
# ============================================================================
def delong_pairwise(y, fm_oof_lr_proba: dict[str, np.ndarray]):
    rows = []
    fms = sorted(fm_oof_lr_proba.keys())
    # Use the cell-set common to all FMs (intersect non-NaN positions)
    valid_mask = np.ones(len(y), dtype=bool)
    for k, p in fm_oof_lr_proba.items():
        valid_mask &= ~np.isnan(p)
    y_v = y[valid_mask]
    log.info("DeLong: %d cells common across FMs (of %d)", valid_mask.sum(), len(y))
    for i, fa in enumerate(fms):
        for fb in fms[i + 1:]:
            pa = fm_oof_lr_proba[fa][valid_mask]
            pb = fm_oof_lr_proba[fb][valid_mask]
            auc_a, auc_b, z, p = delong_test(y_v, pa, pb)
            rows.append(dict(model_a=fa, model_b=fb,
                             auc_a=auc_a, auc_b=auc_b,
                             auc_diff=auc_a - auc_b, z=z, p_value=p,
                             significant=(p < 0.05) if not np.isnan(p) else False))
    return rows


# ============================================================================
# 4. Sample efficiency
# ============================================================================
def sample_efficiency(model_key, X, y, groups):
    rng = np.random.RandomState(42)
    rows = []
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)

    for fold, (tr_full, te) in enumerate(skf.split(X, y, groups)):
        # Group-level subsample within training set
        tr_groups_unique = np.unique(groups[tr_full])
        for frac in SAMPLE_FRACS:
            n_groups_keep = max(2, int(len(tr_groups_unique) * frac))
            keep_groups = rng.choice(tr_groups_unique, size=n_groups_keep, replace=False)
            tr_mask = np.isin(groups[tr_full], keep_groups)
            tr_sub = tr_full[tr_mask]
            if (y[tr_sub] == 1).sum() < 5 or (y[tr_sub] == 0).sum() < 5:
                continue
            proba, pred = fit_classifier(X[tr_sub], y[tr_sub], X[te], kind="lr")
            rows.append(dict(model=model_key, fold=fold, train_frac=frac,
                             n_train=len(tr_sub), n_train_groups=n_groups_keep,
                             auc=roc_auc_score(y[te], proba),
                             balanced_accuracy=balanced_accuracy_score(y[te], pred)))
    return rows


# ============================================================================
# 5. LOBO (Leave-One-Batch-Out, batch = Data_NO)
# ============================================================================
def lobo(model_key, X, y, cohorts, groups):
    rows = []
    cohort_unique = np.unique(cohorts)
    for held_cohort in cohort_unique:
        te_mask = cohorts == held_cohort
        tr_mask = ~te_mask
        n_tr = int(tr_mask.sum())
        n_te = int(te_mask.sum())
        n_te_fatal = int(y[te_mask].sum())
        n_te_survived = int((y[te_mask] == 0).sum())
        if n_te_fatal < 5 or n_te_survived < 5 or n_tr < 100:
            log.info("  LOBO skip cohort=%s (insufficient class balance)", held_cohort)
            continue
        proba, pred = fit_classifier(X[tr_mask], y[tr_mask], X[te_mask], kind="lr")
        rows.append(dict(model=model_key, held_cohort=held_cohort,
                         n_train=n_tr, n_test=n_te,
                         n_test_fatal=n_te_fatal, n_test_survived=n_te_survived,
                         auc=roc_auc_score(y[te_mask], proba),
                         balanced_accuracy=balanced_accuracy_score(y[te_mask], pred)))
    return rows


# ============================================================================
# Main
# ============================================================================
def per_fm(fm_h5ad: Path, output_dir: Path, analyses=None):
    if analyses is None:
        analyses = {"mortality", "clinical", "sample_efficiency", "lobo"}
    model_key = fm_h5ad.stem.split("_")[0].lower()
    log.info("=" * 60)
    log.info("Model: %s (%s)", model_key, fm_h5ad.name)
    log.info("=" * 60)

    adata = sc.read_h5ad(fm_h5ad)
    X, obsm_key = find_embedding(adata, model_key)
    log.info("  obsm=%s shape=%s", obsm_key, X.shape)

    # Sanity check: catches degenerate embeddings (e.g., constant-vector output
    # from a tokenization mismatch — see TF-Metazoa near-miss May 2026).
    try:
        check_embedding_or_raise(X, fm_label=model_key)
    except EmbeddingSanityError as e:
        log.error("Skipping %s — embedding failed sanity check.\n%s", model_key, e)
        return None, None

    if "covid_severity" not in adata.obs.columns:
        return None, None
    if "Sample_ID" not in adata.obs.columns:
        return None, None
    cohort_col = "Data_NO" if "Data_NO" in adata.obs.columns else None

    y_full, valid = build_mortality(adata.obs["covid_severity"].astype(str).values)
    groups = adata.obs["Sample_ID"].astype(str).values
    cohorts = adata.obs[cohort_col].astype(str).values if cohort_col else None
    X_v, y_v, g_v = X[valid], y_full[valid], groups[valid]
    c_v = cohorts[valid] if cohorts is not None else None
    log.info("  CV input: n=%d, fatal=%d, survived=%d, groups=%d",
             len(y_v), int(y_v.sum()), int((y_v == 0).sum()), len(np.unique(g_v)))
    if c_v is not None:
        log.info("  cohorts: %d unique (%s)", len(np.unique(c_v)), cohort_col)

    # 1. Mortality benchmark (also feeds clinical metrics and DeLong OOF)
    fold_rows, oof = [], None
    if analyses & {"mortality", "clinical", "delong"}:
        log.info("[1/5] Mortality benchmark")
        fold_rows, oof = benchmark_mortality(model_key, X_v, y_v, g_v)

    # 2. Clinical metrics (from OOF)
    clin_rows = []
    if "clinical" in analyses:
        log.info("[2/5] Clinical metrics")
        for kind, (proba, pred) in oof.items():
            clin_rows.append(clinical_metrics(model_key, kind, y_v, proba, pred))

    # 3. Sample efficiency
    se_rows = []
    if "sample_efficiency" in analyses:
        log.info("[3/5] Sample efficiency")
        se_rows = sample_efficiency(model_key, X_v, y_v, g_v)

    # 4. LOBO
    lobo_rows = []
    if "lobo" in analyses:
        log.info("[4/5] LOBO")
        lobo_rows = lobo(model_key, X_v, y_v, c_v, g_v) if c_v is not None else []

    lr_oof_proba = oof["lr"][0] if oof is not None else np.full(len(y_v), np.nan)
    return {
        "fold": fold_rows,
        "clinical": clin_rows,
        "sample_eff": se_rows,
        "lobo": lobo_rows,
    }, {"y": y_v, "lr_oof_proba": lr_oof_proba}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--models", nargs="+", default=None,
                   help="Subset of FMs to run, by model_key (e.g. --models scgpt). "
                        "Default: all *_covid_canonical.h5ad in --embedding-dir.")
    p.add_argument("--analyses", nargs="+", default=None,
                   choices=["mortality", "clinical", "delong", "sample_efficiency", "lobo"],
                   help="Subset of analyses to run (default: all five). CSVs are "
                        "written only for analyses that produced rows.")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    analyses = set(args.analyses) if args.analyses else {"mortality", "clinical", "delong",
                                                         "sample_efficiency", "lobo"}

    fm_files = sorted(args.embedding_dir.glob("*_covid_canonical.h5ad"))
    if args.models:
        wanted = {m.lower() for m in args.models}
        fm_files = [f for f in fm_files if f.stem.split("_")[0].lower() in wanted]
    if not fm_files:
        log.error("No *_covid_canonical.h5ad in %s", args.embedding_dir)
        sys.exit(1)
    log.info("Found %d FM files", len(fm_files))

    all_fold = []
    all_clin = []
    all_se = []
    all_lobo = []
    fm_oof = {}  # model_key -> {y, lr_oof_proba}

    for fm in fm_files:
        try:
            tables, oof_data = per_fm(fm, args.output_dir, analyses=analyses)
            if tables is None:
                continue
            all_fold.extend(tables["fold"])
            all_clin.extend(tables["clinical"])
            all_se.extend(tables["sample_eff"])
            all_lobo.extend(tables["lobo"])
            mk = fm.stem.split("_")[0].lower()
            fm_oof[mk] = oof_data
        except Exception as e:
            log.exception("Failed on %s: %s", fm.name, e)

    # 3 (across FMs). DeLong pairwise on LR OOF preds — align via Cell_ID
    delong_rows = []
    if "delong" in analyses:
        log.info("[5/5] DeLong pairwise (across FMs, Cell_ID-aligned)")
    if "delong" in analyses and fm_oof:
        # Reload each FM's adata briefly to get the Cell_IDs corresponding to
        # the valid (non-healthy/unknown) cells used for OOF predictions.
        fm_cellids = {}
        for fm_path in fm_files:
            mk = fm_path.stem.split("_")[0].lower()
            if mk not in fm_oof:
                continue
            try:
                a = sc.read_h5ad(fm_path)
                if "Cell_ID" not in a.obs.columns:
                    log.warning("  %s: no Cell_ID column, skipping for DeLong", mk)
                    continue
                # Align: rebuild valid mask exactly the same way as in per_fm
                y_full, valid = build_mortality(a.obs["covid_severity"].astype(str).values)
                cellids = a.obs["Cell_ID"].astype(str).values[valid]
                fm_cellids[mk] = cellids
                log.info("  %s: %d Cell_IDs for OOF cells", mk, len(cellids))
                del a
            except Exception as e:
                log.warning("  %s: %s", mk, e)

        # Compute intersection of Cell_IDs across all FMs
        if len(fm_cellids) >= 2:
            common_cellids = set(next(iter(fm_cellids.values())))
            for cids in fm_cellids.values():
                common_cellids &= set(cids)
            log.info("  Common Cell_IDs across %d FMs: %d", len(fm_cellids), len(common_cellids))
            common_list = sorted(common_cellids)
            common_set = set(common_list)
            # Build aligned (y_aligned, proba_aligned[mk]) over the common Cell_IDs
            # First need a reference y per common Cell_ID. Take from first FM.
            first_mk = next(iter(fm_cellids))
            first_cellids = fm_cellids[first_mk]
            first_y = fm_oof[first_mk]["y"]
            first_proba = fm_oof[first_mk]["lr_oof_proba"]
            cid_to_pos = {c: i for i, c in enumerate(first_cellids)}
            keep_pos_first = [cid_to_pos[c] for c in common_list if c in cid_to_pos]
            y_aligned = first_y[keep_pos_first]

            aligned = {}
            for mk in fm_cellids:
                cellids = fm_cellids[mk]
                proba = fm_oof[mk]["lr_oof_proba"]
                cid_to_pos_mk = {c: i for i, c in enumerate(cellids)}
                keep_mk = [cid_to_pos_mk[c] for c in common_list if c in cid_to_pos_mk]
                aligned[mk] = proba[keep_mk]
                # also re-derive y to ensure consistency (sanity check)
                y_mk = fm_oof[mk]["y"][keep_mk]
                if not np.array_equal(y_mk, y_aligned):
                    log.warning("  %s: y mismatch on common Cell_IDs (severity disagreement?)", mk)
            log.info("  DeLong on %d FMs with %d aligned cells", len(aligned), len(y_aligned))
            delong_rows = delong_pairwise(y_aligned, aligned)
        else:
            log.warning("  Fewer than 2 FMs with Cell_ID — skipping DeLong")

    # --- Write outputs (skip CSVs for analyses that produced no rows, so a
    # subset run cannot clobber full-run outputs with empty files) ---
    df_fold = pd.DataFrame(all_fold)
    df_clin = pd.DataFrame(all_clin)
    df_se = pd.DataFrame(all_se)
    df_lobo = pd.DataFrame(all_lobo)
    df_delong = pd.DataFrame(delong_rows)
    df_summary = None
    if not df_fold.empty:
        df_summary = (df_fold.groupby(["model", "classifier"])
                      .agg(auc_mean=("auc", "mean"), auc_std=("auc", "std"),
                           balanced_accuracy_mean=("balanced_accuracy", "mean"))
                      .reset_index())

    written = []
    if df_summary is not None:
        df_fold.to_csv(args.output_dir / "mortality_per_fold.csv", index=False)
        df_summary.to_csv(args.output_dir / "mortality_summary.csv", index=False)
        written += ["mortality_per_fold.csv", "mortality_summary.csv"]
    if not df_clin.empty:
        df_clin.to_csv(args.output_dir / "clinical_metrics.csv", index=False)
        written.append("clinical_metrics.csv")
    if not df_se.empty:
        df_se.to_csv(args.output_dir / "sample_efficiency.csv", index=False)
        written.append("sample_efficiency.csv")
    if not df_lobo.empty:
        df_lobo.to_csv(args.output_dir / "lobo.csv", index=False)
        written.append("lobo.csv")
    if not df_delong.empty:
        df_delong.to_csv(args.output_dir / "delong_pairwise.csv", index=False)
        written.append("delong_pairwise.csv")

    log.info("=" * 60)
    log.info("Wrote %d CSVs to %s: %s", len(written), args.output_dir, written)
    if df_summary is not None:
        print("\n=== MORTALITY SUMMARY ===")
        print(df_summary.to_string(index=False))
    if not df_clin.empty:
        print("\n=== CLINICAL METRICS ===")
        print(df_clin.to_string(index=False))
    if not df_delong.empty:
        print("\n=== DELONG ===")
        print(df_delong.to_string(index=False))


if __name__ == "__main__":
    main()
