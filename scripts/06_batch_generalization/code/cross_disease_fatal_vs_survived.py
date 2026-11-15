#!/usr/bin/env python3
"""Cross-Disease Generalization — fatal vs survived (true outcome label).

Trains on COVID with binary label (fatal=1 vs survived=0) and tests on Sepsis
with the analogous fatal/survived label, and vice versa. Replaces the buggy
`severe_classes = ['severe', 'fatal']` lumping in the original
`cross_disease_generalization.py`, which conflated severe-survivors with fatal
cases and made the cross-disease comparison meaningless.

Usage (current canonical invocation, see run_cross_disease_fatal.sh):
    python cross_disease_fatal_vs_survived.py \
        --models STATE UCE TranscriptFormer scPRINT2 STACK Geneformer \
        --covid-label-col covid_severity \
        --covid-fatal-tokens fatal \
        --covid-survived-tokens mild,moderate,recovered,severe \
        --sepsis-label-col covid_severity \
        --sepsis-fatal-tokens fatal \
        --sepsis-survived-tokens severe \
        --group-col Sample_ID \
        --n-bootstrap 200

Notes:
- COVID 'severe' is by default INCLUDED in the survived set (severe = ICU
  survivor, matching sepsis-severe): the default --covid-survived-tokens is
  mild,moderate,recovered,severe. Pass --exclude-covid-severe for the legacy
  pre-2026-05-01 behavior (severe excluded as ambiguous). NOTE: the
  manuscript's tab:crossdisease predates this default — see the provenance
  pin in run_cross_disease_fatal.sh.
- Per project policy: no synthetic data fallback. Errors out if labels missing.
- Includes patient-level bootstrap CI (resample over Sample_ID groups).
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("cross_disease_fatal")

BASE = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB = BASE / "02_EMBEDDINGS"
RESULTS = BASE / "05_RESULTS" / "cross_disease_fatal_vs_survived"

# Option C canonical-input embeddings (2026-05-01). Each FM was given the
# identical raw input AnnData and applied its own native preprocessing.
MODEL_PATHS = {
    "STATE": EMB / "option_c/state_covid_canonical.h5ad",
    "UCE": EMB / "option_c/uce_covid_canonical.h5ad",
    "TranscriptFormer": EMB / "option_c/transcriptformer_covid_canonical.h5ad",
    "scPRINT2": EMB / "option_c/scprint2_covid_canonical.h5ad",
    "STACK": EMB / "option_c/stack_covid_canonical.h5ad",
    "Geneformer": EMB / "option_c/geneformer_covid_canonical.h5ad",
    # Legacy paths (kept commented for cross-comparison):
    # "scGPT": EMB / "scgpt/covid_scgpt_clean.h5ad",  # dropped from main pipeline
}
OBSM_KEYS = {
    "STATE": ("X_state", "X_emb", "X_STATE"),
    "UCE": ("X_uce", "X_uce_4layer"),
    "TranscriptFormer": ("X_transcriptformer", "X_TranscriptFormer"),
    "scPRINT2": ("X_scprint2", "scprint_emb"),
    "STACK": ("X_stack",),
    # Legacy:
    "scGPT": ("X_scGPT", "X_scgpt"),
    "Geneformer": ("X_geneformer", "X_Geneformer"),
}
DISEASE_ALIASES = {
    "COVID": ("COVID-19", "covid-19", "COVID", "covid", "COVID19"),
    "Sepsis": ("Sepsis", "sepsis", "SEPSIS"),
}


def load_embeddings(model_name: str):
    p = MODEL_PATHS[model_name]
    if not p.exists():
        raise FileNotFoundError(p)
    log.info("Loading %s from %s", model_name, p)
    adata = sc.read_h5ad(p)
    X = None
    used_key = None
    for key in OBSM_KEYS.get(model_name, ()):
        if key in adata.obsm:
            X = adata.obsm[key]
            used_key = key
            break
    if X is None:
        raise ValueError(
            f"No embedding key for {model_name}. Tried {OBSM_KEYS[model_name]}; "
            f"available: {list(adata.obsm.keys())}"
        )
    if hasattr(X, "toarray"):
        X = X.toarray()
    log.info("  obsm key=%s  shape=%s", used_key, X.shape)
    return np.asarray(X), adata


def get_disease_array(adata) -> np.ndarray:
    obs = adata.obs
    for col in ("Disease", "disease"):
        if col in obs.columns:
            return obs[col].astype(str).values
    raise ValueError("No Disease/disease column on adata.obs")


def disease_mask(disease: np.ndarray, name: str) -> np.ndarray:
    aliases = DISEASE_ALIASES.get(name, (name,))
    return np.isin(disease, aliases)


def build_binary(
    label_values: np.ndarray,
    fatal_tokens: list[str],
    survived_tokens: list[str] | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (binary_labels, valid_mask).

    binary[i] = 1 if value matches fatal_tokens (case-insensitive)
              = 0 if value matches survived_tokens
    valid[i]  = True if value matched either set; False otherwise (excluded)
    """
    fatal_set = {t.strip().lower() for t in fatal_tokens if t.strip()}
    survived_set = (
        {t.strip().lower() for t in survived_tokens if t.strip()}
        if survived_tokens
        else None
    )

    n = len(label_values)
    binary = np.zeros(n, dtype=np.int8)
    valid = np.zeros(n, dtype=bool)
    for i, raw in enumerate(label_values):
        s = str(raw).strip().lower()
        if s in fatal_set:
            binary[i] = 1
            valid[i] = True
        elif survived_set is None:
            # No explicit survived tokens — anything not unknown/nan/empty counts as 0
            if s not in ("unknown", "nan", "none", "", "na"):
                binary[i] = 0
                valid[i] = True
        elif s in survived_set:
            binary[i] = 0
            valid[i] = True
    return binary, valid


def patient_bootstrap_auc(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    groups: np.ndarray,
    n_iter: int = 200,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Bootstrap AUC by resampling patients (groups) with replacement."""
    rng = np.random.default_rng(seed)
    unique_groups = np.unique(groups)
    aucs = []
    for _ in range(n_iter):
        sampled = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        idx = np.concatenate([np.where(groups == g)[0] for g in sampled])
        if len(np.unique(y_true[idx])) < 2:
            continue
        try:
            aucs.append(roc_auc_score(y_true[idx], y_prob[idx]))
        except ValueError:
            continue
    if not aucs:
        return float("nan"), float("nan"), float("nan")
    arr = np.asarray(aucs)
    return float(arr.mean()), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def transfer(
    X: np.ndarray,
    y: np.ndarray,
    valid: np.ndarray,
    train_mask: np.ndarray,
    test_mask: np.ndarray,
    test_groups: np.ndarray,
    n_bootstrap: int = 200,
) -> dict:
    train_sel = train_mask & valid
    test_sel = test_mask & valid

    Xtr = X[train_sel]
    ytr = y[train_sel]
    Xte = X[test_sel]
    yte = y[test_sel]
    gte = test_groups[test_sel]

    log.info(
        "  Train n=%d (fatal=%d, survived=%d) | Test n=%d (fatal=%d, survived=%d)",
        len(ytr), int(ytr.sum()), int((1 - ytr).sum()),
        len(yte), int(yte.sum()), int((1 - yte).sum()),
    )

    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        log.warning("  Skipping: insufficient class diversity")
        return {}

    scaler = StandardScaler()
    Xtr_s = scaler.fit_transform(Xtr)
    Xte_s = scaler.transform(Xte)

    out = {"n_train": int(len(ytr)), "n_test": int(len(yte)),
           "fatal_train": int(ytr.sum()), "fatal_test": int(yte.sum())}

    # Logistic regression
    clf_lr = LogisticRegression(max_iter=2000, random_state=42, n_jobs=-1, class_weight="balanced")
    clf_lr.fit(Xtr_s, ytr)
    p_lr = clf_lr.predict_proba(Xte_s)[:, 1]
    auc_lr = float(roc_auc_score(yte, p_lr))
    bal_lr = float(balanced_accuracy_score(yte, clf_lr.predict(Xte_s)))
    mean_lr, lo_lr, hi_lr = patient_bootstrap_auc(yte, p_lr, gte, n_iter=n_bootstrap)
    out["LogReg"] = {"auc": auc_lr, "balanced_accuracy": bal_lr,
                     "boot_auc_mean": mean_lr, "boot_ci_lo": lo_lr, "boot_ci_hi": hi_lr}

    # Random forest
    clf_rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1, class_weight="balanced")
    clf_rf.fit(Xtr_s, ytr)
    p_rf = clf_rf.predict_proba(Xte_s)[:, 1]
    auc_rf = float(roc_auc_score(yte, p_rf))
    bal_rf = float(balanced_accuracy_score(yte, clf_rf.predict(Xte_s)))
    mean_rf, lo_rf, hi_rf = patient_bootstrap_auc(yte, p_rf, gte, n_iter=n_bootstrap)
    out["RandomForest"] = {"auc": auc_rf, "balanced_accuracy": bal_rf,
                           "boot_auc_mean": mean_rf, "boot_ci_lo": lo_rf, "boot_ci_hi": hi_rf}

    return out


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", nargs="+", default=list(MODEL_PATHS.keys()))
    ap.add_argument("--covid-label-col", default="covid_severity")
    ap.add_argument("--covid-fatal-tokens", default="fatal",
                    help="Comma-separated, case-insensitive")
    ap.add_argument("--covid-survived-tokens",
                    default="mild,moderate,recovered,severe",
                    help="Default INCLUDES 'severe' (COVID-severe = ICU survivor, "
                         "matches sepsis-severe=ICU survivor). 'healthy' and 'unknown' "
                         "are EXCLUDED from mortality prediction.")
    ap.add_argument("--exclude-covid-severe", action="store_true",
                    help="Legacy mode: exclude COVID 'severe' as ambiguous "
                         "(pre-2026-05-01 default behavior).")
    ap.add_argument("--sepsis-label-col", required=True,
                    help="obs column that contains fatal/survived for sepsis cells")
    ap.add_argument("--sepsis-fatal-tokens", default="fatal,died,deceased,non-survivor")
    ap.add_argument("--sepsis-survived-tokens", default="survived,survivor,alive,discharged")
    ap.add_argument("--group-col", default="Sample_ID",
                    help="obs column for patient-level bootstrap grouping")
    ap.add_argument("--n-bootstrap", type=int, default=200)
    ap.add_argument("--directions", nargs="+",
                    default=["COVID->Sepsis", "Sepsis->COVID"],
                    choices=["COVID->Sepsis", "Sepsis->COVID"])
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    covid_fatal = args.covid_fatal_tokens.split(",")
    covid_surv = args.covid_survived_tokens.split(",")
    if args.exclude_covid_severe:
        covid_surv = [t for t in covid_surv if t.strip().lower() != "severe"]
        log.info("Legacy mode: excluding COVID 'severe' from survived set")

    sep_fatal = args.sepsis_fatal_tokens.split(",")
    sep_surv = args.sepsis_survived_tokens.split(",")

    rows = []
    full = {}

    for model in args.models:
        log.info("=" * 70)
        log.info("Model: %s", model)
        log.info("=" * 70)
        try:
            X, adata = load_embeddings(model)
        except Exception as e:
            log.error("Skip %s: %s", model, e)
            continue

        disease = get_disease_array(adata)
        cov_mask = disease_mask(disease, "COVID")
        sep_mask = disease_mask(disease, "Sepsis")
        log.info("COVID cells: %d  Sepsis cells: %d", int(cov_mask.sum()), int(sep_mask.sum()))

        if args.covid_label_col not in adata.obs.columns:
            raise ValueError(f"--covid-label-col '{args.covid_label_col}' missing from obs")
        if args.sepsis_label_col not in adata.obs.columns:
            raise ValueError(f"--sepsis-label-col '{args.sepsis_label_col}' missing from obs")
        if args.group_col not in adata.obs.columns:
            raise ValueError(f"--group-col '{args.group_col}' missing from obs")

        cov_lab = adata.obs[args.covid_label_col].astype(str).values
        sep_lab = adata.obs[args.sepsis_label_col].astype(str).values
        groups = adata.obs[args.group_col].astype(str).values

        cov_y, cov_valid = build_binary(cov_lab, covid_fatal, covid_surv)
        sep_y, sep_valid = build_binary(sep_lab, sep_fatal, sep_surv)

        # Combine: a cell uses COVID label if it's COVID, sepsis label if sepsis.
        y = np.zeros(len(disease), dtype=np.int8)
        valid = np.zeros(len(disease), dtype=bool)
        y[cov_mask] = cov_y[cov_mask]
        valid[cov_mask] = cov_valid[cov_mask]
        y[sep_mask] = sep_y[sep_mask]
        valid[sep_mask] = sep_valid[sep_mask]

        log.info("Valid+labeled cells: COVID=%d  Sepsis=%d",
                 int((cov_mask & valid).sum()), int((sep_mask & valid).sum()))

        full[model] = {}
        for direction in args.directions:
            train, test = direction.split("->")
            tm = cov_mask if train == "COVID" else sep_mask
            te = sep_mask if test == "Sepsis" else cov_mask
            log.info("\nDirection: %s", direction)
            res = transfer(X, y, valid, tm, te, groups, args.n_bootstrap)
            full[model][direction] = res
            for clf, m in res.items():
                if clf in ("n_train", "n_test", "fatal_train", "fatal_test"):
                    continue
                rows.append({
                    "model": model,
                    "direction": direction,
                    "classifier": clf,
                    "auc": m["auc"],
                    "balanced_accuracy": m["balanced_accuracy"],
                    "boot_auc_mean": m["boot_auc_mean"],
                    "boot_ci_lo": m["boot_ci_lo"],
                    "boot_ci_hi": m["boot_ci_hi"],
                    "n_train": res["n_train"],
                    "n_test": res["n_test"],
                    "fatal_train": res["fatal_train"],
                    "fatal_test": res["fatal_test"],
                })
                log.info("  %s: AUC=%.4f BalAcc=%.4f BootMean=%.4f [%.4f, %.4f]",
                         clf, m["auc"], m["balanced_accuracy"],
                         m["boot_auc_mean"], m["boot_ci_lo"], m["boot_ci_hi"])

    if not rows:
        log.warning("No results generated")
        return 1

    df = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / f"cross_disease_fatal_vs_survived_{timestamp}.csv"
    df.to_csv(csv_path, index=False)
    json_path = RESULTS / f"cross_disease_fatal_vs_survived_{timestamp}.json"
    json_path.write_text(json.dumps(full, indent=2))
    log.info("Saved: %s", csv_path)
    log.info("Saved: %s", json_path)

    log.info("\n" + "=" * 70)
    log.info("SUMMARY (LogReg)")
    log.info("=" * 70)
    lr = df[df["classifier"] == "LogReg"]
    for direction in args.directions:
        sub = lr[lr["direction"] == direction]
        log.info("\n%s:", direction)
        for _, r in sub.iterrows():
            log.info("  %-18s AUC=%.4f  BootCI=[%.4f, %.4f]",
                     r["model"], r["auc"], r["boot_ci_lo"], r["boot_ci_hi"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
