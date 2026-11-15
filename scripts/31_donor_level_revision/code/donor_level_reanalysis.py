#!/usr/bin/env python3
"""
Revision P1-B / P1-C: donor-level mortality evaluation + cell-count association.

Adam 2026-08-31: "All the training is on the cell level ... then averaging over
patients to make a person-level prediction" and "You could almost predict
outcome from cell counts alone". The headline metric is moved to the DONOR
level: each donor's out-of-fold (OOF) cell predictions are aggregated to a
single score, and ROC-AUC is computed across donors (n = 322, the correct
unit of inference). The cell-weighted OOF AUC is retained as a secondary,
representation-level metric.

Inputs (default: ../../30_review_response_stats/data/, canonical OOF
predictions from the stage-27 COVID+Sepsis patient-level CV rerun):
  oof_<FM>.csv    cell_id, sample_id, severity, y, p_lr, p_rf

Outputs (../results/):
  donor_level_auc.csv         per-FM donor-level AUC under mean/median/trimmed-mean
                              aggregation + donor bootstrap 95% CI (primary table)
  donor_level_pairwise.csv    paired dAUC on common donors, bootstrap CI, BH-FDR q
  cellcount_association.csv   Spearman(n_cells, donor score) per FM +
                              AUC within cell-count tertiles (P1-C)
  donor_predictions_<FM>.csv  per-donor aggregated scores (for downstream joins)

Usage:
  python donor_level_reanalysis.py [--data DIR] [--out DIR] [--B 5000]
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from repeated_cv import (DEFAULT_SEEDS, aggregate_patient_level, paired_delta_auc,
                          repeated_patient_auc, summarize_repeated)

FMS = ["STATE", "UCE", "Geneformer", "scPRINT-2", "scGPT"]
SCORES = ["p_lr", "p_rf"]
RNG = np.random.default_rng(42)


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, p)


def boot_auc_ci(y, p, B, rng):
    """Percentile bootstrap CI, resampling donors (rows) with replacement."""
    n = len(y)
    aucs = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, n)
        aucs[b] = safe_auc(y[idx], p[idx])
    aucs = aucs[~np.isnan(aucs)]
    return np.percentile(aucs, [2.5, 97.5])


def aggregate_donors(df, score):
    """One row per donor: aggregated score, label, cell count."""
    def trim(x):
        return stats.trim_mean(x, 0.1) if len(x) > 2 else x.mean()

    g = df.groupby("sample_id")
    out = pd.DataFrame({
        "y": g["y"].max(),               # fatal if any fatal cell (verified below)
        "y_nunique": g["y"].nunique(),
        "n_cells": g.size(),
        "mean": g[score].mean(),
        "median": g[score].median(),
        "trim": g[score].apply(trim),
    })
    n_bad = int((out["y_nunique"] > 1).sum())
    if n_bad:
        print(f"  WARNING: {n_bad} donors have mixed labels; using max(y)")
    return out.drop(columns="y_nunique")


def run_repeated(args, fms):
    """CV Protocol v2 (2026-09-04): repeated-CV reporting.

    Expects oof_<fm>_seed<NN>.csv (NN = 42..42+n_seeds-1) in --data, each with
    a patient_id column already attached (the CV-fitting scripts -- e.g.
    rerun_covid_sepsis_mortality_patient_v2.py --n-seeds N -- write this
    directly; no --patient-map merge needed here).

    Outputs (--out):
      donor_level_auc_repeated.csv      per-fm: repeat-mean +/- sd AUC (LR/RF,
                                         mean/median patient aggregation)
      donor_level_pairwise_repeated.csv same-seed-paired dAUC: mean + 95% CI
                                         (t-based) + Wilcoxon signed-rank p
      per_seed_auc_repeated.csv         the raw per-seed AUC matrix (for audit)
    """
    seeds = DEFAULT_SEEDS[:args.n_seeds]
    oof_by_seed_by_fm = {}
    y_by_fm = {}
    pid_by_fm = {}

    for fm in fms:
        oof_by_seed = {}
        for seed in seeds:
            path = os.path.join(args.data, f"oof_{fm}_seed{seed}.csv")
            df = pd.read_csv(path)
            assert "patient_id" in df.columns, f"{fm} seed{seed}: no patient_id column"
            oof_by_seed[seed] = {"p_lr": df["p_lr"].to_numpy(), "p_rf": df["p_rf"].to_numpy()}
            if seed == seeds[0]:
                y_by_fm[fm] = df["y"].to_numpy()
                pid_by_fm[fm] = df["patient_id"].to_numpy()
        oof_by_seed_by_fm[fm] = oof_by_seed
        print(f"{fm}: loaded {len(seeds)} seed files, "
              f"{len(np.unique(pid_by_fm[fm]))} patients, "
              f"{int(pd.DataFrame({'p': pid_by_fm[fm], 'y': y_by_fm[fm]}).groupby('p')['y'].max().sum())} "
              f"fatal patients", flush=True)

    # ------------------------------------- headline: repeat-mean +/- sd AUC
    auc_rows = []
    per_seed_rows = []
    auc_cache = {}  # (fm, agg, classifier) -> (aucs array, seeds)
    for fm in fms:
        for agg in ["mean", "median"]:
            for score_key, clf_name in [("p_lr", "LR"), ("p_rf", "RF")]:
                aucs, s = repeated_patient_auc(oof_by_seed_by_fm[fm], pid_by_fm[fm],
                                                y_by_fm[fm], score_key=score_key, agg=agg)
                auc_cache[(fm, agg, clf_name)] = aucs
                summ = summarize_repeated(aucs)
                n_pat = len(np.unique(pid_by_fm[fm]))
                n_fatal_pat = int(pd.DataFrame({'p': pid_by_fm[fm], 'y': y_by_fm[fm]})
                                   .groupby('p')['y'].max().sum())
                auc_rows.append({"fm": fm, "aggregation": agg, "classifier": clf_name,
                                  "n_patients": n_pat, "n_fatal_patients": n_fatal_pat,
                                  "n_seeds": summ["n"], "auc_repeat_mean": round(summ["mean"], 4),
                                  "auc_repeat_sd": round(summ["sd"], 4)})
                for seed, auc in zip(s, aucs):
                    per_seed_rows.append({"fm": fm, "aggregation": agg, "classifier": clf_name,
                                           "seed": seed, "auc": round(auc, 4)})
    tab = pd.DataFrame(auc_rows)
    tab.to_csv(os.path.join(args.out, "donor_level_auc_repeated.csv"), index=False)
    pd.DataFrame(per_seed_rows).to_csv(os.path.join(args.out, "per_seed_auc_repeated.csv"),
                                        index=False)
    print("\n=== CV Protocol v2: repeat-mean +/- sd AUC (LR, mean aggregation) ===")
    print(tab.query("aggregation == 'mean' and classifier == 'LR'")
          [["fm", "n_patients", "n_fatal_patients", "n_seeds", "auc_repeat_mean",
            "auc_repeat_sd"]].to_string(index=False))

    # seed-42 single-split view, clearly labeled secondary (per reporting spec item 3)
    print("\n=== secondary: seed-42 single-split AUC (LR, mean aggregation) ===")
    for fm in fms:
        aucs = auc_cache[(fm, "mean", "LR")]
        idx42 = seeds.index(42) if 42 in seeds else 0
        print(f"  {fm}: {aucs[idx42]:.4f} (seed {seeds[idx42]})")

    # ------------------------------------- pairwise: same-seed-paired dAUC
    pair_rows = []
    for i in range(len(fms)):
        for j in range(i + 1, len(fms)):
            a, b = fms[i], fms[j]
            aucs_a = auc_cache[(a, "mean", "LR")]
            aucs_b = auc_cache[(b, "mean", "LR")]
            res = paired_delta_auc(aucs_a, aucs_b)
            pair_rows.append({"pair": f"{a} - {b}", "n_seeds": res["n_seeds"],
                               "mean_delta": round(res["mean_delta"], 4),
                               "ci_lo": round(res["ci_lo"], 4), "ci_hi": round(res["ci_hi"], 4),
                               "wilcoxon_p": round(res["wilcoxon_p"], 5)
                               if not np.isnan(res["wilcoxon_p"]) else float("nan")})
    pw = pd.DataFrame(pair_rows)
    pw["q_bh"] = stats.false_discovery_control(pw["wilcoxon_p"].to_numpy())
    pw.to_csv(os.path.join(args.out, "donor_level_pairwise_repeated.csv"), index=False)
    print("\n=== CV Protocol v2: same-seed-paired dAUC (LR, mean aggregation) ===")
    print(pw.to_string(index=False))
    print(f"\nSAVED -> {args.out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--data", default=os.path.join(
        here, "..", "..", "30_review_response_stats", "data"))
    ap.add_argument("--out", default=os.path.join(here, "..", "results"))
    ap.add_argument("--B", type=int, default=5000)
    ap.add_argument("--group-col", default="sample_id",
                    help="column to aggregate/group by; use patient_id after the "
                         "2026-09-02 multi-draw finding (requires --patient-map)")
    ap.add_argument("--patient-map", default=os.path.join(
        here, "..", "results", "donor_meta_with_patient_id.csv"),
                    help="CSV with sample_id,patient_id (from the revision finding)")
    ap.add_argument("--extra-methods", default="",
                    help="comma-separated extra method names with their own "
                         "oof_<name>.csv in --data (e.g. PCA-200, for the Job H2 "
                         "within-COVID table); appended after the 5 primary FMs, "
                         "does not affect default behavior when empty")
    ap.add_argument("--exclude-methods", default="",
                    help="comma-separated method names to drop from the default "
                         "5 (e.g. scGPT, while its embedding is still pending)")
    ap.add_argument("--n-seeds", type=int, default=1,
                    help="CV Protocol v2 (2026-09-04): read oof_<fm>_seed<NN>.csv "
                         "(NN in 42..42+N-1) instead of a single oof_<fm>.csv, and "
                         "report per the repeated-CV spec (repeat-mean +/- sd "
                         "headline, same-seed-paired dAUC). Default 1 = legacy "
                         "single-split behavior, unchanged.")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    excluded = {m.strip() for m in args.exclude_methods.split(",") if m.strip()}
    fms = [m for m in FMS if m not in excluded] + \
        [m.strip() for m in args.extra_methods.split(",") if m.strip()]

    if args.n_seeds > 1:
        run_repeated(args, fms)
        return

    donor = {}
    for fm in fms:
        path = os.path.join(args.data, f"oof_{fm}.csv")
        df = pd.read_csv(path)
        if args.group_col != "sample_id":
            pmap = pd.read_csv(args.patient_map)[["sample_id", args.group_col]]
            df = df.merge(pmap, on="sample_id", how="left")
            n_missing = int(df[args.group_col].isna().sum())
            assert n_missing == 0, f"{fm}: {n_missing} cells lack {args.group_col}"
            df["sample_id"] = df[args.group_col]
        d = aggregate_donors(df, "p_lr")
        d_rf = aggregate_donors(df, "p_rf")[["mean"]].rename(columns={"mean": "mean_rf"})
        d = d.join(d_rf)
        donor[fm] = d
        d.to_csv(os.path.join(args.out, f"donor_predictions_{fm}.csv"))
        print(f"{fm}: {len(df)} cells -> {len(d)} donors "
              f"({int(d['y'].sum())} fatal)", flush=True)

    # ------------------------------------------------ P1-B: donor-level AUC
    rows = []
    for fm in fms:
        d = donor[fm]
        y = d["y"].to_numpy()
        # secondary: cell-weighted OOF AUC on the same file (sanity vs stage 30)
        cell = pd.read_csv(os.path.join(args.data, f"oof_{fm}.csv"))
        auc_cell = safe_auc(cell["y"].to_numpy(), cell["p_lr"].to_numpy())
        for agg in ["mean", "median", "trim"]:
            p = d[agg].to_numpy()
            lo, hi = boot_auc_ci(y, p, args.B, RNG)
            rows.append({"fm": fm, "aggregation": agg, "classifier": "LR",
                         "n_donors": len(d), "n_fatal_donors": int(y.sum()),
                         "auc_donor": round(safe_auc(y, p), 4),
                         "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                         "auc_cell_weighted_ref": round(auc_cell, 4)})
        p_rf = d["mean_rf"].to_numpy()
        lo, hi = boot_auc_ci(y, p_rf, args.B, RNG)
        rows.append({"fm": fm, "aggregation": "mean", "classifier": "RF",
                     "n_donors": len(d), "n_fatal_donors": int(y.sum()),
                     "auc_donor": round(safe_auc(y, p_rf), 4),
                     "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                     "auc_cell_weighted_ref": round(auc_cell, 4)})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(args.out, "donor_level_auc.csv"), index=False)
    print("\n=== P1-B: donor-level AUC (LR, mean aggregation) ===")
    print(tab.query("aggregation == 'mean' and classifier == 'LR'")
            [["fm", "n_donors", "n_fatal_donors", "auc_donor", "ci_lo", "ci_hi",
              "auc_cell_weighted_ref"]].to_string(index=False))

    # --------------------------------------- P1-B: paired pairwise dAUC (donor level)
    pairs = []
    for i in range(len(fms)):
        for j in range(i + 1, len(fms)):
            a, b = fms[i], fms[j]
            common = donor[a].join(donor[b], how="inner", lsuffix="_a", rsuffix="_b")
            y = common["y_a"].to_numpy()
            pa = common["mean_a"].to_numpy()
            pb = common["mean_b"].to_numpy()
            d_auc = safe_auc(y, pa) - safe_auc(y, pb)
            n = len(y)
            boots = np.empty(args.B)
            for k in range(args.B):
                idx = RNG.integers(0, n, n)
                boots[k] = safe_auc(y[idx], pa[idx]) - safe_auc(y[idx], pb[idx])
            boots = boots[~np.isnan(boots)]
            lo, hi = np.percentile(boots, [2.5, 97.5])
            p_two = 2 * min((boots <= 0).mean(), (boots >= 0).mean())
            pairs.append({"pair": f"{a} - {b}", "n_common_donors": n,
                          "d_auc": round(d_auc, 4),
                          "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                          "p_boot": round(p_two, 5)})
    pw = pd.DataFrame(pairs)
    pw["q_bh"] = stats.false_discovery_control(pw["p_boot"].to_numpy())
    pw.to_csv(os.path.join(args.out, "donor_level_pairwise.csv"), index=False)
    print("\n=== P1-B: paired donor-level dAUC ===")
    print(pw.to_string(index=False))

    # --------------------------------------- P1-C: cell-count association
    rows = []
    for fm in fms:
        d = donor[fm]
        rho, pval = stats.spearmanr(d["n_cells"], d["mean"])
        d = d.copy()
        d["tertile"] = pd.qcut(d["n_cells"], 3, labels=["low", "mid", "high"])
        rec = {"fm": fm, "spearman_rho": round(rho, 4),
               "spearman_p": f"{pval:.2e}"}
        for t in ["low", "mid", "high"]:
            sub = d[d["tertile"] == t]
            rec[f"auc_{t}_count"] = round(
                safe_auc(sub["y"].to_numpy(), sub["mean"].to_numpy()), 4)
            rec[f"n_{t}"] = len(sub)
            rec[f"n_fatal_{t}"] = int(sub["y"].sum())
        rows.append(rec)
    cc = pd.DataFrame(rows)
    cc.to_csv(os.path.join(args.out, "cellcount_association.csv"), index=False)
    print("\n=== P1-C: cell-count association ===")
    print(cc.to_string(index=False))
    print(f"\nSAVED -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
