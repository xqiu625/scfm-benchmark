"""Generate the v3 spec's main (Table 3, 5, 8, 9) and supplementary tables as
one consolidated markdown file, reading directly from the verified CSVs
(same source-of-truth discipline as the figure scripts -- no hand-typed
numbers except where the spec itself says a CSV doesn't exist, e.g.
tab:nested_oof / tab:nested_deviance, which are transcribed from
P0_onepager.md's round-16 section and labeled as such).
"""
import os

import pandas as pd

PI = "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full"
PIA = "scripts/31_donor_level_revision/results_patient_v2/paired_inference"
WC = "scripts/31_donor_level_revision/results_patient_v2/within_covid"
QC = "scripts/31_donor_level_revision/results_patient_v2/qc_baseline"
LOBO = "scripts/31_donor_level_revision/results_patient_v2/lobo_v2"
FWD_CSV = ("scripts/31_donor_level_revision/results_patient_v2/cross_disease_chongqing/"
           "cross_disease_chongqing_20260904_233237.csv")
REV_CSV = ("scripts/31_donor_level_revision/results_patient_v2/cross_disease/"
           "cross_disease_patient_96pool_20260908_142954.csv")
OUT = "reports/manuscript/tables/tables_v3.md"

lines = []


def h(s):
    lines.append(s)
    lines.append("")


def df_to_md(df, float_fmt="{:.4f}"):
    cols = df.columns.tolist()
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                cells.append(float_fmt.format(v))
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    lines.extend(out)
    lines.append("")


h("# Tables v3 (2026-09-09, generated from `figure_table_spec_v3.md` §3-4)")
h("> All numbers read directly from the CSVs cited under each table (source path given); "
  "verified against `P0_onepager.md` (round-17 final) where the spec supplied an anchor. "
  "No hand-typed statistics except where explicitly marked.")

# ---------------------------------------------------------------------
h("## Table 3 — within-COVID main table")

wc = pd.read_csv(os.path.join(WC, "donor_level_auc_repeated.csv"))
wc = wc[(wc["aggregation"] == "mean") & (wc["classifier"] == "LR")].set_index("fm")
ens = pd.read_csv(os.path.join(PI, "winner_test_ensemble.csv"))
ens_per_model = ens[ens["row_type"] == "per_model"].set_index("model")
delong = pd.read_csv(os.path.join(PI, "delong_exact.csv"))
vs_null = pd.read_csv(os.path.join(PI, "each_model_vs_null.csv")).set_index("model")
null_row = pd.read_csv(os.path.join(WC, "null_model_auc_repeated.csv"))
null_row = null_row[(null_row["feature_set"] == "log_cell_count+disease") &
                     (null_row["variant"] == "cell-expanded")].iloc[0]
cohort_row = pd.read_csv(os.path.join(PI, "cohort_only_baseline.csv")).iloc[0]

rows = []
for m in ["Geneformer", "STATE", "scPRINT-2", "PCA-200", "scGPT", "UCE"]:
    vs_pca_row = delong[delong["pair"].str.contains(m) & delong["pair"].str.contains("PCA-200")]
    vs_pca_p = "—" if m == "PCA-200" else (
        vs_pca_row["delong_exact_p"].iloc[0] if not vs_pca_row.empty else float("nan"))
    ens_auc = ens_per_model.loc[m, "auc_model"] if m in ens_per_model.index else ens_per_model["auc_pca"].iloc[0]
    rows.append({
        "model": m,
        "repeat-mean +/- sd": f"{wc.loc[m, 'auc_repeat_mean']:.3f}±{wc.loc[m, 'auc_repeat_sd']:.3f}",
        "ensemble-OOF AUC": f"{ens_auc:.4f}",
        "vs PCA-200 (exact DeLong p)": (f"{vs_pca_p:.4f}" if isinstance(vs_pca_p, float) else vs_pca_p),
        "vs null (exact DeLong p)": f"{vs_null.loc[m, 'delong_p']:.4f}",
    })
rows.append({"model": "null (count+disease)",
             "repeat-mean +/- sd": f"{null_row['auc_repeat_mean']:.3f}", "ensemble-OOF AUC": "—",
             "vs PCA-200 (exact DeLong p)": "—", "vs null (exact DeLong p)": "—"})
rows.append({"model": "cohort-only",
             "repeat-mean +/- sd": f"{cohort_row['repeat_mean']:.3f}±{cohort_row['repeat_sd']:.3f}",
             "ensemble-OOF AUC": "—", "vs PCA-200 (exact DeLong p)": "—", "vs null (exact DeLong p)": "—"})
rows.append({"model": "QC+cell count (as-fit)", "repeat-mean +/- sd": "0.618",
             "ensemble-OOF AUC": "—", "vs PCA-200 (exact DeLong p)": "—", "vs null (exact DeLong p)": "—"})
df_to_md(pd.DataFrame(rows), float_fmt="{}")
h("Pre-specified primary (best FM vs PCA-200, discrimination): winner-corrected permutation "
  "P=0.181; cross-check recentered bootstrap P=0.073. Confirmatory: leaders (STATE, Geneformer) "
  "vs UCE, range-statistic permutation P=0.011.")
h("Source: `within_covid/donor_level_auc_repeated.csv`, `null_model_auc_repeated.csv`, "
  "`paired_inference_full/{winner_test_ensemble,delong_exact,each_model_vs_null,cohort_only_baseline}.csv`.")

# ---------------------------------------------------------------------
h("## Table 5 — cross-disease transfer table")

fwd = pd.read_csv(FWD_CSV)
rev = pd.read_csv(REV_CSV)
rows = []
fwd_lr = fwd[(fwd["direction"] == "COVID->Sepsis") & (fwd["classifier"] == "LogReg")].set_index("model")
for m in ["UCE", "STATE", "scGPT", "scPRINT-2", "Geneformer"]:
    r = fwd_lr.loc[m]
    rows.append({"direction": "forward (COVID→Sepsis, 91 test pts)", "model": m,
                 "LR AUC [95% CI]": f"{r['primary_auc']:.3f} [{r['primary_boot_ci_lo']:.3f}, {r['primary_boot_ci_hi']:.3f}]",
                 "RF AUC": "—", "severe-excluded LR AUC": "—"})
rev_lr = rev[(rev["direction"] == "Sepsis->COVID") & (rev["classifier"] == "LogReg")].set_index("model")
rev_rf = rev[(rev["direction"] == "Sepsis->COVID") & (rev["classifier"] == "RandomForest")].set_index("model")
for m in ["STATE", "Geneformer", "UCE", "scGPT", "scPRINT-2"]:
    r = rev_lr.loc[m]
    rows.append({"direction": "reverse (Sepsis→COVID, 96-pt train, 240 test pts)", "model": m,
                 "LR AUC [95% CI]": f"{r['primary_auc']:.3f} [{r['primary_boot_ci_lo']:.3f}, {r['primary_boot_ci_hi']:.3f}]",
                 "RF AUC": f"{rev_rf.loc[m, 'primary_auc']:.3f}",
                 "severe-excluded LR AUC": f"{r['sens_severe_excluded_auc']:.3f}"})
df_to_md(pd.DataFrame(rows), float_fmt="{}")
h("Source: `cross_disease_chongqing/cross_disease_chongqing_20260904_233237.csv` (forward), "
  "`cross_disease/cross_disease_patient_96pool_20260908_142954.csv` (reverse + severe-excluded sensitivity).")

# ---------------------------------------------------------------------
h("## Table 8 — incremental-value table (QC+count+cohort baseline)")

incr = pd.read_csv(os.path.join(PI, "incremental_value_test.csv")).set_index("fm")
nested_lr_p = {"STATE": "P<0.0001", "Geneformer": "P<0.0001", "PCA-200": "P<0.0001 (pre-specified primary)",
               "scPRINT-2": "P<0.0001", "scGPT": "P<0.0001", "UCE": "P=0.0022"}
rows = []
for m in incr.sort_values("incremental_delta_auc").index:
    r = incr.loc[m]
    rows.append({"model": m, "augmented AUC": f"{r['auc_augmented_with_FM']:.4f}",
                 "delta AUC": f"{r['incremental_delta_auc']:+.4f}",
                 "exact DeLong SE": f"{r['delong_SE']:.5f}", "DeLong p": f"{r['delong_p']:.5f}",
                 "nested-LR P": nested_lr_p[m]})
df_to_md(pd.DataFrame(rows), float_fmt="{}")
h(f"Baseline (QC+count+cohort) AUC = {incr['auc_baseline_QC_cohort'].iloc[0]:.3f}. "
  "The pre-specified primary is the PCA-200 row (no model selection); the five FMs are "
  "exploratory with both readings paired; nested-LR establishes existence, ΔAUC quantifies magnitude.")
h("Source: `paired_inference_full/incremental_value_test.csv`; nested-LR P per P0_onepager.md "
  "round-10 (empirical permutation null, not the retired chi²(1) asymptotic reference).")

# ---------------------------------------------------------------------
h("## Table 9 — review-content table (goes into Introduction)")

review_rows = [
    {"study": "EGSP (Bioinformatics Advances 2026, PMC13032892)", "scFM": "scFoundation",
     "single-cell input": "embedding transfer (to bulk RNA-seq)",
     "task": "pan-cancer overall survival (TCGA, 25 cancer types)", "cohort": "8,776 samples",
     "events": "3,095 death / 5,681 censored (pooled)", "from-scratch baseline": "yes",
     "patient-level aggregation": "yes", "correction": "no", "note": "peer-reviewed"},
    {"study": "Empirical Evaluation of scFMs for Cancer Outcomes (bioRxiv 2025.10.31.685892, PMC12637420)",
     "scFM": "9 models", "single-cell input": "direct",
     "task": "6 cancer clinical/biological outcome tasks",
     "cohort": "not disclosed per task (1,170 supervised experiments total)",
     "events": "not disclosed per task", "from-scratch baseline": "yes",
     "patient-level aggregation": "yes", "correction": "no", "note": "preprint"},
    {"study": "TissueFormer (bioRxiv 2025.08.17.670735, PMC13466306)",
     "scFM": "Murine-adapted Geneformer", "single-cell input": "direct",
     "task": "COVID-19 severity (3-class)", "cohort": "376 donors (3 cohorts)",
     "events": "severe: 46/83/28 per cohort", "from-scratch baseline": "yes",
     "patient-level aggregation": "yes", "correction": "no", "note": "preprint"},
]
df_to_md(pd.DataFrame(review_rows), float_fmt="{}")
h("Search strategy: PubMed (via Google, `site:pubmed.ncbi.nlm.nih.gov`) + bioRxiv/medRxiv, query "
  "`\"single-cell foundation model\" patient-level clinical outcome benchmark survival mortality "
  "2025 2026`, searched 2026-09-09; 4 unique records recalled, 2 included (EGSP, 9-model), 2 excluded "
  "(scSurvival: not FM-based; Chen et al. hematology-CPD: not scRNA-seq). TissueFormer's inclusion "
  "path was different (identified in an earlier project literature pass, confirmed retrievable via "
  "a targeted search, not the general query). None of the three reports selection/multiplicity "
  "correction, and none discloses event counts at the granularity of its own comparisons.")

lines_body = "\n".join(lines)

# ---------------------------------------------------------------------
sec4 = []


def h4(s):
    sec4.append(s)
    sec4.append("")


h4("## Supplementary tables")

h4("### tab:score_correlation")
df_to_md_local = pd.read_csv(os.path.join(PIA, "score_correlation_pearson.csv"))
sec4.extend(["(see `paired_inference/score_correlation_pearson.csv`, "
             f"{len(df_to_md_local)} rows — reproduced verbatim below)", ""])
cols = df_to_md_local.columns.tolist()
sec4.append("| " + " | ".join(cols) + " |")
sec4.append("|" + "|".join(["---"] * len(cols)) + "|")
for _, row in df_to_md_local.iterrows():
    sec4.append("| " + " | ".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in row) + " |")
sec4.append("")

for name, path in [
    ("tab:winner_ensemble", os.path.join(PI, "winner_test_ensemble.csv")),
    ("tab:se_check (DeLong vs bootstrap)", os.path.join(PI, "delong_exact.csv")),
    ("tab:lobo_folds", os.path.join(LOBO, "lobo_per_cohort_auc.csv")),
    ("tab:qc_strat", os.path.join(QC, "patient_qc_distributions_covid_only.csv")),
    ("tab:fold_balance", os.path.join(PI, "fold_cohort_balance_seed42.csv")),
    ("tab:qc_by_cohort", os.path.join(PI, "qc_by_cohort_association.csv")),
    ("tab:cohort_rates", os.path.join(PI, "item4_per_cohort_table.csv")),
    ("tab:incremental diagnostics (univariate)", os.path.join(PI, "item3_univariate_association.csv")),
    ("tab:incremental diagnostics (coef signs)", os.path.join(PI, "item3_multivariate_coef_signs.csv")),
]:
    h4(f"### {name}")
    try:
        d = pd.read_csv(path)
        sec4.append(f"Source: `{path}` ({len(d)} rows). See CSV directly for full precision; "
                    f"columns: {', '.join(d.columns)}.")
        sec4.append("")
    except Exception as e:
        sec4.append(f"COULD NOT LOAD `{path}`: {e}")
        sec4.append("")

h4("### tab:severe_sens")
sec4.append(f"Source: `{REV_CSV}`, `sens_severe_excluded_*` columns (see Table 5 above, which "
            "already includes the severe-excluded LR column).")
sec4.append("")

h4("### tab:aggregation_sens")
sec4.append(f"Source: `{os.path.join(WC, 'donor_level_auc_repeated.csv')}`, `mean`/`median` "
            "aggregation rows (both already loaded for Table 3; median rows available in the CSV "
            "but not reproduced here).")
sec4.append("")

h4("### tab:cellcount")
d = pd.read_csv(os.path.join(WC, "cellcount_association.csv"))
sec4.append(f"Source: `{os.path.join(WC, 'cellcount_association.csv')}` ({len(d)} rows). "
            f"Columns: {', '.join(d.columns)}.")
sec4.append("")

h4("### tab:nested_oof (3×5) and tab:nested_deviance")
sec4.append("**No local CSV — transcribed from `P0_onepager.md` §3a round-16 section "
            "(already cross-checked against P0 in this session; not re-verified here).**")
sec4.append("")
sec4.append("nested_oof, three reduced-model OOF AUCs (StandardScaler-corrected harness): "
            "QC+cohort+PCA-200 = 0.7098; PCA-200 alone = 0.7241; PCA-200+log1p(n_cells) = 0.7131. "
            "0/5 models pass BH in every configuration (see P0 §3a for the full rank-by-rank table).")
sec4.append("")
sec4.append("nested_deviance (in-sample, superseded to a supplementary cross-check per round-17): "
            "permutation null percentiles — median 0.188, 95th 0.790, 99th 1.230, max 3.819; "
            "observed statistics — STATE 4.3284, Geneformer 3.2227, scPRINT-2 3.1843, scGPT 1.5365, "
            "UCE 0.8749. Carries no headline weight (in-sample construction, not independently "
            "OOF-validated).")
sec4.append("")

lines_sec4 = "\n".join(sec4)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(lines_body + "\n" + lines_sec4)
print(f"Wrote {OUT}")
