# Figure and Table Inventory

## Repository presentation figures

The portable-demo README uses a separate [figure record](assets/readme/README.md):
`workflow.svg` illustrates the current software, and `historical_umap.png`
redraws the committed scPRINT-2 coordinate table with cohort and severity colors.
The latter is an explicitly historical display, not a new post-rebuild result.
Its renderer, source hash, drawing settings and label counts are recorded in
`assets/readme/provenance.json`. The existing synthetic metric preview remains
in `examples/benchmark/`. These figures do not replace the manuscript assets below.

## Historical manuscript assets

> **Refresh required after the 2026-09 rebuild.** This inventory describes the
> pre-rebuild manuscript assets. Input-dependent figures and tables are
> provisional until P1–P3 in `docs/execution-plan.md` complete.

Every final figure/table used in the manuscript/report, with its source data,
generating script, parameters, and manifest. Figures are gitignored locally
(`*.png`, `*.pdf`, `figures/`) and uploaded to Overleaf separately; this inventory
is the durable provenance record.

Status values: `planned` · `needs-evidence` · `verified` · `inconsistent` · `removed`

Rules:
- No figure is manually edited without documenting the edit here.
- Prefer regenerating figures from code over editing image files.
- Every table has machine-readable source data.
- Every manuscript figure/table appears in this inventory.

| ID | File | Panel/Table | Purpose | Source Data | Generating Script | Params/Config | Manifest | Used In | Status |
|---|---|---|---|---|---|---|---|---|---|
| Fig1 | `reports/manuscript/figures/Figure1_overview.png` (gitignored) | overview | Study-design schematic (5-FM platelet-OOD benchmark) | — (schematic; no data) | `scripts/26_figures/code/generate_fig1_overview_5fm.py` → copied by `regenerate_indexed_figures.py` | `configs/hpcc.yaml` | `scripts/26_figures/results/manifest.json` | Manuscript Fig 1 (`fig:overview`, main.tex L146) | **verified** | Diagram, no quantitative claim; generator in-repo. |
| Fig2 | `reports/manuscript/figures/Figure2_performance.png` (gitignored) | A–B | Per-FM mortality LR/RF AUC bars + PCA-200+LR | `scripts/26_figures/data/mortality_lr_rf.csv` (+`pca_baseline.csv`) | `scripts/26_figures/code/regenerate_indexed_figures.py` ✅ | `configs/hpcc.yaml` | `scripts/26_figures/results/manifest.json` | Manuscript Fig 2 (`fig:performance`, L343) | **verified** | Source `mortality_lr_rf.csv` is **byte-identical** to the verified C004 `…/covid_sepsis_only/mortality_lr_rf.csv`. |
| Fig3 | `reports/manuscript/figures/Figure3_cross_disease.png` (gitignored) | a–b | Bidirectional cross-disease transfer + activation axis | `scripts/26_figures/data/cross_disease_severe_excluded.csv` (+`activation_summary_per_fm.csv`) | `scripts/26_figures/code/regenerate_indexed_figures.py` ✅ | `configs/hpcc.yaml` | `scripts/26_figures/results/manifest.json` | Manuscript Fig 3 (`fig:cross_disease`, L486) | **verified** | Snapshot **regenerated 2026-06-08** by `scripts/26_figures/code/build_cross_disease_snapshot.py` to the C002 manuscript set (STATE/UCE/scPRINT-2/Geneformer from `…_20260504_152836.csv` + scGPT from `…_20260430_193915.csv`; STACK/TF dropped). Panel-a LogReg AUCs match manuscript: STATE 0.798/0.707, UCE 0.807/0.646, scPRINT-2 0.855/0.543, scGPT 0.352/0.555. The generator's previously **hardcoded scGPT** values (0.352/0.555) were removed — scGPT is now read from the snapshot (identical output, fully data-traced). |
| Fig4 | `reports/manuscript/figures/Figure4_lobo.png` (gitignored) | heatmap | Leave-one-cohort-out per-cohort AUC heatmap | `scripts/26_figures/data/lobo.csv` | `scripts/26_figures/code/regenerate_indexed_figures.py` ✅ | `configs/hpcc.yaml` | `scripts/26_figures/results/manifest.json` | Manuscript Fig 4 (`fig:lobo`, L511) | **verified** | Source `lobo.csv` is **byte-identical** to the verified C004 `…/covid_sepsis_only/lobo.csv`. |
| Fig5 | `reports/manuscript/figures/Figure5_adaptation.pdf` (tracked vector) | curve | Cross-disease adaptation curve | `scripts/26_figures/data/adaptation_curve_summary.csv` | `scripts/26_figures/code/regenerate_indexed_figures.py` ✅ | `configs/hpcc.yaml` | `scripts/26_figures/results/manifest.json` | Manuscript Fig 5 (`fig:adaptation_curve`, L544) | **verified** | Source is **byte-identical** to the verified C005 `…/adaptation_curve_bootstrap/summary.csv` (PCA asymptote 0.964). |
| TabPerf | `scripts/24_option_c_canonical_input/results/option_c_primary_benchmark/covid_sepsis_only_mortality.csv` | `tab:performance` | primary mortality benchmark (COVID+Sepsis subset, 42,095/322) | `data/processed/covid_sepsis_mortality.h5ad` | `scripts/24_option_c_canonical_input/code/` | `configs/hpcc.yaml` | `scripts/24_option_c_canonical_input/results/manifest.json` | Manuscript | **verified** | All AUCs match to 3 dp incl. per-fold strings (see claim C004). **Re-pinned 2026-06-06:** prior pin `mortality_summary.csv` is a non-disease-matched run that main.tex L317 excludes from headlines. tab:delong / tab:calibration source = `…/covid_sepsis_only/{delong_pairwise,calibration}.csv`. |

> **Rows corrected 2026-06-08.** The previous Fig1/Fig2/Fig4 rows were **mis-seeded** — they
> described figures (model-ranking, diff-diagnosis confusion, sample-eff/dropout) that are **not
> the manuscript's figures**. The manuscript's actual figures are `Figure1_overview` …
> `Figure5_adaptation` (`\includegraphics` at main.tex L146/343/486/511/544), **all produced by
> the single generator `scripts/26_figures/code/regenerate_indexed_figures.py`** reading committed
> snapshots in `scripts/26_figures/data/*.csv`. Verified that **Fig2/Fig4/Fig5 source CSVs are
> byte-identical to the verified C004/C005 claim evidence**; Fig1 is a data-free schematic.
> **Fig3 is `needs-rerun`** — its `cross_disease_severe_excluded.csv` snapshot carries the wrong
> FM set (STACK/TF instead of Geneformer/scGPT) vs the C002 manuscript set.
>
> The generator's hardcoded `ROOT` (old `~/Desktop/302006` path) was fixed to resolve relative to
> the script (survives the relocation) + `OUT.mkdir`, so figures are now regenerable in-place.
> Figure files are gitignored (`*.png`; vector `*.pdf`/`*.svg` tracked) — "verify a figure" =
> regenerate from the ✅ generator against the in-repo `26_figures/data/` snapshot and diff vs the
> Overleaf version. **All 5 manuscript figures are now `verified`** (Fig3 snapshot regenerated from
> C002 evidence 2026-06-08; scGPT de-hardcoded). The figure data snapshots in `26_figures/data/`
> are byte-identical to / derived from the verified C002/C004/C005 claim evidence.
