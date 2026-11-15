# scfm-benchmark

> HISTORICAL SNAPSHOT, retained 2026-11-13. The following README predates the
> portable demo and contains superseded performance, portability and submission
> claims. Later edits remove material maintained in a separate repository.
> Do not use it as current evidence. Start with the repository README.

**Stress-testing single-cell foundation models across disease and species.**

A reproducible benchmark + diagnostic framework for clinical deployment of
single-cell foundation models (scFMs). Evaluates 6 modern FMs — STATE, UCE,
TranscriptFormer, Geneformer, scPRINT-2, STACK — on COVID-19 mortality
(47,808 human platelets, 413 patients).

The main finding is uncomfortable: a PCA-200 + logistic regression baseline
**beats every FM at full cross-disease data** in the adaptation curve (0.964
asymptote) and beats 2 of 6 FMs on patient-level within-disease mortality.
Only STATE is bidirectionally robust across diseases under the corrected
severe=survived label definition (0.790 forward, 0.641 reverse); STACK
collapses to anti-correlated forward (0.123) under the same label, a
concrete label-boundary failure of in-context retrieval.

## Pipeline

```
       ┌─────────────────────────────────────────────────────────────┐
raw  → │  embed  →  evaluate  →  diagnose  →  report                │
h5ad   │  (FM)     (patient-CV)  (suite)      (HTML bundle)         │
       └─────────────────────────────────────────────────────────────┘
                              │
                              ├── adaptation curve  (sample-efficiency
                              │                       vs from-scratch)
                              ├── GSEA per-dim      (biology axes)
                              ├── SHAP attribution  (gene-level)
                              └── Jaccard overlap   (cross-FM agreement)
```

Same canonical raw-counts input ("Option C") goes to every FM; each applies
its own native preprocessing. Patient-level `StratifiedGroupKFold` everywhere
— cell-level CV inflates AUCs by 7–17 points and is the dominant source of
overstated performance in published scFM benchmarks.

## Quickstart

```bash
pip install scfm-benchmark

scfm-bench embed     --h5ad raw.h5ad        --model state    --out emb.h5ad
scfm-bench evaluate  --emb  emb.h5ad        --task mortality --out results/
scfm-bench diagnose  --emb  emb.h5ad        --out diagnostic/
scfm-bench report    --results-dir results/ --out report.html
```

Embedding generation requires the FM-specific environment and weights on a
GPU host; `scfm-bench embed` dispatches to the appropriate SLURM submitter
under `scripts/foundation_models/<MODEL>/`.

## Headline findings

| Result | Number | Implication |
|---|---|---|
| STATE binary AUC (raw LR, patient-CV) | **0.825** | Only FM beating PCA-200 baseline |
| PCA-200 baseline | 0.814 | Beats UCE, TF, Geneformer, scGPT |
| STATE LOBO mean AUC | 0.789 | Best distribution-shift generalisation |
| UCE gene-dropout drop @ 70% | **0.0%** | Most dropout-robust FM |
| STACK cross-disease COVID→Sepsis (severe=survived) | **0.123** | Anti-correlated; label-boundary failure of in-context retrieval |
| STATE bidirectional cross-disease | 0.790 / 0.641 | Only bidirectional generalizer (severe=survived label) |
| PCA-200+LR adaptation-curve asymptote | **0.964** | Beats every FM at full data |
| Late-fusion AUC | 0.849 | Beats any single FM |
| FP16 quantisation | 0% AUC loss, 2× speedup | Production-ready inference |

Full numbers and confidence intervals: `reports/manuscript/main.tex` and
`docs/current-status.md`. See `docs/submission-checklist.md` for the
current submission punch list.

## Repository layout

Follows the Universal Data Science Pipeline Architecture
(`UNIVERSAL_DATA_SCIENCE_PIPELINE_ARCHITECTURE.md`). See `docs/pipeline-index.md`
for the full stage table and `docs/decisions.md` for the restructure rationale.

```
.
├── configs/                # Environment-specific paths/params (defaults/hpcc/local/cloud.yaml)
├── src/scfm_benchmark/      # Pip-installable CLI + API; io/preprocessing/.../submission subpkgs
├── scripts/                # Numbered workflow stages, one folder per ACC number
│   ├── 00_setup/           # HPCC sync, cell-index inspection, ACC17–21 figures
│   ├── 02_embedding_classification/
│   ├── 06_batch_generalization/      # Cross-disease transfer
│   ├── 17_domain_adapted_fm/         # Continued pre-training / LoRA
│   ├── 24_option_c_canonical_input/  # Standardised input + primary mortality benchmark
│   ├── 29_adaptation_curve/          # From-scratch baseline curves
│   ├── 99_audit/                     # Traceability/reproducibility audit (no new science)
│   ├── foundation_models/            # Per-FM install + embedding submitters
│   ├── finetuning/ tracking/ agent/                 # Reusable runners
│   └── README.md · MANUSCRIPT_MAPPING.md
├── data/                   # raw/ interim/ processed/ external/ (large files on HPCC; gitignored)
├── results/                # Cross-stage aggregated outputs (features/models/predictions/…)
├── figures/                # Working figures (gitignored; Overleaf is the upload target)
├── reports/
│   ├── manuscript/         # LaTeX source for Briefings in Bioinformatics
│   └── audit/              # code_review / data_integrity / reproduction / claim_check / …
├── docs/                   # pipeline_index, data_contracts, claim_traceability,
│   │                       # figure_table_inventory, experiment_log, audit_report, decisions
│   ├── runbooks/           # gate1_gate2, hpcc
│   └── superseded_2026-09-21/  # Archived memory-bank — history, not status
├── notebooks/              # Exploration only (promote stable logic into src/)
├── tests/                  # pytest suite
├── logs/                   # Run logs (HPCC is canonical)
├── archive/ · legacy/      # Superseded experiments, retained for reproducibility
├── environment.yml · requirements.txt · pyproject.toml · Makefile
```

## Reproducing the manuscript

### Quick path — single FM via the CLI

```bash
pip install scfm-benchmark
scfm-bench embed    --h5ad raw.h5ad       --model state    --out emb_state.h5ad
scfm-bench evaluate --emb  emb_state.h5ad --task mortality --out results_state/
scfm-bench diagnose --emb  emb_state.h5ad --out diagnostic_state/
```

Reproduces the headline benchmark + the three-view Diagnostic Suite for one FM
in one call. The `embed` step dispatches to the per-FM SLURM submitter on HPCC.

### Full path — every manuscript table and figure

**Prerequisite — generate Option C standardised embeddings** (one job per FM, ~24 GPU-hours total):

```bash
sbatch scripts/24_option_c_canonical_input/code/step2c_human_covid_canonical_input.sh
```

Outputs land in `02_EMBEDDINGS/option_c/<fm>_covid_canonical.h5ad`.

**Then run the per-analysis SLURM submitters:**

| Manuscript outputs | Submitter / script |
|---|---|
| Headline cohort (COVID + Sepsis subset) — `tab:performance`, `tab:baseline`, `tab:calibration`, `tab:delong`, `tab:batch`, `tab:lobo_per_cohort`, `tab:per_cohort_auc` (+ per-FM OOF predictions) | `scripts/27_covid_sepsis_only_rerun/code/submit_27_array.sh` then `submit_27_gather.sh` |
| `tab:seed_stability`, `tab:severe`, `tab:subspace`, donor-clustered DeLong CIs | `scripts/30_review_response_stats/code/` |
| `tab:multiclassifier` (DeepMLP) | `scripts/24_option_c_canonical_input/code/run_multiclassifier_option_c.sh` |
| `tab:ablation` | `.../run_ablation_option_c.sh` |
| `tab:fusion` | `.../run_fusion_option_c.sh` |
| `tab:embedding_quality` (scIB) | `.../run_embedding_quality_option_c.sh` |
| `tab:activation` | `.../run_activation_markers.sh` |
| `tab:per_fm_gsea` | `.../run_gsea_per_fm_dim.sh` |
| `tab:per_fm_shap`, `tab:per_fm_jaccard` | `.../run_shap_attribution.sh` + `run_shap_gene_map.sh` |
| `tab:clinical` | `.../run_calibration_per_cohort.sh` |
| `tab:sample_eff` | `scripts/05_sample_efficiency/code/submit_sample_efficiency.sh` |
| `tab:dropout` | `scripts/07_gene_dropout/code/submit_gene_dropout.sh` |
| `tab:crossdisease`, `fig:cross_disease` | `scripts/06_batch_generalization/code/cross_disease_fatal_vs_survived.py` |
| `tab:adaptation_curve`, `fig:adaptation_curve` | `scripts/29_adaptation_curve/code/run_adaptation_curve_bootstrap.sh` |
| `tab:domainadapt` (LoRA fine-tune) | `scripts/17_domain_adapted_fm/code/` (chain Jobs 23853347–23866507) |
| Fig 1–5, Fig S4–S5 | `scripts/26_figures/code/regenerate_indexed_figures.py` (+ `generate_fig1_overview_5fm.py`, `regenerate_supp_figures_5fm.py`; read committed snapshots in `scripts/26_figures/data/`) |
| Fig S1–S3 (workflow schematics) | hand-drawn Inkscape sources at `reports/manuscript/figures/fig_*_workflow/*.svg` |

### Cross-validation discipline

All benchmarks use 5-fold patient-level `StratifiedGroupKFold` with `Sample_ID` as the group column (322 donors in the disease-matched COVID + Sepsis headline cohort; 413 in the full 5-disease pool). Cell-level CV is deliberately not exposed — it inflates AUCs 7–17 percentage points and is the dominant source of overstated performance in published scFM benchmarks.

### Conda environments

The core downstream analysis environment is pinned in `environments/scfm_benchmark.yml` (repo-root `environment.yml`; a Dockerfile is also provided). Per-FM HPCC environments (`state_env`, `scprint2_env`, `scpeft_env`, …) are documented in `environments/README.md` and `configs/hpcc.yaml`.

## License

MIT — see `LICENSE`.
