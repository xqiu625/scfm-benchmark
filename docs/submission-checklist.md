# Submission Checklist — BIB v1

> **Scope note (2026-09-21):** this is now a submission-logistics and historical
> completion checklist. The active scientific plan is `docs/execution-plan.md`;
> the QC rebuild makes earlier compile/readiness statements provisional.

**Target:** Briefings in Bioinformatics
**Last updated:** 2026-08-05

Living punch list. Items move into the "done" section as they land. Replaces
the archived `TODO_before_publication.md` (severity-task framing).

---

## Done in the 2026-08-05 submission pass

- [x] **Title changed** — now "Stress-Testing Single-Cell Foundation Models: An
  Out-of-Distribution Platelet Benchmark for In-Hospital Mortality Prediction"
  (`main.tex` + `supplementary.tex` titles updated).
- [x] **scGPT table gaps — 4 of 7 filled from existing verified CSVs** (no HPCC
  rerun needed): `tab:lobo_per_cohort` + `tab:per_cohort_auc` (from
  `covid_sepsis_only/{lobo,per_cohort_oof}.csv`), `tab:ablation` + `tab:fusion`
  (from stage-30 `results/{ablation,fusion}_summary_fullpool.csv`, 2026-06-09
  batch; shared keys bit-identical). Also fixed a wrong per-cohort bold
  (scPRINT-2 0.622 > STATE 0.620 on `302005data02` LOBO) and moved scGPT's
  fusion row into sorted position.
- [x] **"Queued / pending / will be reported" framing removed** from all
  captions and Results/Clinical-metrics text; replaced with factual
  provenance statements.
- [x] **Stale-number fixes** — Limitations §2 (322 donors, STATE folds
  0.549–0.849); Limitations §4 (data-overlap caveat correctly attributed to the
  *legacy* continued-pre-training result, not the corrected LoRA analysis);
  abstract states the 42,095/322 cohort is a subset of the 47,977/413 pool;
  contribution 3 + statistical-replicability bullet harmonized (322 headline /
  413 full pool); Kendall-τ cross-ref → `ssec:mortality_benchmark`;
  per-cohort range text updated (data02 0.26–0.45, data16 0.63–0.76);
  Fig 4 caption "three of seven" → "five of seven"; ablation text updated for
  scGPT (PCA-robustness co-leader, largest z-score gain).
- [x] **All figure snapshots refreshed + figures regenerated locally**
  (`reports/manuscript/figures/`, 600-dpi PNG + PDF + SVG): 5 snapshots were
  stale/missing scGPT (`mortality_lr_rf`, `lobo`, `per_cohort_oof`,
  `delong_pairwise`, `calibration`, `activation_summary_per_fm`) — now
  byte-identical to the verified stage-24/stage-30 evidence. Fixed Figure 2
  (scGPT was "data pending"; bars now sorted by LR AUC to match the caption),
  Figure 3 (added Geneformer to panel a, scGPT to panel b, "only
  bidirectional" → "strongest", dropped "(legacy)", legend/x-limits fixed),
  Figure 1 (title overlap fixed, scGPT "legacy reference" label removed,
  claims softened to match the text).
- [x] **HPCC rerun package for the 3 remaining scGPT tables** — SLURM array
  `scripts/30_review_response_stats/code/submit_scgpt_auxiliary_reruns.sh`
  (validated; 5 underlying scripts patched with `--models`/`--analyses` glue).
  Only `--array=0-2` is still required (adaptation curve, sample efficiency,
  gene dropout); tasks 3–4 (ablation, fusion) are optional — data already
  filled from the stage-30 fullpool CSVs. Details:
  `scripts/30_review_response_stats/SCGPT_AUXILIARY_RERUNS.md`.
- [x] **`environments/` created** (README + `scfm_benchmark.yml`); Methods and
  Data-Availability wording aligned (previously claimed manifests that did not
  exist). README "Reproducing the manuscript" table rewritten to match the
  current BIB manuscript (was stale LMRL-workshop mapping; headline tables now
  point to stage 27).
- [x] **BIB author-guideline checks** (via archived OUP pages): no numeric
  abstract cap ("short" — 240 words fine); Key Points required, 3–5 sentences
  (our 5 compliant); word limit 2,000–7,000 by article type (our main text
  ≈4,500 incl. captions — fine); format-free initial submission allowed.
- [x] **Vancouver numbered citations** — BIB uses numbered `[1]` style, not
  author-year: `natbib` switched to `[numbers,sort&compress]` and the 76-entry
  bibliography reordered by first citation (verified 76/76 bijection).
- [x] **5 remaining references verified**: `corpas2025diversity`,
  `pot2019gender` (Somatechnics 9:170–187 — entry already correct),
  `mcdermott2021reproducibility`, `fahsbender2025czi` (arXiv 2507.10502v2
  confirmed via arXiv API). `pearce2025transcriptformer` updated to the
  journal version: *Science* 393(6806):aec8514 (2026), Crossref-verified.

(Pre-2026-08-05 done items: author block, funding, LoRA chain, SHAP mapping,
terminology, docs site, CLI, docs cleanup, abstract/Key Points/CRediT/COI/
Ethics/bios, references 16→76 with corrections, literature positioning —
see git history of this file.)

---

## Done in the 2026-08-05 HPCC top-up pass (post-rerun)

- [x] **scGPT top-up reruns completed on HPCC** (array 0–4, all exit 0;
  4.5–46 min, peak 8.2 GB) — remaining 3 tables filled from the new CSVs:
  `tab:adaptation_curve` (scGPT Y-int 0.874 / n=100 0.882 / n=1000 0.896 /
  asymptote 0.906 / slope +0.008 — flattest curve), `tab:sample_eff`
  (0.541/0.487/0.560/0.654/0.672, +0.130), `tab:dropout`
  (0.795→0.769, −3.3%).
- [x] **Dropout experiment reworded to "embedding-dimension dropout"**
  everywhere (Methods, Results, abstract, Key Points, Discussion, Conclusion,
  Clinical Implications, supp caption, Figure 1). The script zeroes post-hoc
  embedding dims — it measures representational redundancy, NOT missing-gene
  tolerance; Methods now says so explicitly and the clinical-deployment claim
  is downgraded. True input-gene masking + re-inference = future work.
- [x] **scGPT 0.874-vs-0.352 reconciled** — protocol-dependence paragraph in
  `ssec:adaptation_curve` (severe-excluded/full-sepsis bootstrap vs
  severe-included/5-donor-holdout; UCE moves oppositely, so the protocols
  reorder models); "avoid scGPT for cross-disease deployment" softened to
  protocol-dependent guidance; abstract, Key Points, Discussion, Conclusion,
  and cross-disease Results all carry the qualifier.
- [x] **Sample-efficiency claim corrected** — UCE is lowest at every fraction
  except 10% (scGPT 0.487 < UCE 0.489); scGPT low-end non-monotonicity noted.
- [x] **True primary-5-FM fusion refit adopted** (`option_c_fusion_primary5/`):
  `tab:fusion` + Results/Discussion prose switched to the refit (single STATE
  0.694; best fusion late stacking 0.684; early-concat lowest at 0.671). The
  stage-30 6-FM-pool values are superseded.
- [x] **Figure 5 regenerated as a 2-panel figure** matching its caption —
  `make_figure5_adaptation` now includes scGPT in both panels (curves + the
  three-quantities vs PCA-reference bars) and the per-method asymptote lines;
  caption updated (scGPT added, band = 95% bootstrap CI, symlog axis).
- [x] **Environment exports landed** (2026-08-05, HEAD 022802d): `state_env.yml`
  + `scpeft_env.yml` (conda prefixes), `scgpt_finetune.yml` (conda named),
  `scprint2_env-requirements.txt` (uv venv — no `conda env export`; exported via
  pip). `environments/README.md` rewritten with per-environment commands.
  Methods + Data Availability sentences updated to name the manifests.

## Remaining local manuscript work (2026-08-05 pass)

(none — all local items closed)

## Remaining HPCC work

- [x] **E-tier evidence sync + application (2026-08-06, commits 35131cf +
  02b245a)** — all four rerun output sets committed to the evidence tree; the
  four tables applied and verified cell-by-cell (62/62 numeric cells match the
  CSVs across the three numeric tables). Includes the GSEA "---" FDR cells (Focal
  adhesion 0.005/rank 8; all others 0.005) and the three rounding corrections
  (STATE DeepMLP 0.698, scGPT 10% dropout 0.773, Late average Δ −1.2%).
- [x] **Overleaf compile clean (2026-08-06)** — 0 errors, 0 unresolved refs,
  0 overfull, no LaTeX warnings; main 43 pp + supplementary 14 pp. Compile
  order that resolves the xr cross-refs: supplementary → main → supplementary
  ×2 → main ×2.
- [ ] **Zenodo archive（投稿前完成更稳）** — 归档内容:5 个 FM 的 Option C
  预计算 embedding(`02_EMBEDDINGS/option_c/*_covid_canonical.h5ad`,HPCC 上)。
  步骤:zenodo.org 用 GitHub 登录 → New upload → 传 h5ad 文件 → "Reserve DOI"
  先拿 DOI 再发布。拿到 DOI 后把 `main.tex` Data and Code Availability 里的
  "will be archived to Zenodo with a permanent DOI upon acceptance" 改成引用真实
  DOI(期刊对 data availability 的审查在评审阶段就会看)。

## Code audit fixes landed (2026-08-06)

- [x] **Multiclassifier torch-seeded rerun applied (tab:multiclassifier)** —
  MLP columns replaced with the seeded values; both wrong claims corrected:
  "no MLP beats LR for ANY FM (deficits −0.003 to −0.025)" replaces the
  nonexistent "+0.003 STATE", and the RF statement now carries the UCE
  exception (RF 0.636 > LR 0.631). UCE row bold moved to RF. FigureS2 result
  box fixed to "DeepMLP does not beat LR for any FM (seeded rerun)".
- [x] **GSEA full-ranking rewrite (E-tier rerun, job 27260144)** — the 200-gene
  truncation was upstream (`biological_interpretation_option_c.py` TOP_N_GENES),
  not in the GSEA script. Full-ranking (~13,066 genes/dim, 303 KEGG pathways)
  results are now in `tab:per_fm_gsea` + rewritten §per_fm_gsea: STATE
  (cytoskeleton, NES −2.14, FDR 0.005) and Geneformer (translation, +2.52,
  FDR 0.005) corroborated; scPRINT-2 now mixed (pathway-level Platelet
  activation −3.06 rank 1, gene-level ribosomal); **UCE has no significant
  enrichment (FDR 0.661)** — its axis assignment now rests on gene-level
  SHAP/Jaccard evidence, and the spurious "UCE Platelet activation NES −2.20,
  FDR <0.005" claim is deleted everywhere. The Geneformer dim-297 Ca²⁺
  subdominant-dim observation was dropped (unverifiable under full ranking).
  Abstract/Discussion/SHAP paragraphs re-hung on gene-level evidence; the
  truncated variant is acknowledged in caption + text as superseded.
  Optional: sync `option_c_gsea_fullrank/summary.csv` to fill the remaining
  FDR cells currently shown as "---" (scGPT/scPRINT-2 rows).
- [x] **Dropout mortality-endpoint rerun applied (tab:dropout)** — new numbers
  (STATE 0.792→0.782, highest at every level; UCE the only non-degrading FM,
  +0.003 at 70%). Caption now truthfully says mortality endpoint; main text,
  Clinical Implications item 5, Conclusion updated. Drop column now absolute
  AUC change.
- [x] **Fusion nested-CV stacking applied (tab:fusion)** — new row
  `late_stacking_lr_nested` 0.693 sits within 0.001 of single STATE 0.694
  (tie); legacy leaky row retained at 0.684 for transparency (+0.009 delta
  reported in caption + text). Results + Discussion softened from
  "underperform" to "fail to outperform / tie".
- [x] **A-tier number corrections** — `tab:crossdisease` CIs now match the C002
  evidence CSVs (scGPT's already did); narrative corrected (STATE reverse CI
  excludes only scPRINT-2 by pairwise overlap); Figure 3a and table now agree.
- [x] **B-tier Methods disclosures** — LoRA per-fold best-epoch selected on the
  test fold (inflation direction favors LoRA; +0.113 is an upper bound);
  Geneformer's 47,114-cell input is the scPRINT-2-pipeline QC-filtered object;
  `tab:activation` module-score AUROC non-fatal set includes healthy controls.
- [x] **C-tier code hygiene (9 items)** — submit_27_array.sh off-by-one (scGPT
  silently dropped in array mode); scgpt added to 3 rosters; unsafe obsm
  fallback removed in 5 scripts (loud error now); shap NameError branch fixed;
  05/07 default models = primary 5-FM; CLI embed/evaluate/diagnose repaired
  (loud failures instead of silent no-ops); LICENSE added (MIT);
  fig8/broken runners archived under `archive/superseded_figures_2026-08-06/`;
  stale figure docstrings + 06/29 provenance pins; scanpy added to [dev].

## Remaining local / user tasks

- [x] **FigureS1–S3 workflow schematics generated** (2026-08-05). Correction to
  the previous note: these did **not** exist on the HPCC working copy. An
  exhaustive `find` over `/bigdata/godziklab/shared/Xinru/302006` returned no
  `fig_*_workflow/` directories, no schematic PDFs, and no `FigureS*` files
  anywhere on the cluster — only the three generator scripts are tracked. The
  generators are pure matplotlib (no data, no GPU), so they were run locally
  and the `fig_*_workflow/` → `FigureS1/S2/S3` copy step of
  `regenerate_indexed_figures.py` was applied on its own, leaving the
  visually-verified Figure 1 and Figure 5 untouched. All 10 figures referenced
  by `main.tex` + `supplementary.tex` now resolve in
  `reports/manuscript/figures/`.
- [x] **FigureS1–S3 text-overlap defects fixed** (2026-08-06). All listed
  collisions resolved in the three generators and verified by re-render: S1
  (box text fractions, LoRA strip shortened clear of "Evaluate", matrix label
  arrow repointed off the equation); S2 (box text fractions, "Parameter
  status"/legend shifted below the frozen-matrix box, matrix→box arrow moved
  off the `X_fm` label, both inter-row arrows re-pointed onto box edges so
  they no longer cross `p=0.3`); S3 (box text fractions, "Layer N" labels
  drawn above the chips, chip spacing widened so "Dense" fits, head/train box
  stack re-spaced). PNG+PDF+SVG for all three are in
  `reports/manuscript/figures/` under both the `fig_*_workflow/` stems and the
  `FigureS1/S2/S3_*` names.
- [ ] **Upload 10 figures to Overleaf** — `Figure1_overview.png`,
  `Figure2_performance.png`, `Figure3_cross_disease.png`, `Figure4_lobo.png`,
  `Figure5_adaptation.pdf`, `FigureS1_lora_frozen_workflow.pdf`,
  `FigureS2_deepmlp_workflow.pdf`, `FigureS3_geneformer_lora_workflow.pdf`,
  `FigureS4_per_cohort_oof.pdf`, `FigureS5_delong_significance.pdf` (all in
  `reports/manuscript/figures/`, which is gitignored — they transfer by upload,
  not through git). Figure 1 and Figure 5 changed again this pass, so
  re-upload those two even if uploaded earlier.
- [ ] **Author biographies** — review the drafts (written from public
  affiliation only).
- [ ] **ORCID iD** — BIB 投稿系统要求 submitting author 提供 ORCID iD;没有的话
  先去 orcid.org 注册(免费,几分钟),然后在投稿门户里关联。通讯作者邮箱
  xqiu014@ucr.edu 建议与 ORCID 记录一致。
- [ ] **Cover letter + AI/LLM-use disclosure** — cover letter 要点:一句话卖点
  (首个真正 OOD 底物 + 临床终点的 scFM stress test)、与现有评测文献的差异
  (patient-level CV、统一输入、强基线、双向跨疾病、阴性结果如实报告)、
  声明无利益冲突、数据/代码公开。BIB 要求:若写作或分析中使用了 AI/LLM 工具,
  需在 cover letter(及 Methods/Acknowledgements)中披露使用情况。
- [ ] **Overleaf compile** — `main.tex` + `supplementary.tex`; note the `xr`
  cross-document order (supplementary → main ×2). Citations are now numbered;
  expect a different look from the author-year draft.
- [ ] **BIB portal submission** — PDFs, cover letter, declarations.

---

## Deferred to v1 revision (out of submission scope)

- Cosine fairness rerun on the primary 5-FM set — current Discussion text quotes the result qualitatively; full table can be deferred.
- Quantization on the primary 5-FM set — not in v1 critical path.
- True input-gene-masking dropout evaluation (mask raw genes, re-run each FM encoder; GPU) — the published embedding-dimension dropout measures redundancy only; Methods says so explicitly.
