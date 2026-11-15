# scGPT auxiliary-evaluation reruns (stage 30)

The scGPT Option C canonical embedding
(`02_EMBEDDINGS/option_c/scgpt_covid_canonical.h5ad`, obsm `X_scgpt`, 47,977 cells)
landed on the HPCC **2026-06-08**, after five auxiliary analyses had already been run
for the other FMs. scGPT therefore appears as `---` in five manuscript tables. This
package recomputes exactly those five missing scGPT evaluations under the identical
protocol (same scripts, same seeds, same paths) as the other FMs.

All jobs are CPU-only sklearn (logistic regression / random forest on frozen
embeddings); no GPU. Environment follows the project convention: PATH-prepend
activation of `302006_env_clean` (`export PATH=$ENV_DIR/bin:$PATH`,
`PYTHONNOUSERSITE=1`), `--partition=batch`, logs in
`/bigdata/godziklab/shared/Xinru/302006/06_LOGS/`.

## How to run

```bash
# all five evaluations as one SLURM array (5 independent tasks):
sbatch scripts/30_review_response_stats/code/submit_scgpt_auxiliary_reruns.sh

# or a single evaluation (task id from the table below):
sbatch --array=2 scripts/30_review_response_stats/code/submit_scgpt_auxiliary_reruns.sh
```

The submitter pre-flights the existence of `scgpt_covid_canonical.h5ad` and exits
loudly if it is missing. Logs: `06_LOGS/0006_scgpt_aux_<jobid>_<task>.{out,err}`.

## The five evaluations

| Task | Evaluation | Manuscript table | Runtime class | Output CSV(s) |
|---|---|---|---|---|
| 0 | Adaptation curve (bootstrap) | `tab:adaptation_curve` (main.tex) | tens of minutes – ~1 h | `05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/adaptation_curve_bootstrap/{scGPT.csv, all_rows.csv, summary.csv}` |
| 1 | Sample efficiency | `tab:sample_eff` (supplementary.tex) | minutes (~10–15 min) | `05_RESULTS/option_c_primary_benchmark/scgpt_aux/sample_efficiency.csv` |
| 2 | Gene dropout | `tab:dropout` (supplementary.tex) | ~1–3 h (heaviest) | `05_RESULTS/covid_severity_benchmark/gene_dropout_robustness_<timestamp>.csv` (+ `.json`) |
| 3 | Embedding ablation | `tab:ablation` (supplementary.tex) | ~1–2 h | `05_RESULTS/option_c_ablation/{per_fm_scgpt.csv, ablation_summary.csv}` |
| 4 | Multimodal fusion refit | `tab:fusion` (supplementary.tex) | ~1 h | `05_RESULTS/option_c_fusion_primary5/{fusion_summary.csv, per_fm_oof.csv}` |

Runtime classes are scaled from the original multi-FM jobs (e.g. Job 23851383 ran
6 FMs in <6 h for dropout; Job 24018607 ran 7 methods for the adaptation curve);
only scGPT (tasks 0–3) or the 5 primary FMs (task 4) are computed here.

### Task 0 — adaptation curve → `tab:adaptation_curve`

```
python3 -u scripts/29_adaptation_curve/code/run_adaptation_curve_bootstrap.py --models scGPT
```

COVID+Sepsis subset; n_sepsis ∈ {0, 10, 50, 100, 500, 1000, all}; 10 seeds;
200-iteration patient-level bootstrap. The script's `EMBEDDINGS` list gained a
`('scGPT', 'scgpt_covid_canonical.h5ad', 'X_scgpt')` entry, plus a `--models`
option. On a subset run the new rows are **merged** into the existing
`all_rows.csv` / `summary.csv` (recomputed methods replaced, others preserved) —
the other FMs' numbers do not change. Table row to fill: `scGPT` Y-intercept
(n=0), n=100, n=1000, asymptote, slope. Deterministic (seed 42 holdout), so the
values equal what a full rerun would produce.

### Task 1 — sample efficiency → `tab:sample_eff`

```
python3 -u scripts/24_option_c_canonical_input/code/primary_benchmark_option_c.py \
    --embedding-dir 02_EMBEDDINGS/option_c \
    --output-dir 05_RESULTS/option_c_primary_benchmark/scgpt_aux \
    --models scgpt --analyses sample_efficiency
```

**Caveat — which script backs this table.** `tab:sample_eff` is produced by
`primary_benchmark_option_c.py` (analysis 4, `sample_efficiency()`), **not** by
`scripts/05_sample_efficiency/code/subsampling_robustness.py`. Verified: the
manuscript values match `05_RESULTS/option_c_primary_benchmark/sample_efficiency.csv`
to 3 decimals (STATE 100 % = 0.712, Geneformer = 0.699, scPRINT-2 = 0.688,
UCE = 0.631; "mean across available 5-fold patient-level CV folds"). The
`05_sample_efficiency` script computes a *different* contrast (non-severe vs
severe on GroupShuffleSplit repeats) and its repo CSV already contains legacy
scGPT rows from the pre-Option-C embedding — do not use it for this table.
`--models` / `--analyses` flags were added to `primary_benchmark_option_c.py`;
a subset run writes CSVs only for analyses that produced rows, into a separate
`scgpt_aux/` dir so the canonical outputs are untouched. To fill the table,
append the `scgpt` rows to `05_RESULTS/option_c_primary_benchmark/sample_efficiency.csv`
(values at 5/10/20/50/100 %, mean across folds where the fold fit).

### Task 2 — gene dropout → `tab:dropout`

```
python3 -u scripts/07_gene_dropout/code/gene_dropout_robustness.py --tasks binary --models scGPT
```

Identical command to the original run (Job 23851383), restricted to scGPT via the
existing `--models` option; a `scGPT` entry was added to the script's embedding
config (`option_c/scgpt_covid_canonical.h5ad`, `X_scgpt`). Default model list is
unchanged, so other FMs are never touched. Outputs are timestamped — nothing is
overwritten. 10/30/50/70 % dimension masking + 0 % baseline, 5 repeats, 5-fold
`StratifiedGroupKFold` (the script calls `cross_val_predict` twice per repeat —
that is why it is the slowest task). **Caveat:** the script's `binary` task is
non-severe (healthy/mild/moderate/recovered) vs severe (severe+fatal), matching
the numbers in `tab:dropout` (STATE 0 % = 0.812) even though the caption says
"binary mortality"; run as-is for comparability with the other rows.

### Task 3 — embedding ablation → `tab:ablation`

```
python3 -u scripts/24_option_c_canonical_input/code/ablation_option_c.py \
    --embedding-dir 02_EMBEDDINGS/option_c \
    --output-dir 05_RESULTS/option_c_ablation --models scgpt
```

A `--models` filter was added; on a subset run the script merges into the
existing `ablation_summary.csv` (other FMs' rows preserved) and writes
`per_fm_scgpt.csv`. PCA-512 is skipped automatically for scGPT (512-dim
embedding; `tab:ablation` only reports ≤256 anyway).

**May be skippable:** scGPT full-pool ablation values were already computed in
the 2026-06-09 stage-30 batch and are bit-identical-protocol rows in
`scripts/30_review_response_stats/results/ablation_summary_fullpool.csv`
(scGPT baseline 0.649, PCA-256/128/64 = 0.649/0.647/0.645, remove-top-1/5/10 =
0.645/0.637/0.622, L2 = 0.650, z-score = 0.671, random-256/64 = 0.644/0.620).
Run task 3 only if you want the rows in the canonical `05_RESULTS/option_c_ablation/`
location (or to re-verify); otherwise fill the scGPT column from that CSV.

### Task 4 — multimodal fusion refit → `tab:fusion`

```
python3 -u scripts/24_option_c_canonical_input/code/fusion_option_c.py \
    --embedding-dir 02_EMBEDDINGS/option_c \
    --output-dir 05_RESULTS/option_c_fusion_primary5 \
    --models geneformer scgpt scprint2 state uce
```

Fusion **cannot** add a single model: every combination strategy (average, rank,
max, stacking, weighted, early-concat) uses all FMs, so the primary 5-FM set is
refit together on a fresh Cell_ID intersection. A `--models` option was added
for this; the legacy 6-FM outputs in `05_RESULTS/option_c_fusion/` are left
untouched (new output dir `option_c_fusion_primary5/`).

**Values that will change in `tab:fusion` (read before filling the table):**

- The current combination rows are actually **6-FM** values (geneformer,
  scprint2, stack, state, transcriptformer, uce — `n_fms=6` in the old
  `fusion_summary.csv`) even though the caption says "primary 5-FM". After this
  refit all combination rows become true 5-FM values and **will shift**.
- The Cell_ID intersection may grow beyond the current 43,798 cells (STACK and
  TranscriptFormer no longer limit it), so even the `single_state` /
  `single_geneformer` / `single_scprint2` / `single_uce` rows can move in the
  third decimal. **Refresh the whole table body from the new
  `fusion_summary.csv`, not just the `---` row.**
- Alternative, if you want to keep the existing combination rows and only fill
  the single-scGPT `---`: scGPT OOF on the **same 43,798-cell aligned set** was
  already computed in the 2026-06-09 stage-30 batch —
  `scripts/30_review_response_stats/results/fusion_summary_fullpool.csv` gives
  `single_scgpt = 0.638` (and 7-FM combination values; the other `single_*`
  rows there are bit-identical to the current table). In that case skip task 4
  (`sbatch --array=0-3 ...`).

## What changes in the manuscript tables

| Table | Values that change |
|---|---|
| `tab:adaptation_curve` | only the scGPT row (`---` → numbers) |
| `tab:sample_eff` | only the scGPT row |
| `tab:dropout` | only the scGPT row |
| `tab:ablation` | only the scGPT column |
| `tab:fusion` | whole body if task 4 is used (see above); only the scGPT row if filled from the existing stage-30 fullpool CSV |

All protocols are deterministic (fixed seeds: 42 everywhere; per-FM RNG in
sample efficiency), so recomputed scGPT rows are exactly what a full rerun
would have produced, and no other FM's published number is affected by tasks
0–3.

## Glue edits made (all minimal, default behavior for other FMs unchanged)

- `scripts/29_adaptation_curve/code/run_adaptation_curve_bootstrap.py`
  — scGPT added to `EMBEDDINGS`; new `--models` option; subset runs merge into
  `all_rows.csv`/`summary.csv` instead of clobbering.
- `scripts/24_option_c_canonical_input/code/primary_benchmark_option_c.py`
  — new `--models` and `--analyses` options; empty analyses no longer write
  (empty) CSVs.
- `scripts/07_gene_dropout/code/gene_dropout_robustness.py`
  — scGPT entry added to the embedding config (default model list untouched).
- `scripts/24_option_c_canonical_input/code/ablation_option_c.py`
  — new `--models` option; subset runs merge into `ablation_summary.csv`.
- `scripts/24_option_c_canonical_input/code/fusion_option_c.py`
  — new `--models` option (fusion previously had no way to restrict FMs).
- `scripts/30_review_response_stats/code/submit_scgpt_auxiliary_reruns.sh`
  — new SLURM array submitter (tasks 0–4 above).

`scgpt` was already registered in `OBSM_KEYS` of the three
`24_option_c_canonical_input` scripts, and all glob `*_covid_canonical.h5ad`,
so any future *full* rerun of these scripts now includes scGPT automatically.

## Assumptions to verify on the HPCC (first run)

- `scgpt_covid_canonical.h5ad` carries the same obs schema as the other
  canonical files (`covid_severity`, `Sample_ID`, `Cell_ID`, `Data_NO`). The
  scripts fail loudly with clear messages if a column is missing.
- `05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/adaptation_curve_bootstrap/all_rows.csv`
  and `05_RESULTS/option_c_ablation/ablation_summary.csv` exist on the HPCC
  (they back the current tables); without them the merge falls back to writing
  subset-only pooled files with a loud WARNING.
