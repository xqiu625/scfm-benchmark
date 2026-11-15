# HPCC Runbook — Revision Analyses (P1) for the Platelet scFM Benchmark

**Audience:** agent executing on UCR HPCC. Zero prior context assumed.
**Working copy on HPCC:** `/bigdata/godziklab/shared/Xinru/302006/302006_updategit`
(mirrors this repo's layout 1:1 — `scripts/NN_stage/...` etc.); the bare
`/bigdata/godziklab/shared/Xinru/302006/` above it holds the numbered data dirs
(`01_DATA`, `02_EMBEDDINGS`, `05_RESULTS`, ...). Paths below are HPCC unless marked LOCAL.
**Goal:** produce the new numbers required by the 2026-08-31 Adam review: a null
model (cell count + disease label), canonical donor-level evaluation, and an
in-distribution positive control on monocytes.

Ground rules (from `todo.md`, hard constraints):
- Do **not** modify or overwrite any existing file under `05_RESULTS/` or
  `scripts/*/results/`. All new outputs go to `scripts/31_donor_level_revision/results/`
  (small CSVs) and `05_RESULTS/revision_donor_level/` (larger artifacts).
- No resampling / SMOTE / synthetic minority oversampling anywhere.
- Do not "fix" contradictory results; report both protocols' numbers.
- SLURM conventions: see `configs/hpcc.yaml` (partition `gpu`, qos `gpu`,
  `PYTHONNOUSERSITE=1`; per-FM conda envs: Geneformer→`scgpt_finetune`,
  STATE/UCE→`state_env`, scGPT→`scpeft_env`).

Step 0 — get the new code onto HPCC:

```bash
# LOCAL (already done): scripts live in the git mirror under
#   scripts/31_donor_level_revision/code/{donor_level_reanalysis.py,null_model_cellcount.py}
# Copy them into the HPCC working copy:
scp -r /Users/xinruqiu/Desktop/302006/302006_updategit/scripts/31_donor_level_revision \
  xinruqiu@cluster.hpcc.ucr.edu:/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/
```

---

## Job A — Null model: cell count + disease label (P1-A)  [CPU, <10 min]

**Question (Adam):** "You could almost predict outcome from cell counts alone or
from identifying a disease." Quantify it.

```bash
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision
export PYTHONNOUSERSITE=1
python code/null_model_cellcount.py \
  --make-donor-meta /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c/state_covid_canonical.h5ad
```

- obs-only h5py read; applies the same disease/severity filters as stage 27
  (COVID-19+Sepsis; fatal/mild/moderate/recovered/severe) → 322 donors, 33 fatal.
  **Also export `Data_NO` (source-study) and any timepoint column if present in obs**
  (e.g. Qiu sepsis T0/T+6h) into `null_model_donor_meta.csv` — needed to auto-fill
  the P1-D cohort-characteristics table (see
  `reports/manuscript/revision/P1D_cohort_table_skeleton.md`).
- Trains LR (`C=1.0, lbfgs, balanced`) under the identical
  `StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)` on:
  (a) `log1p(n_cells)` only, (b) `log1p(n_cells)` + disease one-hot; each in
  cell-expanded and donor-level variants. Donor bootstrap B=5000, seed 42.
- Output: `results/null_model_auc.csv` (4 rows), `results/null_model_donor_meta.csv`.

**Decision rule (report back immediately, do not proceed to framing work):**
if feature set (b) cell-expanded AUC ≥ 0.7 → the null becomes a headline result;
flag it. LOCAL preliminary (cell count only): 0.659 [0.551, 0.762] / 0.665 donor-level.

**STATUS (2026-09-02): DONE.** (b) cell-expanded = 0.6385 [0.532, 0.740] — rule NOT
triggered; null stays a supplementary robustness check. Results in
`results/null_model_auc.csv` (+ `null_model_donor_meta.csv` with Data_NO).
Note: no timepoint column exists in the h5ad obs; see Job H for why that matters.

---

## Job H — Patient-level regrouping: multi-draw leakage fix  [CPU, ~1–2 h]  ⚑ PRIORITY

**Finding (2026-09-02, see `reports/manuscript/revision/patient_overlap_finding.md`):**
6 of 11 cohorts encode multiple draws of one patient as separate `Sample_ID`s
(e.g. Qiu sepsis `NS1_T0/NS1_T6` = one patient). All existing analyses group by
`Sample_ID` → same patient can straddle train/test folds. Correct counts:
**322 Sample_IDs → 238 patients; sepsis 17→12 (fatal 4→2); COVID 305→226 (fatal 29→27).**
Three Stephenson patients have draw-level mixed labels; interim rule = patient label
`max(y)` (fatal if any draw fatal), pending metadata verification.

Everything downstream of the CV split must be re-run with `patient_id` groups.
The mapping ships with this revision:
`scripts/31_donor_level_revision/results/donor_meta_with_patient_id.csv`
(columns: sample_id, disease, n_cells, y, Data_NO, patient_id).

1. **Verify the mapping** against the canonical h5ad
   (`02_EMBEDDINGS/option_c/state_covid_canonical.h5ad`, stage-27 filters):
   every post-filter Sample_ID must have exactly one patient_id; report
   n Sample_IDs / n patient_ids / n fatal patients per disease.
2. **Re-run the stage-27 mortality CV** for all 5 FMs with
   `groups=patient_id` (copy `scripts/27_covid_sepsis_only_rerun/code/rerun_covid_sepsis_mortality.py`,
   add the patient_id join; keep filters/classifiers/seeds identical) →
   OOF predictions `05_RESULTS/revision_donor_level/oof_patient_<FM>.csv`
   (same schema as before: cell_id, sample_id, severity, y, p_lr, p_rf).
3. **Re-run cross-disease transfer** (both directions, severe-excluded bootstrap
   protocol as in the manuscript) with patient-level splits; the sepsis arm is
   now 12 patients / 2 fatal — state n inline wherever results are reported.
4. **Re-run the adaptation curve** (stage 29) with the sepsis test split drawn at
   patient level; only 2 fatal sepsis patients exist, so the split must keep ≥1
   fatal patient in train and test — document the split and run 10 seeds.
5. **Re-run the revision analyses on the new OOF predictions, structured BY DISEASE
   GROUP (2026-09-02 decision — this is the manuscript's results structure):**
   - **within-COVID (PRIMARY):** patient-level OOF AUC per FM on the 226 COVID
     patients (27 fatal), mean/median/trimmed aggregation + bootstrap CI +
     pairwise paired tests. This is the headline table.
   - **within-sepsis: NOT computed** (2 fatal patients cannot support CV).
     Descriptive stats only.
   - **cross-disease transfer:** both directions, patient-level, n inline
     (sepsis: 12 patients / 2 fatal).
   - **null model + PCA-200 baseline:** same within-COVID protocol.
   ```bash
   python code/donor_level_reanalysis.py \
     --data /bigdata/.../05_RESULTS/revision_donor_level --group-col patient_id \
     --out results_patient
   python code/null_model_cellcount.py --donor-meta results/null_model_donor_meta.csv \
     --out results_patient   # collapse meta by patient_id first (as in step 2)
   ```
   (donor_level_reanalysis.py already supports `--group-col patient_id` +
   `--patient-map`; filter to COVID rows for the primary table.)
6. **Healthy-control dedup check:** apply the same patient_id rules to the HC arm
   (3,205 cells / 71 donors, excluded from the mortality task) and report the
   healthy patient count for the cohort description
   ("71 donors" → N patients). No analysis on healthy beyond this.
7. **Report back:** new OOF CSVs, the within-COVID table, transfer numbers,
   adaptation curve, per-disease patient/fatal-patient counts, healthy patient count.

**STATUS (2026-09-03): DONE (steps 1-3, 5-6 complete; step 5's PCA-200 baseline
outstanding).** sbatch jobs 28019290 (step 2, 14m30s), 28026809 (step 3, second
attempt — first attempt 28019328 failed on a real bug: off-disease cells had no
patient_id mapping; fixed by restricting the mapping assertion to used cells
only), 28019329 (step 4, 7h52m). Outputs pulled to
`scripts/31_donor_level_revision/results_patient/` (within_covid/, oof_patient/,
cross_disease/, adaptation_curve/) and left in place on HPCC under
`05_RESULTS/revision_donor_level/`.

- Step 1: 322 Sample_IDs -> 238 patients, 0 missing/duplicate. COVID-19 226/27
  fatal, Sepsis 12/2 fatal.
- Step 2: `oof_patient_<FM>.csv` for all 5 FMs, cell-level, patient_id-joined.
- Step 3: cross-disease transfer, patient-level, both directions, LR + RF,
  bootstrap CI, n inline — see `results_patient/cross_disease/`.
  scGPT forward (COVID->Sepsis) LR=0.858 is NOT comparable to the manuscript's
  0.352 — different embedding source (canonical Option C vs. the legacy
  `02_EMBEDDINGS/scgpt/covid_scgpt_clean.h5ad` path used for the manuscript
  number), unrelated to the patient-id fix.
- Step 4: adaptation curve, patient-level, 10 seeds, sepsis test split =
  1 fatal / 3 survived (train = 1 fatal / 7 survived); see
  `results_patient/adaptation_curve/summary.csv` + `split_documentation.txt`.
- Step 5: within-COVID PRIMARY table done (`results_patient/within_covid/`,
  226 patients/27 fatal); null model done (within-COVID, cell-expanded
  AUC=0.5702, rule not triggered). **PCA-200+LR baseline NOT done** —
  `pca_baseline_mortality.py` only emits fold-level AUC, not per-cell OOF, so
  it doesn't drop into `donor_level_reanalysis.py`'s aggregation without more
  surgery. Left outstanding.
- Step 6: healthy arm 3,205 cells / 71 Sample_IDs -> 70 patients (1 collision,
  Stephenson `302004data17`). See `results_patient/healthy_patient_mapping.csv`.

**STATUS CORRECTED (2026-09-03, later same day):** the STATUS above used a
buggy Stephenson patient mapping — `patient_id` was set to the `BGCVnn` 10x
sequencing pool (verified against E-MTAB-10026's SDRF: `Characteristics[sequence
pool]`), not the real patient (`Characteristics[individual]` = the `CVxxxx`
token after the pool prefix). Pools mix 4 unrelated patients each; 13 of the 15
pools mix outcomes, which is what produced the "mixed-label" warnings in the
original run. Found by a separate Job F agent while cross-referencing Stephenson
IDs for the monocyte object; independently confirmed here against the live SDRF
(not assumed) before fixing. `healthy_patient_count.py` had the identical bug
(same fix applied). Two Sample_IDs (`BGCV06_CV0201`, `BGCV13_CV0326`) carry a
pool-prefix that doesn't match either of that individual's real SDRF pools, but
the individual code itself is unambiguous either way — resolved without
guessing. Corrected `donor_meta_with_patient_id.csv` in place (overwritten;
old buggy version not preserved — git history has it if needed).

Corrected counts: 322 Sample_IDs -> **259 patients** (was 238). COVID-19:
**247/27 fatal** (was 226/27 — fatal count unchanged, only survivor count
corrected). Sepsis 12/2 fatal, unchanged (bug was Stephenson-only). Healthy
arm: 71 Sample_IDs -> 70 patients, **unchanged count** (Stephenson healthy
subset breaks down as 24->23 patients either way — the old code coincidentally
landed on the same total via a different, wrong grouping; identity is now
correct even though the number wasn't wrong).

Re-ran steps 1 (re-verify, PASS), 2 (sbatch 28052867, 13m29s), 3 (sbatch
28052868, 5m21s), 5 (within-COVID table + null model), 6 (rerun with fixed
script). Step 4 (adaptation curve) re-submitted as sbatch 28052869 (~8h ETA,
same as before) — **still running as of this note**; results not yet folded
in below. Diff vs. the buggy numbers (both fatal counts unchanged, so deltas
are all from corrected survivor grouping):

within-COVID PRIMARY, donor-level AUC (LR, mean), buggy -> corrected:
| FM | buggy (226 pts) | corrected (247 pts) |
|---|---|---|
| STATE | 0.7755 | 0.7875 |
| Geneformer | 0.7117 | 0.7079 |
| scPRINT-2 | 0.6989 | 0.6783 |
| scGPT | 0.6609 | **0.7387** |
| UCE | 0.6469 | **0.7155** |

Ranking changed: buggy had STATE > Geneformer > scPRINT-2 > scGPT > UCE;
corrected has STATE > scGPT > UCE > Geneformer > scPRINT-2. scGPT and UCE
moved up the most. Null model (within-COVID, cell-expanded): 0.5702 -> 0.5301
(still far below the 0.7 headline threshold either way). Cross-disease
COVID->Sepsis direction is numerically unchanged (Stephenson isn't in the
12-patient sepsis test set); Sepsis->COVID direction's point AUCs are also
~unchanged (computed cell-level, and the corrected mapping doesn't change
which cells exist — it only relabels patient identity) but n_test_patients
rose from 148 to 162 and bootstrap CIs narrowed accordingly (more patient
groups for the same cells). Full corrected tables in `results_patient/`.

**Also found while fixing this:** the exact same STATUS block above sat next
to a much larger "Job H2" section (added to this file after the original Job
H run, referencing `reports/manuscript/revision/literature_audit/T3_patient_sample_mapping/`) that
independently reached the same Stephenson verdict AND found three more
mapping bugs — Combes pooled-library GSMs that must be dropped, plus unresolved
Wilk/Lee repeat-draw pairs — none of which this fix touches. That audit
already computed a 9-cohort, 806-row corrected patient-sample map
(`patient_sample_map.csv`). This fix is Stephenson-only and
should be treated as superseded by a proper Job H2 run, not as the final word.

⚠ Job B's "expected values" are **superseded** — all current donor-level numbers
are pre-fix (sample-level folds). Job B as originally specified is cancelled;
Job H step 5 replaces it. Run Job F (monocytes) only AFTER Job H, and use
patient-level groups there from the start.

---

## CV Protocol v2 — repeated-CV is now the primary protocol everywhere (2026-09-04, user decision)

**Why:** the fold-composition investigation (see Job H2 STATUS below) confirmed
`StratifiedGroupKFold(shuffle=True, random_state=42)` reassigns ~80% of patients
to different folds from a 2-3 label correction — a single-split ranking is noise,
not signal. Full spec: `reports/manuscript/revision/P0_onepager.md` §2 (Estimand).

**Applies to every CV-dependent step from here forward:** Job H2 steps 2-4, Job S
step 5, Job F evaluation, the PCA-200 baseline, the null model (Job A). Each
becomes **10 seeds × 5-fold patient-level StratifiedGroupKFold = 50 fits/model**.
Cost is trivial (LR/RF on frozen embeddings, no GPU) — do not skip for speed.

**Reporting, per task:**
1. Per-repeat pooled-OOF patient-level AUC → mean ± sd across the 10 repeats
   (the headline number).
2. Pairwise ΔAUC, same seed paired across models → mean Δ + 95% CI across the
   10 seeds, Wilcoxon signed-rank as secondary.
3. The seed-42 single split (with its donor-clustered bootstrap) is kept only as
   a clearly-labeled secondary view — not the headline.

The post-fix (Step-0 + v2 mapping) labels stay; only the AUC-magnitude reporting
changes. **The single-split Job H table (STATE 0.7755 etc.) and the pre-repeated-CV
Job H2 tables (both the 27-fatal and 25-fatal versions) are superseded — don't
report them further**, keep only as historical record in this file.

One Methods sentence (protocol justification, not a discovery narrative):
"Fold assignments under shuffled stratified group CV are unstable to small label
changes (~80% of patients change folds after a two-label correction); all results
are therefore reported as mean ± sd over 10 repeated 5-fold splits."

---

## Job H2 — v2 patient mapping + PCA-200 baseline + within-COVID refresh  [CPU, ~2 h]  ⚑ PRIORITY

**Why (2026-09-03):** Job H ran on the v1 mapping. The Task-3 audit
(`reports/manuscript/revision/literature_audit/T3_patient_sample_mapping/`, verdicts in
`patient_overlap_finding.md` §Task-3) proved v1 defective in four places:
(a) Stephenson rule was INVERTED — `BGCVnn` is the 10x sequencing pool (4 patients
per pool), the trailing `CVxxxx` is the patient; (b) 3 Combes GSMs are
multi-patient pooled libraries and must be excluded; (c) Wilk has one repeat-draw
pair; (d) Lee has three repeat-draw pairs. v1's Stephenson error was conservative
(over-grouping, no leakage), but Combes pools are a data-validity issue.

**Step 0 — label corrections (Task-4 audit, 2026-09-03; do BEFORE rebuilding anything):**
three Combes samples carry wrong labels in the canonical h5ad obs
(`covid_severity`), verified against Combes 2021 Supplementary Table 1 p.3
(see `reports/manuscript/revision/severity_label_audit.md` §Task 4):
- `GSM4995436` (Patient 29 / COMET 1029): paper = Mild/Moderate, no ICU,
  0 vent days, Discharged=Yes (SURVIVED) → `covid_severity`: fatal → **moderate**
  (y flips 1→0; 324 cells)
- `GSM4995442` (Patient 47 / COMET 1047): paper = Discharged=Yes (SURVIVED),
  SEVERE arm → `covid_severity`: fatal → **severe** (y flips 1→0; 159 cells)
- `GSM4995432` (Patient 14 / COMET 1014): paper = Mild/Moderate arm →
  `covid_severity`: severe → **moderate** (y stays 0)
Apply the obs edit wherever `covid_severity` is stored for these Sample_IDs
(canonical option_c h5ads and any derived subset), then recompute y downstream.
New fatal totals after Step 0: COVID 27 fatal sid / 25 fatal patients
(was 29 sid / 27 patients); whole-cohort fatal patients 27 (COVID 25 + sepsis 2).
Document the before/after in the diff report.

**Step 1 — rebuild the mapping → `donor_meta_with_patient_id_v2.csv` (+ diff report):**
- `302004data17` (Stephenson): for IDs matching `BGCV\d+_.*`, patient_id = the
  trailing token after the last `_` (e.g. `BGCV02_CV0052` → `CV0052`); same trailing
  id across different BGCV prefixes = same patient (D0/D28 resample). MH*/newcastle*/
  AP* ids unchanged.
- `302005data02` (Combes): **drop GSM4995425, GSM4995426, GSM4995430 entirely**
  (pooled libraries; 2,597 cells) — remove from the mapping AND from all downstream
  cell sets.
- `302004data01` (Wilk): one repeat pair exists (9d/11d, straddling intubation) but
  the `covid_5xx`↔C-label mapping is unstated; if not resolvable from the
  integration scripts within 1 h, keep 7 and FLAG — do not guess.
- `302004data02` (Lee): merge pairs (nCoV 3&4, 7&8, 9&10) → 8 patients.
- Report per-disease patients/fatal before-vs-after (expect COVID ≈ 240 patients,
  sepsis unchanged 12/2, healthy recount with the Stephenson fix).

**PROTOCOL CHANGE (2026-09-04, CV-instability finding):** the single-split CV is
no longer the primary protocol anywhere. Confirmed on HPCC: flipping 2–3 labels
reshuffles ~80% of patients across folds (StratifiedGroupKFold(shuffle=True) is
unstable to stratification-target changes), so single-split rankings are noise.
**All CV-dependent steps below use REPEATED-CV: 10 seeds × 5-fold patient-level
StratifiedGroupKFold = 50 fits per model/task.** Report per task:
(i) per-repeat pooled-OOF AUC distribution → mean ± sd;
(ii) pairwise ΔAUC per seed with the SAME seed paired across models → mean Δ +
95% CI across the 10 seeds (+ Wilcoxon signed-rank as secondary);
(iii) keep the donor-clustered bootstrap on the seed-42 repeat only as a
secondary, clearly-labeled single-split view. Cost is trivial (LR/RF on frozen
embeddings; no GPU). The existing 8-seed stability analysis (stage 30) is the
precedent — this makes it the primary protocol.

**Scope (2026-09-04, agent pushback accepted):** repeated-CV applies to
**fold-based** evaluations only (within-COVID, within-sepsis, PCA/null baselines
under StratifiedGroupKFold). Two analyses are NOT fold-based and keep their own,
already-adequate uncertainty treatments — do not force 10 seeds onto them:
- **cross-disease transfer**: fixed train/test split; uncertainty = existing
  200-iteration donor bootstrap (test-side resampling). LR is deterministic
  given the split, so seeds would be a literal no-op. Train-side refit variance
  is second-order and not included (one runbook/M&M sentence). OPTIONAL, nearly
  free: 10 RF `random_state`s for the RF column only (adds real, if minor,
  variance where a stochastic classifier is actually reported).
- **adaptation curve**: already runs a 10-seed loop over which patients enter
  the training pool — that IS its repeated design (subsampling variance). No
  changes; do not re-run (7–8 h for zero gain).
Manuscript protocol language: three uncertainty treatments, each matched to its
analysis design — repeated-CV (fold-based), donor bootstrap (fixed-split
transfer), training-pool subsampling (adaptation curve).

**Step 2 — rerun stage-27-style CV** with `groups=patient_id_v2` for all 5 FMs on
the pool-excluded cohort, **under repeated-CV (10 seeds × 5-fold, seeds
42..51)** → per-seed OOF preds `oof_patient_<FM>_seed<NN>.csv` to
`05_RESULTS/revision_donor_level_v2/` (new dir; do NOT touch
`revision_donor_level/`).

**Step 3 — PCA-200+LR baseline fix (the outstanding Job H item):** modify
`scripts/04_pca_baseline/code/pca_baseline_mortality.py` to emit per-cell OOF
probabilities (same schema as `oof_patient_<FM>.csv`: cell_id, sample_id,
severity, y, p_lr) under the identical repeated-CV protocol (PCA-200 on raw
counts; LR C=1.0 balanced; same 10 seeds).

**Step 4 — rerun `code/donor_level_reanalysis.py`** on the v2 OOFs
(`--group-col patient_id --patient-map donor_meta_with_patient_id_v2.csv`,
aggregating per repeat) → within-COVID primary table: 5 FMs **+ PCA-200 + null
model** in one table, reported as repeat-mean ± sd with per-seed paired
pairwise tests. If cheap in the same pass, also refresh cross-disease transfer
and the adaptation curve on v2 under repeated-CV (COVID-side pools changed;
sepsis arm unchanged).

**Report back:** v2 mapping CSV + diff, the within-COVID table incl. PCA-200
(repeat-mean ± sd), per-seed OOF files, healthy patient recount, and (if run)
v2 transfer/adaptation numbers.
**Then the title/abstract decision unblocks.**

**Job H2 STEP 0 + RERUN STATUS (2026-09-03, later same day): DONE.** Applied
the 3 Combes label corrections at read-time in every v2 script (not by
mutating the shared canonical h5ads, which Job F may be reading concurrently)
via a small `SEVERITY_FIX` dict in `rerun_covid_sepsis_mortality_patient_v2.py`,
`pca_baseline_patient_oof_v2.py`, `cross_disease_fatal_vs_survived_patient_v2.py`;
also hand-patched `donor_meta_with_patient_id_v2.csv`'s `y` column for the 2
flipped-fatal sample_ids (GSM4995436, GSM4995442: 1→0). Verified count:
COVID 27→**25** fatal patients, whole-cohort 29→**27** — matches the runbook's
expected numbers exactly. Reran steps 2 (sbatch 28060400, 14m38s), 3 (sbatch
28060401, 4m10s), 4.

**⚠ Ranking changed materially, not just numbers — report this prominently:**
within-COVID PRIMARY (240 patients/**25 fatal**), donor-level AUC (LR/mean):
**Geneformer 0.7849 (now best)**, scPRINT-2 0.7358, STATE 0.7113 (down from
0.8025, now 3rd), PCA-200 0.6608, scGPT 0.6199, **UCE 0.5953 (down from 0.7054,
now worst, barely above chance)**. STATE vs PCA-200: dAUC 0.0504, q=0.107 (still
not significant, but the point estimates moved). Null model: 0.5517 (still far
below 0.7). Full tables: `results_patient_v2/within_covid/`. **This is a genuine
result of fixing 2 wrong fatal labels, not noise** — flipping 2 of 27 fatal
patients (7.4% of the fatal class) in a 240-patient/25-fatal regime is exactly
the kind of small-n perturbation that can reorder close AUCs; STATE and UCE
were evidently more sensitive to those 2 specific patients than Geneformer/
scPRINT-2 were. `cross_disease_fatal_vs_survived_patient_v2.py` re-fixed and
re-run (sbatch 28121440, 4m11s) — relabeled cross-disease numbers now current:
COVID->Sepsis LR: scGPT 0.854, scPRINT-2 0.840, STATE 0.813, UCE 0.778,
Geneformer 0.767. Sepsis->COVID LR: STATE 0.723, scGPT 0.666, Geneformer 0.626,
UCE 0.654, scPRINT-2 0.545. Full table:
`results_patient_v2/cross_disease/cross_disease_patient_20260904_110222.csv`.

**Job H2 STATUS (2026-09-03): DONE except the adaptation curve (deliberately
deferred, not "cheap" — v1's took 7h52m; not re-run here).** v1's Stephenson-only
STATUS CORRECTED block above stays as historical record.

- **Step 1:** `donor_meta_with_patient_id_v2.csv` built and verified against the
  canonical h5ad (319/319 mapped sample_ids = 322 minus the 3 intentionally-
  dropped Combes GSMs, exactly; 0 duplicate patient_id). Wilk turned out simpler
  than expected — our Sample_IDs already use the paper's own C-labels (`C1 A`,
  `C1 B`, `C2`, ...) rather than ambiguous `covid_5xx` GEO ids, so the repeat
  pair merges with a plain suffix strip; nothing to flag. Lee merged as specified
  (nCoV 3&4, 7&8, 9&10; no y conflicts). Combes dropped GSM4995425/26/30
  (392+711+1494 = 2,597 cells exactly). Diff report (8 rows):
  `results_patient/patient_map_v1_to_v2_diff.csv`.
  Counts, v1→v2: total 259→**252**; COVID 247→**240** (fatal unchanged, 27);
  Sepsis 12/2 unchanged. Matches the "≈240" expectation exactly.
  Healthy recount: **unchanged at 70** — checked directly (not assumed): Wilk
  and Combes have zero healthy-arm cells, Lee's healthy samples use `Normal_1..4`
  naming (not `nCoV_*`), so none of the v2 rules touch the healthy arm.
- **Step 2:** sbatch 28058454 (13m58s) → `oof_patient_<FM>.csv` for all 5 FMs,
  `05_RESULTS/revision_donor_level_v2/` (new dir, v1's dir untouched). New
  script `code/rerun_covid_sepsis_mortality_patient_v2.py` (copy of the v1
  script with `PATIENT_MAP`→v2 and an inner-join instead of an asserted-zero-
  missing left-join, so v2-excluded cells drop cleanly instead of crashing).
- **Step 3:** the runbook's cited `scripts/04_pca_baseline/code/pca_baseline_mortality.py`
  doesn't exist — the real file is `scripts/24_option_c_canonical_input/code/pca_baseline_mortality.py`,
  and (separately from the OOF gap) its `valid` mask has no Disease filter, so
  without a fix it would silently pull in Flu/Lung_infect/Non_covid/SLE cells
  sharing a severity label — added the same Disease∈{COVID-19,Sepsis} filter the
  FM runs use. New script `code/pca_baseline_patient_oof_v2.py`: PCA-200 on raw
  counts (`01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad`), scaler+PCA
  fit per-fold on train only, StratifiedGroupKFold(patient_id_v2). sbatch 28058493
  (3m48s) → `oof_PCA-200.csv`. p_rf column is a duplicate of p_lr (LR-only
  baseline, no RF variant was requested) — flagged in the script docstring so it
  isn't mistaken for a real RF fit later.
- **Step 4:** `donor_level_reanalysis.py` patched with a new `--extra-methods`
  flag (comma-separated, appended after the 5 primary FMs; empty by default, so
  Job H/Job B's existing behavior is unchanged) to fold PCA-200 into the same
  table. Within-COVID PRIMARY (240 patients/27 fatal), donor-level AUC (LR/mean):
  STATE 0.8025, Geneformer 0.7816, PCA-200 0.7700, scPRINT-2 0.7574, scGPT
  0.7338, UCE 0.7054. **STATE vs PCA-200: dAUC 0.0325, q=0.287 — NOT
  significant** (directly relevant to the manuscript's OOD claim). Null model
  (within-COVID, cell-expanded): 0.5015 (essentially chance; lower than v1's
  0.5301 and further below the 0.7 threshold). Full tables:
  `results_patient_v2/within_covid/`.
- **Cross-disease (opportunistic, "if cheap" — v1 took 5m so v2 was run too):**
  new `code/cross_disease_fatal_vs_survived_patient_v2.py` (v2 map + explicit
  exclusion of the 3 Combes sample_ids from the `valid` mask, since dropping
  them from the mapping alone isn't enough — the original script's `valid` is
  severity-based and doesn't know about the v2 exclusion list). sbatch 28059149
  (4m47s). Point AUCs shifted modestly from v1 (unlike the Stephenson-only fix,
  this one changes the actual COVID training-cell set, not just patient labels)
  — e.g. STATE COVID→Sepsis 0.7975→0.7760. Full table:
  `results_patient_v2/cross_disease/`.
- **Adaptation curve: NOT re-run.** v1's took 7h52m; explicitly deferred per
  the "if cheap" hedge rather than committing another ~8h without a fresh
  go-ahead. v1's Stephenson-only numbers stand as the most recent adaptation-
  curve result; flag as stale pending a dedicated v2 run if needed.
- **Severity/disease-consistency audit (a coordinator side-request during this
  work, using the v1 mapping):** 259 patients, 207 single-draw, 52 multi-draw;
  49 consistent, 3 flagged (all Bernardes 004/005/012, `TA/TA2/TB/TC/TE`
  visit-code draws showing moderate→severe→recovered-type sequences) — all
  read as genuine clinical trajectory per Bernardes' documented visit schedule,
  none fatal, none Disease-mismatched. `results_patient/severity_consistency_audit.csv`.
  **Independently corroborated** by a separate, more thorough audit that
  appeared in `reports/manuscript/revision/severity_label_audit.md` during this
  same work session (not authored by this fix): confirms the same 3 Bernardes
  patients as benign draw-level variation, and goes further — cross-referencing
  21/33 fatal labels directly against GEO/repository records (Ren 9, Stephenson
  8, Qiu 4, all 100% match), with 12/33 fatal labels still open pending source-
  paper supplement checks (Schulte-Schrepping 4, Liu 4, Bernardes 2, Combes 2),
  and one open severity-label discrepancy unrelated to Job H2's scope
  (Combes `GSM4995432`: labeled severe here, `M/M` in the repository record).
  Worth reading directly rather than re-summarized further here.

**CV PROTOCOL v2 RERUN STATUS (2026-09-04): DONE.** New shared module
`code/repeated_cv.py` (importable, no CLI -- `run_repeated_kfold`,
`repeated_patient_auc`, `summarize_repeated`, `paired_delta_auc`,
`DEFAULT_SEEDS=42..51`); `rerun_covid_sepsis_mortality_patient_v2.py`,
`pca_baseline_patient_oof_v2.py`, `null_model_cellcount.py`,
`donor_level_reanalysis.py` all gained `--n-seeds N` (default 1 = old
single-split behavior, unchanged). Job F reuses this same module for its own
evaluation.

**⚠ The single-split numbers immediately above (STATE 0.7113, Geneformer 0.7856
best, etc.) are SUPERSEDED — do not report them further, historical record
only.** The repeated-CV numbers below are the new headline, and they change
the story, which is exactly why this protocol exists:

**Within-COVID PRIMARY, repeated-CV (240 patients/25 fatal, 10 seeds x 5-fold),
repeat-mean +/- sd AUC (LR, mean aggregation):**
| FM | repeat-mean | sd | seed-42 (secondary) |
|---|---|---|---|
| Geneformer | 0.7856 | 0.0264 | 0.7849 |
| **STATE** | **0.7778** | 0.0329 | 0.7113 |
| scPRINT-2 | 0.7402 | 0.0270 | 0.7358 |
| PCA-200 | 0.7251 | 0.0326 | 0.6608 |
| scGPT | 0.7196 | 0.0474 | 0.6199 |
| UCE | 0.6611 | 0.0461 | 0.5953 |

**STATE moves from 3rd (single-split, noise-dominated) to a statistical tie
for 1st under repeated-CV** -- STATE vs Geneformer: mean-delta -0.0079,
95% CI [-0.0361, 0.0204], Wilcoxon p=0.48 (NOT significant; q=0.51). This is
the single-split instability the protocol change was built to catch: seed 42
happened to land STATE in an unlucky fold configuration.

**STATE vs PCA-200 is now SIGNIFICANT** (it was not, under any single-split
version so far): mean-delta 0.0526, 95% CI [0.0450, 0.0602], Wilcoxon
p=0.00195, q=0.00325. Directly relevant to the manuscript's OOD claim --
under the robust protocol, STATE does beat the simple baseline.

Full pairwise table (all 15 pairs, mean delta + 95% CI + Wilcoxon p + BH-q):
`results_patient_v2/within_covid/donor_level_pairwise_repeated.csv`. Per-seed
raw AUCs (audit trail): `per_seed_auc_repeated.csv`. Headline table:
`donor_level_auc_repeated.csv`.

**Null model, repeated-CV:** 0.5544 +/- 0.0458 (cell-expanded), 0.5542 +/- 0.0115
(donor-level) -- consistent with the single-split 0.5517; decision rule still
not triggered. `results_patient_v2/within_covid/null_model_auc_repeated.csv`.

**Fix-1/Fix-2/Fix-3 STATUS (2026-09-04, code review): DONE.**

- **Fix-1:** `build_patient_map_v2.py` now bakes the Task-4 SEVERITY_FIX
  directly into `donor_meta_with_patient_id_v2.csv` (previously applied only
  at read-time in the CV/PCA scripts, which meant the CSV itself still had
  y=1 for GSM4995436/442). Regenerated: COVID 25 fatal patients, whole-cohort
  27 -- matches expected. **The within-COVID FM+PCA-200 repeated-CV table
  above needed NO rerun** -- it computes y from h5ad severity with the
  read-time fix already applied, independent of the CSV. Only the null model
  read y from the CSV directly; rebuilt its collapsed meta and reran:
  0.5544±0.0458 (cell-expanded), numerically unchanged from before this fix.
- **Fix-2/Fix-3:** new shared `transfer_utils.py` -- patient-level primary AUC
  (mean cell score per patient -> AUC; bootstrap on the patient-level array)
  with cell-weighted kept as a secondary column; severe=survived is now the
  primary COVID protocol, severe-excluded a labeled sensitivity variant.
  Both `cross_disease_fatal_vs_survived_patient_v2.py` (legacy 12-patient
  sepsis arm) and `jobS_cross_disease.py` (Chongqing) rewritten to use it.

**Cross-disease transfer, patient-level primary (Fix-2/3), legacy 12-patient
sepsis arm, ALL 5 FMs** (`05_RESULTS/revision_donor_level_v2/cross_disease_patient_20260904_233237.csv`),
LR, severe-included primary / severe-excluded sensitivity:
| FM | COVID→Sepsis primary (n=12/2 fatal) | sens (severe-excl) | Sepsis→COVID primary (n=240/25 fatal) | sens |
|---|---|---|---|---|
| STATE | 1.000 [1.0,1.0] | 1.000 | 0.781 [0.697,0.865] | 0.862 |
| Geneformer | 1.000 [1.0,1.0] | 1.000 | 0.794 [0.718,0.878] | 0.863 |
| scGPT | 0.950 [0.727,1.0] | 1.000 | 0.751 [0.659,0.836] | 0.819 |
| scPRINT-2 | 0.900 [0.549,1.0] | 1.000 | 0.601 [0.487,0.724] | 0.646 |
| UCE | 0.600 [0.200,0.909] | 1.000 | 0.682 [0.584,0.798] | 0.755 |

**COVID→Sepsis with only 12 test patients (2 fatal) is extremely small-n --
patient-level aggregation collapsing to 0.60-1.00 across FMs (several exactly
1.0) is a genuine artifact of tiny-n instability, not a bug** (cell-weighted
numbers on the same data were the more moderate 0.53-0.86 range reported
earlier) -- reported as-is per the "no adjudication" writing rule.
Sepsis→COVID (n=240, both classes well-populated) is far more stable.

**Cross-disease transfer, patient-level primary (Fix-2/3), Chongqing sepsis
arm, ALL 5 FMs** (`05_RESULTS/revision_donor_level_v2/cross_disease_chongqing_20260904_233237.csv`):
| FM | COVID→Sepsis primary (n=91/17) | sens | Sepsis→COVID primary (n=240/25) | sens |
|---|---|---|---|---|
| UCE | 0.684 [0.543,0.805] | 0.671 | 0.607 [0.488,0.728] | 0.653 |
| STATE | 0.653 [0.490,0.808] | 0.657 | 0.689 [0.583,0.785] | 0.735 |
| Geneformer | 0.644 [0.486,0.824] | 0.600 | 0.674 [0.557,0.784] | 0.728 |
| scGPT | 0.650 [0.495,0.780] | 0.617 | 0.624 [0.494,0.753] | 0.641 |
| scPRINT-2 | 0.647 [0.496,0.791] | 0.595 | 0.590 [0.479,0.713] | 0.606 |

With real n (91 vs 12), no extreme values -- consistent, modest transfer
signal across all 5 FMs, in line with the within-sepsis null finding above.

**Fix-4 (platelet-monocyte co-expression, code review 2026-09-04): DONE.**
Computed post-hoc on the already-saved platelet subset (no integration
rerun needed): CD14/LST1/FCN1 co-expression in **3,181/22,261 platelet cells
(14.3%)**, consistent with the design doc's single-sample pre-check. Merged
into `results_patient/chongqing_platelet_counts_per_donor.csv`.

**Fix-5 (number sync):** done throughout this file and `sepsis_arm_expansion.md`
-- "108 donors/20 fatal" -> **91 donors/17 fatal** everywhere; adaptation-curve
pool 120 -> **103** (12+91).

**UPDATE (2026-09-05): scGPT landed.** The ~12h GPU-queue wait was avoidable --
scGPT was queued on the congested `gpu` partition (4/9 A100 nodes down) when
`short_gpu` (already used successfully by the other 4 FMs) was nearly empty.
Cancelled and resubmitted on `short_gpu` with a generic `--gres=gpu:1`;
completed in 40s (20,873/28,841 genes matched). scGPT is now included in
within-sepsis CV, PCA-200, and both cross-disease tables above/below --
none of these sections say "4/5 FMs" or "pending" anymore as of this update.
Job F independently hit the same `gpu`-partition congestion for its own
scGPT job and applied the same `short_gpu` fix.

**Adaptation curve: DONE** (see full results below the within-sepsis table).

---

## Job B — Canonical donor-level evaluation (P1-B) + cell-count association (P1-C)  [CPU, ~15 min]

A LOCAL run against the canonical OOF predictions already produced the target
numbers (see §Expected). This job is the canonical reproduction on HPCC.

1. Verify the OOF inputs on HPCC match the LOCAL canonical copies (md5 must match;
   if the HPCC copies differ or are missing, scp the LOCAL ones up — they were
   originally produced by stage 27 and collected in stage 30):

```
LOCAL scripts/30_review_response_stats/data/:
oof_STATE.csv       8d4df22e1f237ee0842ada16982b5900   (42,095 cells + header)
oof_UCE.csv         10876659a06d6d2d19eeedf2cb4a3ff4   (42,095)
oof_scGPT.csv       673f81664d4b63a1a12a27692542b4f0   (42,095)
oof_Geneformer.csv  f0e6652a5fa67509bba127ada68b3b64   (41,320)
oof_scPRINT-2.csv   fd8c7439bc7be6b575186ed812f8f419   (41,320)
```

2. Run:

```bash
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision
python code/donor_level_reanalysis.py \
  --data ../30_review_response_stats/data --out results
```

**Expected (must reproduce exactly; same inputs, seed 42, B=5000):**

donor-level AUC, LR, mean aggregation (primary metric):
| FM | AUC | 95% CI |
|---|---|---|
| STATE | 0.7944 | [0.6988, 0.8763] |
| Geneformer | 0.7729 | [0.6813, 0.8563] |
| scPRINT-2 | 0.7435 | [0.6472, 0.8333] |
| scGPT | 0.7431 | [0.6415, 0.8386] |
| UCE | 0.6780 | [0.5689, 0.7785] |

Paired donor-level ΔAUC: STATE−Geneformer 0.0215 (n.s., q=0.34);
STATE−scPRINT-2 0.0509 (q=0.11); STATE−scGPT 0.0513 (q=0.044);
STATE−UCE 0.1164 (q<0.001). I.e. **top-3 are a statistical tie at donor level;
STATE separates only from scGPT and UCE.**

P1-C: Spearman(n_cells, donor score) = 0.22–0.29 (p ≤ 5e-5) for every FM —
a real but weak count association; AUC within count tertiles stays 0.62–0.81.

Outputs: `results/donor_level_auc.csv`, `donor_level_pairwise.csv`,
`cellcount_association.csv`, `donor_predictions_<FM>.csv`.

---

## Job F — In-distribution positive control: monocytes (P1-F)  [GPU, multi-day]

**Why:** the OOD claim ("FMs never saw platelets, so matching them with PCA is
meaningful") is refutable unless the SAME pipeline shows FMs beat the PCA-200+LR
baseline on a cell type they WERE trained on. This is the single most important
new experiment.

### Phase 1 — inventory (go/no-go), CPU only, ~1 h

The platelet cohort was integrated from PBMC studies by Qiu et al. 2024 IJMS
(`scripts/data/integrate_platelet_data_harmony.R`); the full objects may still
exist. Find monocyte-containing data from the SAME source studies (accessions in
`docs/cohort-sources.md`: Ren GSE158055, Stephenson E-MTAB-10026,
Schulte-Schrepping EGAS00001004571 [controlled access — flag if unavailable],
Combes GSE163668, Bernardes FastGenomics, Qiu sepsis GSE167363, ...).

```bash
# candidate locations — adjust after ls
find /bigdata/godziklab/shared/Xinru/302006/data \
     /bigdata/godziklab/shared/Xinru/302006/01_DATA* -name '*.h5ad' 2>/dev/null
# for each candidate: print obs columns + any cell-type column crosstab, e.g.
python - <<'PY'
import h5py, sys
p = sys.argv[1]
with h5py.File(p, 'r') as h:
    print(p, list(h['obs'].keys()))
PY
```

GO criteria: monocytes (or CD14+ PBMC) identifiable for donors overlapping the
322-donor platelet cohort, with both fatal and survivor donors represented,
≥ ~20k cells total. If the original pre-platelet-subset objects are gone:
re-download from GEO/ArrayExpress (Ren, Stephenson, Combes are public;
Schulte-Schrepping is controlled — skip and document). **Report the inventory
table before Phase 2.**

### Phase 2 — build monocyte canonical object, CPU

Subset to monocytes; keep the SAME obs schema as option C
(`Sample_ID, Disease, covid_severity, Data_NO`, raw counts in X). Labels come
from the same donor metadata as the platelet object. Write
`02_EMBEDDINGS/option_c_monocyte/monocyte_canonical_raw.h5ad`.

### Phase 3 — embed 5 FMs, GPU (sbatch)

Use the existing per-FM entry points under `scripts/foundation_models/<FM>/`
(e.g. `Geneformer/run_geneformer_covid_option_c.sh`, `UCE/uce_human_platelet.sh`)
with input/output paths switched to `option_c_monocyte`. Same checkpoints as the
platelet run (Geneformer-V2-316M; scPRINT-2 small-v2; etc.). Per-FM envs per
`configs/hpcc.yaml`. Outputs:
`02_EMBEDDINGS/option_c_monocyte/<fm>_monocyte_canonical.h5ad`.

### Phase 4 — evaluation, CPU

1. Stage-27-style mortality rerun on monocytes
   (`scripts/27_covid_sepsis_only_rerun/code/rerun_covid_sepsis_mortality.py`,
   paths switched) → OOF predictions `oof_monocyte_<FM>.csv`.
2. PCA-200 + LR baseline under the identical protocol (see
   `scripts/04_pca_baseline/`), including a from-scratch PCA on monocyte counts.
3. `donor_level_reanalysis.py --data <monocyte oof dir>` → donor-level AUCs.

### Phase 5 — deliverable

`05_RESULTS/revision_donor_level/indistribution_vs_platelet.csv`:
rows = FM + PCA-200 baseline, columns = donor-level AUC on monocytes
(in-distribution) vs platelets (OOD), with CIs. **Send back** this file plus all
`scripts/31_donor_level_revision/results/*.csv`.

**STATUS (2026-09-07): DONE, all 5 FMs.** Monocyte canonical object built from
Ren GSE158055 + Stephenson E-MTAB-10026 (curated_atlas_label, majorType/
full_clustering monocyte categories) + Qiu sepsis GSE167363 (marker_gene_threshold:
CD14/LYZ/FCN1/S100A8/S100A9 vs. competing lineage scores — flagged distinctly via
`obs['monocyte_source_provenance']`, not blended silently). Combes/Bernardes/
Schulte-Schrepping skipped (no usable annotation or controlled access — see
checkpoint-1 inventory). 248,240 cells, 148 patients (19 fatal): Ren 53pt/9fatal,
Stephenson 90pt/8fatal, Qiu-sepsis 5pt/2fatal. Object + all 5 FM embeddings
(STATE/UCE/scPRINT-2/scGPT/Geneformer) + raw counts in
`02_EMBEDDINGS/option_c_monocyte/`. Evaluated under CV Protocol v2 (10 seeds x
5-fold patient-level StratifiedGroupKFold), same classifier config as the
platelet-side v2 scripts. Full tables:
`scripts/31_donor_level_revision/results_patient_v2/monocyte/` (also
`05_RESULTS/revision_donor_level_v2/monocyte/`).

Geneformer needed its own fix: the `geneformer_py310` conda env used by the
original platelet embedding is broken (package not installed); the
`geneformer.sif` container is a known-deprecated dead end (scipy/h5py too old
for backed-mode h5ad reads — confirmed after 5 failed attempts, matches
`docs/runbooks/hpcc.md`). Working approach (same as Job S's Chongqing
Geneformer embedding): native `scpeft_env` + `BertModel.from_pretrained()`
against `03_MODELS/geneformer/geneformer-12L-30M`, custom rank-value
tokenization against `03_MODELS/scPEFT/geneformer_peft/geneformer/
token_dictionary.pkl` (Ensembl-ID-keyed — needed the symbol->Ensembl mapping
already built via mygene.info; Chongqing's raw input happened to already have
Ensembl var_names, which is why its script didn't need this step). Script:
`code/geneformer_monocyte_native.sbatch`.

**In-distribution (monocyte) vs. OOD (platelet, within-COVID v2, 240pt/25fatal)
— LR/mean, repeat-mean +/- sd, COVID-only cut (Ren+Stephenson, 143pt/17fatal,
disease-matched to the platelet comparison) as primary:**

| FM | Monocyte (COVID-only) | Monocyte (all, incl. 5 sepsis pt) | Platelet (OOD) | Delta (COVID-only) |
|---|---|---|---|---|
| STATE | 0.6864 +/- 0.0396 | 0.6964 +/- 0.0370 | 0.7778 +/- 0.0329 | -0.0914 |
| Geneformer | 0.6665 +/- 0.0421 | 0.6741 +/- 0.0465 | 0.7856 +/- 0.0264 | -0.1191 |
| scPRINT-2 | 0.6552 +/- 0.0647 | 0.6723 +/- 0.0577 | 0.7402 +/- 0.0270 | -0.0850 |
| scGPT | 0.6642 +/- 0.0340 | 0.6765 +/- 0.0370 | 0.7196 +/- 0.0474 | -0.0554 |
| UCE | 0.6564 +/- 0.0328 | 0.6475 +/- 0.0561 | 0.6611 +/- 0.0461 | -0.0047 |
| PCA-200 | 0.6537 +/- 0.0300 | 0.6535 +/- 0.0432 | 0.7251 +/- 0.0326 | -0.0714 |

**(a) Designed question (FM vs. PCA-200, in-distribution) has a model-specific
answer, cut-dependent for one FM:** paired vs. PCA-200, COVID-only —
STATE +0.0327 (p=0.00195, sig), scGPT +0.0105 (p=0.0371, sig), Geneformer
+0.0128 (p=0.084, n.s.), scPRINT-2 +0.0014 (n.s.), UCE +0.0027 (n.s.). In the
all-cohort cut Geneformer's delta vs PCA-200 IS significant (+0.0205, q=0.015)
— the extra 5 sepsis patients tip it over the line; flagged as cut-dependent
rather than picking the answer that looks cleaner. STATE and scGPT robustly
beat PCA-200 in-distribution under both cuts; scPRINT-2 and UCE do not under
either.

**(b) The monocyte<platelet inversion hits PCA-200 equally — the key internal
control.** PCA-200's own delta (-0.0714) sits in the same range as most FMs'
deltas (STATE -0.0914, Geneformer -0.1191, scPRINT-2 -0.0850) and exceeds
scGPT's (-0.0554) and UCE's (-0.0047). PCA-200 has no pretraining at all, so it
cannot have an "OOD penalty" — its presence in this inversion, at a magnitude
comparable to or exceeding several FMs, means the monocyte<platelet inversion
is a property of the task/cohort (this monocyte assembly vs. this platelet
cohort), not a property of FM pretraining mismatch.

**(c) Composition caveat addressed, not explained away.** COVID-only deltas
match the all-cohort deltas closely (STATE -0.0914 vs -0.0814, Geneformer
-0.1191 vs -0.1115, scPRINT-2 -0.0850 vs -0.0679, scGPT -0.0554 vs -0.0431,
UCE -0.0047 vs -0.0136, PCA-200 -0.0714 vs -0.0716) — restricting to a
disease-matched pure-COVID monocyte cohort does not make the inversion go away.
The 5-patient sepsis admixture was not responsible for it.

**Source-stratified AUC** (re-aggregation only, no new CV;
`monocyte_source_stratified_auc.csv`): Ren (n=53/9fatal) and Stephenson
(n=90/8fatal) — both curated_atlas_label — diverge substantially per FM (e.g.
scGPT: 0.769 Ren vs 0.541 Stephenson; UCE: 0.705 Ren vs 0.563 Stephenson) —
real cohort/batch heterogeneity between the two curated sources themselves,
unrelated to label-provenance quality. Qiu-sepsis (marker_gene_threshold,
n=5/2fatal) has enormous sd (0.27-0.38 across FMs, PCA-200 even below chance
at 0.333) — consistent with pure small-n noise, not a distinguishable
provenance signal; **not informative on its own**, kept for descriptive
completeness only.

---

## Job S — Sepsis arm expansion: GSE279452 (Chongqing) integration  [network ~5 GB + CPU + GPU, 1–2 days]

**Decisions (user, 2026-09-03):** endpoint = **28-day mortality (mort28)** primary,
mort90 sensitivity (one Methods sentence on the in-hospital vs 28-day difference);
cohort = **third sepsis source** (`Data_NO=chongqing2025`, donors kept distinct from
Qiu/Jiang); adult scRNA arm only — pediatric, ICU-control, healthy-control all excluded.
Design doc: `reports/manuscript/revision/sepsis_arm_expansion.md`.
Scout + mortality table: `reports/manuscript/revision/literature_audit/sepsis_fatal_scout/`.

**Step 1 — download (can run in parallel with Job H2).** From GSE279452 pull the
MTX triplets (barcodes/features/matrix) for the **91 adult-sepsis scRNA donors**
(donor list: `chongqing_sepsis_donor_mortality.csv` filtered `cohort=adult &
group=='Adult Sepsis' & scRNAseq=1`, 93 rows, minus 2 with missing mort28 -> 91;
GSM↔donor_id via GEO sample titles, e.g. "scRNA Res-S22"). Skip
VDJ/CITE/pediatric/control samples **and the 8 ICU-control + 7 healthy-control
donors** that the naive `cohort=adult` filter would otherwise silently include.

**Step 2 — integrate to option-C schema.** Per-sample QC (same thresholds as the
existing integration), merge into one h5ad: `Sample_ID`=donor_id (e.g. Res-S22),
`Disease`=Sepsis, `Data_NO`=chongqing2025, `covid_severity`: mort28==1 → `fatal`,
else `severe` (survivor). Drop the 2 donors with missing mort28. Record the
endpoint as 28-day mortality in obs (`mort14/mort28/mort90` columns preserved).

**Step 3 — platelet annotation** with the SAME procedure as the current
integration (clustering + platelet marker module; report the platelet–monocyte
aggregate fraction). Output: platelet-subset h5ad + per-donor platelet counts.
(Pre-check on GSM8571191: 18.3% of cells carry platelet markers ≥2 counts.)

**Step 4 — embed the 5 FMs** (per-FM envs per `configs/hpcc.yaml`; platelet
subset estimated 20–80k cells → hours, not days) →
`02_EMBEDDINGS/option_c_chongqing/`.

**Step 5 — evaluation (patient-level, v2 rules; REPEATED-CV 10 seeds × 5-fold
everywhere — see the PROTOCOL CHANGE block in Job H2):**
- **within-sepsis mortality CV** on the 91 patients (17 fatal) — NEW primary
  result; repeated-CV (10 seeds × 5-fold patient-level StratifiedGroupKFold),
  LR+RF, per-seed OOF predictions emitted; report repeat-mean ± sd + per-seed
  paired pairwise tests.
- **COVID→sepsis transfer** (test = this cohort) and **sepsis→COVID transfer**
  (train includes this cohort), patient-level, n inline.
- **adaptation curve rerun** — sepsis training pool is now 103 patients
  (12 legacy + 91), test split patient-level; 10 seeds; note the protocol is
  otherwise unchanged, so pre/pos-expansion curves are directly comparable.
- PCA-200+LR baseline on raw counts for all of the above.

**Step 6 — report back:** per-donor platelet counts + aggregate fraction,
fatal-patient list with mort28 join check, within-sepsis table (5 FMs + PCA-200
+ null), both transfer directions, the new adaptation curve.

**Job S STATUS (2026-09-04 → 2026-09-05): ALL STEPS DONE, all 5 FMs + PCA-200 +
null, including the adaptation curve.** §五 report-back checklist (see
`HPCC_next_steps.md`) is fully closed for Job S.

- **Step 1 (download):** fixed a real GEO bug -- GSM8571120/121's donor
  identity is swapped between GEO's own sample-title metadata and the actual
  tar archive filenames (would have silently swapped two patients' expression
  data under the wrong mortality label); fixed via the tar's own filelist.txt
  as ground truth. Also corrected the cohort scope: the "108 donors/20 fatal"
  figure in `sepsis_arm_expansion.md` silently included 8 ICU-control + 7
  healthy-control donors under `cohort=adult`; **used `group=='Adult Sepsis'`
  only -> 93 donors, 2 missing mort28 dropped -> 91 donors, 17 fatal.**
- **Step 2-3 (integrate + platelet-annotate):** original annotation script not
  recoverable (same gap as Job F's Phase 1) -- used a standard, documented
  cluster+marker pipeline instead (thresholds in
  `code/integrate_chongqing_sepsis.py`'s docstring). 91 donors -> 626,471 cells
  merged -> **22,261 platelet cells (3.8%)**, all 91 donors represented (min 3,
  median 145, max 1549 platelet cells/donor). Per-donor counts:
  `results_patient/chongqing_platelet_counts_per_donor.csv`.
- **Step 4 (embed):** all 5 FMs targeted `02_EMBEDDINGS/option_c_chongqing/`.
  Hit 3 more real bugs (beyond the gene-ID issue already documented above this
  block): (a) UCE's own code string-concatenates output-dir+filename with no
  path separator, landing the file one directory up from where the wrapper
  looked -- twice, on two different attempts -- not a UCE failure, verified
  the embedding was valid both times before fixing the path; (b) scPRINT-2's
  CLI writes `<output>_embedded.h5ad` and stores the embedding under
  `scprint_emb`, not `X_scprint2` -- fixed via the same rename the existing
  `normalize_scprint2_outputs.py` does for the canonical file; (c) had to
  redo STATE/UCE/scGPT's first successful run because Leiden clustering isn't
  perfectly reproducible run-to-run (23,973 vs 22,261 cells across two
  identical-code runs) -- moved to `_stale_first_integration/`, resubmitted
  all 5 against the second run's output so every FM embeds the same cell set.
  All 5 done. scGPT initially queued ~12h on the congested `gpu` partition
  (4/9 A100 nodes down); resubmitted on `short_gpu` (already proven by the
  other 4 FMs) with generic `--gres=gpu:1` -- completed in 40s. See the
  "UPDATE (2026-09-05)" note below the Fix-4/5 block for detail.
- **Step 5 (within-sepsis CV, patient-level, repeated-CV, ALL 5 FMs):**
  91 patients/17 fatal, repeat-mean +/- sd AUC (LR, mean aggregation):

  | FM | repeat-mean±sd | seed-42 (secondary) |
  |---|---|---|
  | UCE | 0.5477±0.0376 | 0.5946 |
  | scPRINT-2 | 0.5234±0.0492 | 0.5676 |
  | Geneformer | 0.5151±0.0431 | 0.5421 |
  | STATE | 0.5112±0.0399 | 0.5564 |
  | scGPT | 0.4490±0.0523 | -- |
  | PCA-200 | **0.4284±0.0817** | 0.5397 |

  **All 5 FMs are at or near chance level.** scGPT and PCA-200 are BOTH
  significantly *worse* than every one of the other 4 FMs (all pairwise
  q≤0.0037 -- e.g. STATE−PCA-200 Δ=0.0828 q=0.0049, STATE−scGPT Δ=0.0622
  q=0.0037; scGPT vs PCA-200 itself: NOT significant, q=0.40, i.e. the two
  worst performers are statistically indistinguishable from each other). Among
  the top 4, the only individually-significant FM-vs-FM pair remains STATE vs
  UCE (Δ=−0.0365, q=0.0098, UCE modestly ahead). This is a genuine result, not
  a pipeline bug -- reported as-is, not explained away. The Timing gate note
  below (fatal patients drawn earlier than survivors, median 10 vs 18.5 d) is
  a plausible contributing confound worth citing in Discussion, not a reason
  to dismiss the null finding. Full tables: `results_patient_v2/within_sepsis/`
  (`donor_level_pairwise_repeated.csv` has all 15 pairwise tests).
- **Step 6: SUPERSEDED IN FULL** by the Fix-2/3 patient-level primary +
  severe-included tables in the Fix-1/2/3 STATUS block above (both the legacy
  12-patient and Chongqing 91-patient cross-disease tables there now include
  all 5 FMs). This cell-weighted/severe-excluded-only version (point AUC (LR):
  STATE 0.5802/0.6458, UCE 0.5653/0.6181, Geneformer 0.5303/0.5942, scPRINT-2
  0.5300/0.5190) is kept below only as historical record, not current:
  `results_patient_v2/cross_disease_chongqing/cross_disease_chongqing_20260904_202252.csv`.
- **Adaptation curve: DONE (2026-09-05).** Chongqing-expanded sepsis pool =
  103 patients (12 legacy + 91 Chongqing), v2 mapping (Task-4 label fix
  included), test split 34 patients (6 fatal/28 survived per
  `make_sepsis_holdout`'s 1/3 stratified split -- a different, larger holdout
  than within-sepsis CV's, since this is a fixed single split not repeated-CV),
  train pool up to 16,248 sepsis cells (69 patients) at the full grid point.
  sbatch 28148282 ran the 7 FM/method curves successfully (7h+ for STATE
  through TranscriptFormer) but crashed on the PCA-200+LR step:
  `ValueError: cannot reindex on an axis with duplicate labels` in
  `load_pca_raw()` -- the Chongqing raw-counts gene-symbol index has duplicate
  labels (multiple Ensembl IDs sharing one symbol), which broke a plain
  `pd.Series.reindex()`. **Fixed** in `run_adaptation_curve_chongqing.py` via
  `_first_occurrence_index()` (dict-based first-occurrence lookup instead of
  reindex); PCA-200 rerun via `rerun_pca200_and_compile_chongqing.py`
  (sbatch 28152662, 1h34m -- had to drop the mem request from 128G to 24G
  after the first submission hit `AssocGrpMemLimit`, the account's own group
  memory quota, saturated by an unrelated array job also running under this
  account, not general cluster congestion) then compiled all 8 methods'
  per-method CSVs from disk into `all_rows.csv`/`summary.csv` (no rerun of the
  7 already-completed methods).

  **point_auc_mean at zero-shot (n=0) / n=1000 / full sepsis pool, the 6
  methods WITH Chongqing embeddings** (STACK/TranscriptFormer below trained
  on the tiny legacy-only 2,745-2,751-cell pool the whole curve, NOT
  comparable to these 6):
  | Method | n=0 | n=1000 | full (n≈16,245-16,248) |
  |---|---|---|---|
  | Geneformer | 0.532 | 0.487 | 0.599 |
  | STATE | 0.530 | 0.428 | 0.582 |
  | UCE | 0.274 | 0.298 | 0.577 |
  | scGPT | 0.492 | 0.430 | 0.569 |
  | scPRINT-2 | 0.369 | 0.413 | 0.532 |
  | PCA-200+LR | 0.365 | 0.294 | 0.410 |

  **Every method with a Chongqing embedding is at or below chance across most
  of the grid, and the trajectory is often non-monotonic** (several methods
  DIP below their zero-shot AUC at intermediate n before a partial recovery
  at the full pool; PCA-200+LR never recovers above 0.41). This is the third
  independent confirmation of the same signal (within-sepsis CV null, weak
  cross-disease transfer, now a weak/sub-chance adaptation curve): adding
  Chongqing sepsis data to a COVID-trained classifier does not rescue
  cross-disease transfer under this protocol -- reported as-is, no
  explanation given. STACK (0.398→0.824→0.858) and TranscriptFormer
  (0.732→0.882→0.897) look strong but are NOT evidence against this, since
  they never saw the harder expanded pool -- their curve is identical in kind
  to the pre-expansion legacy-only result. Full per-seed rows and all 8
  methods' complete grids: `results_patient_v2/adaptation_curve/all_rows.csv`
  + `summary.csv` (local) / `05_RESULTS/revision_donor_level_v2/
  adaptation_curve_bootstrap_chongqing/` (HPCC). New code:
  `code/rerun_pca200_and_compile_chongqing.py` +
  `code/rerun_pca200_and_compile_chongqing.sh`.

**Comparability-gap fix (2026-09-06, code review):** the "STACK/TranscriptFormer
below trained on the tiny legacy-only pool, NOT comparable" caveat above was,
until this fix, prose-only -- `all_rows.csv`/`summary.csv` and the per-method
CSVs had no machine-readable marker of which test set a row belongs to, and
`make_sepsis_holdout()`'s legacy-only split composition was silently discarded
(`_` in the return unpacking) rather than saved. Fixed:
(1) `run_adaptation_curve_chongqing.py`'s `adaptation_curve()` now takes and
stamps a `test_set` column (`expanded_103pt` / `legacy_12pt`) into every
output row, and the legacy-only split is captured and written to
`split_documentation_legacy_only.txt` (**4 test patients: 1 fatal, 3
survived** -- identical composition to the pre-expansion v1 split, confirmed
independently before rerunning anything);
(2) the delivered `all_rows.csv`/`summary.csv` (produced before this fix
existed) were regenerated post-hoc from the already-complete per-method CSVs
with the `test_set` column added and `summarize()`'s groupby extended to
include it -- no curve refitting needed, this was a labeling gap, not a
computation error. The full-pool-point contrast now on record demonstrates
exactly why this matters: TranscriptFormer 0.897 / STACK 0.858
(`legacy_12pt`) vs. 0.41-0.60 for all 6 `expanded_103pt` methods at the same
`n_sepsis_req=-1` grid point -- a naive same-table comparison would have
looked like TranscriptFormer/STACK are the best methods by a wide margin,
when it's entirely a smaller/easier-test-set artifact.

**Timing gate (Task 5 answered 2026-09-03):** NOT an early-sepsis cohort — draws
are median 16 d from admission (range 3–121; 56% >14 d) and fatal patients were
drawn earlier than survivors (median 10 vs 18.5 d) — an outcome-associated
timing confound to report in the within-sepsis analysis (state + stratified
sensitivity, no explaining). `hosp_days` semantics (draw-day vs total LOS) needs
one supplement confirmation; dates file covers 93/108 donors — use the mortality
CSV as label authority, dates for timing only. Cohort-table sampling-window cell:
"median 16 d from admission (range 3–121)".

---

## Job J — GSE151263 (Jiang) exclusion rerun  [CPU, <1h]  (2026-09-07 user decision B)

**Reason:** GSE151263 has no per-patient outcome record in the source paper (only
literature-level mortality was reported); its 7 "survivor" labels (`Data_NO`
`302004data18`, 476 cells) were a study-design inference, not verifiable.
Decision: exclude entirely from the mortality task so every remaining label is
traceable to a real outcome.

**Scope note — the spec's own "reverse transfer 103→96" framing was wrong, flagged
not silently resolved:** `jobS_cross_disease.py` (which produces the headline
§3d Sepsis→COVID numbers, Chongqing-only 91-patient sepsis arm) never touches
Jiang cells at all — it's entirely unaffected by this exclusion. The script
actually affected is `cross_disease_fatal_vs_survived_patient_v2.py` (the
legacy-only sepsis arm, which does include Jiang): 12→5 patients. Both scripts
and `run_adaptation_curve_chongqing.py` (both its `expanded_103pt`→`expanded_96pt`
and `legacy_12pt`→`legacy_5pt` pools) were patched with a `Data_NO != '302004data18'`
filter; the adaptation curve's `test_set` column values were relabeled accordingly.

**STATUS (2026-09-08): DONE — legacy cross-disease rerun AND adaptation curve
rerun both complete and verified.** (Adaptation curve submission was briefly
interrupted by HPCC filesystem slowness on 2026-09-07; pre-Job-J output backed
up first, then cleanly resubmitted as sbatch 28187218 once the outage cleared
— completed 8h17m, exit 0.)

- Legacy reverse-transfer rerun: sbatch 28187143 (~9 min, via
  `run_jobH2_crossdesease.sh`'s existing wrapper after an interactive attempt hung
  on the login node). Output: `05_RESULTS/revision_donor_level_v2/cross_disease_patient_jobJ_20260907_154957.csv`
  (old non-`_jobJ` file preserved for before/after diffing).

**Sepsis→COVID (primary protocol), before (12-patient legacy pool) vs. after
(5-patient pool, Jiang excluded):**

| FM | Old LR | New LR | Old RF | New RF |
|---|---|---|---|---|
| STATE | 0.7230 | 0.7993 [0.718, 0.876] | 0.6396 | 0.6527 [0.555, 0.758] |
| UCE | — | 0.6534 [0.536, 0.755] | — | 0.6313 [0.527, 0.732] |
| scPRINT-2 | — | 0.6954 [0.608, 0.797] | — | 0.7014 [0.624, 0.793] |
| Geneformer | — | 0.6677 [0.574, 0.766] | — | 0.5603 [0.466, 0.659] |
| scGPT | — | 0.6608 [0.536, 0.765] | — | 0.5763 [0.461, 0.694] |

(Old UCE/scPRINT-2/Geneformer/scGPT values not pulled this pass — diff against
the preserved non-`_jobJ` CSV directly if needed.) COVID→Sepsis direction
(test=5 patients, extremely small) shows degenerate 1.000/0.667 AUCs for
several FMs — expected small-n instability from a 5-patient test set,
consistent with this project's established pattern; report as-is.

**Adaptation curve, Jiang-excluded — endpoint (n=all sepsis training cells)
vs. pre-Job-J:**

| Method | test_set | n_test | New endpoint (point/boot-mean) | Old endpoint |
|---|---|---|---|---|
| Geneformer | expanded_96pt | 32 (6 fatal) | 0.616 / 0.587 | 0.599 |
| STATE | expanded_96pt | 32 | 0.574 / 0.559 | 0.582 |
| scPRINT-2 | expanded_96pt | 32 | 0.572 / 0.569 | 0.532 |
| UCE | expanded_96pt | 32 | 0.559 / 0.540 | 0.577 |
| scGPT | expanded_96pt | 32 | 0.532 / 0.518 | 0.569 |
| PCA-200+LR | expanded_96pt | 32 | **0.404 / 0.416** | 0.410 |
| TranscriptFormer | legacy_5pt | **2** | 0.765 / 0.765 | 0.897 |
| STACK | legacy_5pt | **2** | 0.729 / 0.729 | 0.858 |

**Null persists, if anything slightly weaker across the board** (all 6
expanded_96pt methods land 0.40–0.62, same order of magnitude as the pre-Job-J
0.41–0.60 range, no rescue). **Data-quality note, not a bug:** `legacy_5pt`'s
test set is now only 2 patients (1 fatal `NS2` + 1 survived `S1` —
`split_documentation_legacy_only.txt`), so its bootstrap CI is degenerate
(ci_lo=ci_hi=point estimate identically) — with 2 test items, patient-level
bootstrap resampling has only one non-degenerate outcome available. TF/STACK's
endpoint dropped further (0.897→0.765, 0.858→0.729) on this now-even-smaller
easy test set — reinforces, doesn't weaken, the "not comparable, don't put in
main text" verdict on these two rows. §3d's Chongqing-based forward/reverse
transfer numbers and the `expanded_96pt`-family adaptation-curve numbers
supersede the pre-Job-J `expanded_103pt` ones; full per-seed data in
`results_patient_v2/adaptation_curve_jobJ/`. Pre-Job-J reference preserved at
`05_RESULTS/revision_donor_level_v2/adaptation_curve_bootstrap_chongqing_pre_jobJ/`
on HPCC.

**Job J is now fully complete — both reverse-transfer and adaptation-curve
reruns done.** Note for whoever picks up Jobs K-O next: the "Jobs K-O" section
below records that "Job P" (purge legacy 4-patient adaptation-curve rows) was
deliberately held pending Job J landing — it's unblocked now.

---

## Not HPCC jobs (handled separately, listed for completeness)

- P1-D cohort-characteristics table with sampling time (desk research on the 11
  source studies; skeleton from `docs/cohort-sources.md`).
- P1-E platelet/megakaryocyte fraction in each FM's pretraining corpus
  (literature + CELLxGENE census query).
- Task 5: GSE279452 sampling-time post-diagnosis (desk research on the paper's
  supplement; feeds the Job S cohort-table cell and Step-5 interpretation gate).

## Job P — purge legacy 4/5-patient adaptation-curve numbers (2026-09-08) -- STATUS: DONE

**Directory-name correction before executing**: the spec's literal target,
`results_patient_v2/adaptation_curve/`, was stale -- Job J's rerun writes to
`adaptation_curve_bootstrap_chongqing/` on HPCC (mirrored locally as
`adaptation_curve_chongqing/`), and had already relabeled `test_set` values
from `expanded_103pt`/`legacy_12pt` to `expanded_96pt`/`legacy_5pt` after the
Jiang exclusion. **There was no `legacy_12pt` value left anywhere to delete**
-- confirmed directly against the live data before touching anything. Also:
the specific old numbers named in the spec (TF 0.897 / STACK 0.858) no
longer exist either -- Job J's rerun already produced different values for
the new 5-patient legacy test set (**TF 0.765, STACK 0.729** at the n=all
grid point). Purged the CURRENT legacy_5pt rows regardless, since Job P's
actual intent (no small-n legacy-test-set numbers in the quotable record)
applies to whatever is currently there, not just the specific stale values
named when the spec was written.

**Backup**: full pre-purge directory copied to
`05_RESULTS/attic/adaptation_curve_bootstrap_chongqing_pre_jobP_backup/`
(HPCC only, not synced to git -- it's a full duplicate of already-tracked
data, kept only as an HPCC-side undo point).

**Purge, on the live HPCC directory** (`adaptation_curve_bootstrap_chongqing/`):
`all_rows.csv` 488->366 rows, `summary.csv` 56->42 rows (both: dropped every
`test_set==legacy_5pt` row, kept `expanded_96pt`). `STACK.csv` and
`TranscriptFormer.csv` deleted entirely -- both files were 100% legacy_5pt
rows (neither method has a Chongqing embedding, so they had no
`expanded_96pt` data to keep). `split_documentation_legacy_only.txt`
deleted (documented only the now-purged split). Verified clean afterward:
no `legacy_5pt`/`legacy_12pt`/`TranscriptFormer`/`STACK` string anywhere in
the remaining CSVs.

**Local git mirror**: `results_patient_v2/adaptation_curve_chongqing/`
refreshed to the purged state (stale local copy still had `legacy_12pt`
rows from before Job J's rerun -- replaced entirely, not merged).
`results_patient_v2/adaptation_curve_jobJ/` -- a redundant duplicate
snapshot of the same post-Job-J/pre-purge state, serving no distinct
historical purpose the way v1-vs-v2 splits do elsewhere in this
project -- deleted rather than kept as clutter with the same legacy numbers
this job exists to remove. The true v1 (pre-Chongqing, 12-patient-only)
curve at `results_patient_v2/adaptation_curve/` is untouched -- that is a
different, legitimately-historical job (Job H's original run), not a
duplicate of anything Job P targets.

**TF/STACK numbers now exist ONLY in**: (a) this runbook's own historical
STATUS prose (Job S's original completion report, kept as a documented
record of what was found and fixed, same convention as every other
correction in this file), and (b) the HPCC-only attic backup. Neither is
part of the actively-quotable results record.

---

## Job J-supplement — reverse-transfer 96-pool completion (2026-09-08)

Per `HPCC_next_steps.md` §11. Job J's own reverse-transfer number (in
`cross_disease_patient_jobJ_20260907_154957.csv`) turned out to be the
5-patient GSE167363-only pool, not the intended 96-patient pool (91
GSE279452 + 5 GSE167363, Jiang excluded) -- caught before it became the
manuscript's canonical reverse-transfer value.

**STATUS: DONE.** New script `code/cross_disease_96pool.py` (sbatch 28205555,
15m26s, exit 0) combines `jobS_cross_disease.py`'s Chongqing-loading logic
with `cross_disease_fatal_vs_survived_patient_v2.py`'s Jiang-excluded
legacy-Qiu-loading logic. Verified pool: Qiu 2,946 cells/5 patients +
Chongqing 22,261 cells/91 patients = 25,207 cells/96 patients (exact match).

Sepsis->COVID (96-pool, 240 patients/25 fatal test set), primary protocol LR/RF:
STATE 0.7209/0.7113, Geneformer 0.7118/0.6358, UCE 0.6928/0.6074, scGPT
0.6417/0.6221, scPRINT-2 0.6186/0.5167. vs. the old 91-pool value (STATE-LR
0.689): diff = 0.032 >= the 0.02 threshold -> triggered a manuscript text
update, not just a table swap. STATE and Geneformer's LR point estimates are
now only 0.009 apart -- `P0_onepager.md` section 3d rewritten to state both
rather than "STATE alone," CIs noted as fully overlapping.

Legacy 5-patient pool's adaptation-curve endpoint values confirmed unchanged
(STACK 0.7289, TranscriptFormer 0.7646, from the pre-Job-P attic backup --
only these two lack a Chongqing embedding, the 5 primary FMs+PCA-200 all ran
on the full expanded pool instead, no separate legacy-only run exists for them).

Output: `05_RESULTS/revision_donor_level_v2/cross_disease_patient_96pool_20260908_142954.{csv,json}`,
synced to `results_patient_v2/cross_disease/` locally.

---

## Jobs K-O — external pre-submission review follow-ups (2026-09-07)

Per `reports/manuscript/revision/HPCC_next_steps.md` §八. Job P (purge legacy
4-patient adaptation-curve rows) deliberately NOT run here -- held until Job J
is confirmed fully landed (concurrent file access conflict).

**Job K (QC-metrics-only baseline, highest priority) -- STATUS: DONE, decision
rule NOT triggered.** LR on log1p(nCount_RNA)/log1p(nFeature_RNA)/pct_mito per
cell, identical repeated-CV (10x5, patient-grouped) to the main analysis, same
v2 filters/patient-map/severity-fix. **AUC = 0.5282 +/- 0.0260 (10 seeds)** --
essentially chance. The main FM result is NOT a QC/technical artifact.
n=252 patients (240 COVID + 12 legacy sepsis pool, same scope as the null
model), 27 fatal. Per-patient QC by disease x outcome (median of medians):
COVID survived 1105 nUMI/474 nGene/2.33% mito; COVID fatal 1554/652/2.79%;
sepsis survived 804/365/3.82%; sepsis fatal 1945/611/5.21% -- fatal patients
have systematically higher nUMI/mito than survivors in both diseases, which
is exactly why a QC baseline could in principle have picked up signal, and
notably still didn't reach 0.7. Platelet-calling rules documented: legacy =
2024 IJMS cluster+marker; GSE279452 = cluster+marker 12-gene panel
(`integrate_chongqing_sepsis.py`); Fix-4 aggregate co-expression report
already on record at `scripts/31_donor_level_revision/chongqing_platelet_counts_per_donor.csv`.
Code: `code/jobK_qc_baseline.py`. Results:
`results_patient_v2/qc_baseline/{qc_baseline_auc,qc_baseline_per_seed,patient_qc_distributions}.csv`.

**Job L (definitive patient-level paired inference) -- STATUS: submitted
(sbatch 28193089), running at last check.** All 10 seeds' oof_patient
files, within-COVID-only (COVID-19 disease patients per the v2 map), paired
patient bootstrap B=10,000 for every FM-PCA200 and STATE-Geneformer, plus a
seed-ensemble version, Hanley-McNeil SEs, and pairwise score correlations.
This determines the actual headline wording (seed-42 pilot: STATE-PCA
p=0.059, Geneformer-PCA p=0.002 -- not yet confirmed on the full 10-seed
version at the time of this STATUS write-up). Code: `code/jobL_paired_inference.py`.
Results land in `results_patient_v2/paired_inference/` once complete --
**pull and report before treating any headline significance wording as
final.**

**Job M (GSE279452 age/site stratification) -- STATUS: PARTIAL, blocked on
paywalled source.** Confirmed via GEO's own plain-text series record
(`acc.cgi?...&form=text`) that GSE279452 explicitly separates Adult and
Pediatric sepsis cohorts as distinct groups (133 total PBMC samples across
both) -- our 91-donor scope, filtered to `group=='Adult Sepsis'`, is
therefore drawing from the correct age-stratified group at the study-design
level, not an unstratified pool. Could NOT obtain per-patient exact ages:
the source paper (Nature Immunology 2025, "Single-cell multi-omic landscape
reveals anatomical-specific immune features in adult and pediatric sepsis")
is paywalled, no free PMC mirror found, and GEO's own per-sample
characteristics don't expose age directly (checked the brief record; a full
SOFT/MINiML parse of all 91 GSM records was not attempted given time
budget). **Age distribution, age-stratified fatal counts, and explicit
per-patient age>=18 confirmation remain open** -- the group-label evidence is
reassuring but not a substitute for the actual per-patient numbers this job
asked for. Site (single-center, Chongqing) and chemistry (10x) can be
confirmed from the existing `gse279452_donor_sampling_dates.csv` scout file
already in the repo.

**Job N (LOBO on v2, secondary) -- STATUS: DONE.** Per-cohort (Data_NO)
leave-one-out within-COVID, restricted to cohorts with >=1 fatal and >=1
survivor patient (5 of the pool's cohorts qualify); LR only. **Per-cohort
mean AUC runs notably HIGHER than the pooled repeated-CV estimate for every
FM** -- scPRINT-2 0.875+/-0.135, Geneformer 0.869+/-0.141, STATE 0.859+/-0.132,
UCE 0.816+/-0.152, scGPT 0.779+/-0.123 (vs. pooled 0.740/0.786/0.778/0.661/0.720
respectively) -- a reassuring cross-cohort generalization signal, though n=5
cohorts per FM makes the sigma wide. Code: `code/jobN_lobo_v2.py`. Results:
`results_patient_v2/lobo_v2/{lobo_per_cohort_auc,lobo_fm_summary}.csv`.

**Job O (monocyte patient-matched comparison) -- STATUS: DONE. Earlier
"blocked, no crosswalk" assessment in this file was WRONG -- corrected by
the coordinator, who caught that Ren's actual platelet-side Data_NO is
`302004data15` (patient_ids `S-Mxxx`/`S-Sxxx`), not `302005data01` (that
Data_NO is Bernardes, numeric ids like `001`) -- Ren, Stephenson
(`302004data17`, e.g. `MH9143327`), and Qiu sepsis (`302005data03`, e.g.
`NS1`/`S1`) all use IDENTICAL patient_id strings on the platelet and
monocyte sides already; no crosswalk was ever needed.** Intersection:
**148 of 148 monocyte patients also have platelet data** (53 Ren/9 fatal,
90 Stephenson/8 fatal, 5 Qiu-sepsis/2 fatal) -- essentially the entire
monocyte cohort, not a small subset. On this exact patient-matched,
cell-count-matched comparison (seed-42 single split, all-disease pool since
that's this comparison's natural scope): **the official mono<platelet
inversion does NOT hold uniformly** -- STATE (mono 0.691 vs plt 0.660), UCE
(0.669 vs 0.539), scGPT (0.698 vs 0.575), and PCA-200 (0.686 vs 0.667) all
favor monocyte here; only Geneformer (0.683 vs 0.750) and scPRINT-2 (0.596
vs 0.650) favor platelet. **Reporting this plainly as a real complication,
not explaining it away**: this comparison differs from the official
repeated-CV/COVID-only framing in three ways (single seed-42 split vs.
10-seed repeat; all-disease pool vs. COVID-only; n=148 patients on both
sides vs. n=240 platelet-side). Code: `code/jobO_monocyte_matched.py`. Results:
`results_patient_v2/monocyte_matched/{intersection_patients,cell_count_matching,matched_comparison_auc}.csv`.

**Follow-up -- root cause of the flip, INVESTIGATED (not left untested).**
Ran the same matched comparison across all 4 combinations of {seed-42 vs.
10-seed ensemble score} x {all-disease n=148 vs. COVID-only n=143 pool}
(`code/jobO_investigate_flip.py`, results
`results_patient_v2/monocyte_matched/flip_investigation.csv`). **The flip was
primarily seed-42 single-split noise** -- exactly the same fold-instability
mechanism that motivated switching to repeated-CV in the first place, now
demonstrated again on a completely different comparison. Ensembling over 10
seeds flips STATE (+0.031 to +0.041 seed42 -> **-0.040/-0.032 ensemble**),
PCA-200 (+0.019 to +0.032 -> **-0.076/-0.061**), and scGPT (+0.123 to +0.138
-> ~0.00) back to matching or being consistent with the official mono<
platelet direction; the COVID-only vs. all-disease restriction changes
magnitudes only slightly and never flips a sign, so patient-pool composition
is a minor secondary factor, not the driver. Geneformer stayed
platelet-favoring in all 4 configs (-0.067 to -0.073). **UCE is the one
genuine, robust exception**: monocyte-favoring in all 4 configs (+0.060 to
+0.155), including both ensemble configs where every other flipped FM
reverted -- this is not noise, UCE really does perform better on monocytes
than platelets even under a fair, patient-matched, seed-averaged comparison.
scPRINT-2 shows a small (+0.009 to +0.012) sign reversal under ensembling,
likely noise-level given the tiny magnitude, worth a footnote not a claim.
**Manuscript implication**: the "task/cohort property, not FM property"
framing for the mono<platelet inversion holds robustly for 5 of 6 methods
once evaluated properly (matched patients + repeated-CV/ensemble scoring),
but UCE is a genuine, real exception that should be called out by name
rather than folded into the aggregate claim.

**Job M (GSE279452 age/site) -- STATUS: DONE.** The 10MB single-fetch cap
was a limitation of the coordinator's own fetch tool, not of GEO -- HPCC's
login node has direct internet egress (verified) and pulled the FULL
266-sample series record in one shot (`acc.cgi?acc=GSE279452&targ=gsm&
view=brief&form=text`, 891KB, no limit hit). Parsed all 266 samples' title/
age/group/gender via regex on `Sample_characteristics_ch1` lines. Joined
against the 93 `cohort==adult & scRNAseq==1 & group=='Adult Sepsis'`
donor_ids from the existing scout file (case/format-normalized to handle
title variants like `RES-S118` vs `Res-S118`, `Skin-152` vs `Skin-S152`):
**92 of 93 matched (1 residual title-format miss, `Skin-S152`, not chased
further).**

- **Age: n=92, median=67, range=19-93. Zero pediatric contamination**
  (min age 19, well clear of 18) -- the `group=='Adult Sepsis'` filter is
  confirmed correct at the per-patient level, not just the study-design
  level.
- **Age-stratified fatal counts: <40 n=6/0 fatal; 40-60 n=30/7 fatal;
  >60 n=56/10 fatal.** 7+10=17 fatal sums exactly to the established
  91-patient/17-fatal figure used throughout this revision -- a clean
  internal consistency check.
- Site (single-center, Chongqing) and chemistry (10x) already confirmed
  via the existing `gse279452_donor_sampling_dates.csv` scout file.

Code: none saved as a standalone script (ad-hoc parse of
`/tmp/gse279452_gsm_list.txt` on HPCC + local join) -- results:
`results_patient_v2/monocyte_matched/{gse279452_age_table_raw,
gse279452_age_final}.csv`.

---

## Job Q — paired bootstrap completion (2026-09-08) -- STATUS: DONE

Same seed-ensemble patient-level paired bootstrap method as Job L (B=10,000,
average each patient's 10 per-seed scores, one bootstrap on the ensembled
scores). Scripts: `code/jobQ_paired_inference_full.py` (modes: within_covid,
within_sepsis, monocyte), `code/jobQ_qc_covid_only.py`. Total runtime well
under the spec's 30-min estimate (all 4 sub-tasks ran in a few minutes each,
interactively, no sbatch needed).

**1. Within-COVID, 15 pairs** -> `paired_inference_full/within_covid_15pairs.csv`
(240 patients). Cross-checks Job L's original 6 pairs closely (e.g.
STATE-PCA200 p=0.037 here vs. 0.038 in Job L; tiny RNG-seed difference,
20260908 vs 20260907 -- expected, not a discrepancy).

**2. Within-sepsis, 15 pairs** -> `paired_inference_full/within_sepsis_15pairs.csv`
(91 patients). **Real finding, not a labeling change**: under honest
seed-ensemble inference, NO pair reaches significance after accounting for
15 comparisons -- the sole nominal (uncorrected) exception is UCE-scGPT
(p=0.047). In particular NONE of the FM-vs-PCA-200 deltas are significant
(p=0.079-0.575). This directly contradicts the retired per-seed Wilcoxon
result ("PCA significantly worse than 4/5 FMs, q=0.0037") -- confirms that
method was anti-conservative, exactly as suspected. Figure 3 updated with
this note.

**3. Monocyte, 5 FM-PCA-200 pairs** -> `paired_inference_full/monocyte_vs_pca_5pairs.csv`
(143 patients, COVID-only cut = Ren+Stephenson, matching Figure 4's primary
framing -- verified `source_study != 'Qiu_sepsis_GSE167363'` gives exactly
53+90=143). **Real finding**: under honest inference, NO FM-vs-PCA-200
comparison is significant in-distribution (STATE p=0.23, UCE p=0.80,
Geneformer p=0.57, scPRINT-2 p=0.68, scGPT p=0.94) -- reverses the retired
per-seed result (STATE p=0.002, scGPT p=0.037, both "significant"). Figure 4
panel B rebuilt on this data with per-point p-value annotations.

**4. QC baseline, COVID-only variant** -> `qc_baseline/qc_baseline_auc_covid_only.csv`
(240 patients/25 fatal). **AUC = 0.3762 +/- 0.0429** -- notably BELOW chance,
and further from the 0.7 decision threshold than the pooled 252-patient
cohort's 0.528. Reinforces (does not merely repeat) the "main result is not
a QC artifact" conclusion: within COVID alone, QC metrics carry no
outcome-predictive signal at all, let alone enough to explain the FM/PCA
performance.

**Downstream figure updates applied** (per figure_generation_spec.md's
"待 Job Q 落地后按 bootstrap p 重标" instruction for each): Figure 3 (caption
note, no stars added since nothing survives correction), Figure 4 panel B
(rebuilt from Job Q data, per-point p-values, title changed to "no pair
significant"), FigureS2 (significance asterisks restored from Job Q's
within-COVID 15-pair p-values, magnitudes unchanged from the existing
repeated-CV pairwise CSV per spec). All three regenerated and verified
(assertions in each script guard against silently stale significance
claims if the underlying data changes again). Local matplotlib venv used
for regeneration (`/tmp/figvenv`, not part of the repo) since this session
had no local Python data-science stack -- HPCC was only used for the actual
bootstrap computation.

---

## Job R — seed-ensemble winner-corrected max-statistic test + fragility + consensus (2026-09-08) -- STATUS: DONE

External round-3 review decision item: the primary within-COVID FM-vs-PCA-200
contrast must be a winner-corrected max-statistic test (does ANY of the 5 FMs
beat PCA-200, corrected for testing 5 models at once), not 5 separate per-FM
tests. New script `code/jobR_winner_test_ensemble.py` (adapted from
`jobL_paired_inference.py`, same seed-ensemble construction, B=10,000,
RNG seed fixed at 20260908, sanity-gated on n_patients==240/n_fatal==25 with
no output on failure). sbatch 28203972, 16m54s, exit 0, gate passed.

**familywise P = 0.0045** (local seed-42 pilot: 0.0003 — same order of
magnitude, <0.01 → per the pre-agreed escalation rule, no manuscript text
changes needed, just fill in the placeholder `tab:winner_ensemble` table).
obs max delta = 0.0644 [0.0145, 0.1283] (Geneformer is the argmax).

Per-model vs. PCA (uncorrected point estimates, reference only): STATE
+0.0523 (p=0.0328), Geneformer +0.0644 (p=0.0264), scPRINT-2 +0.0138 (n.s.),
scGPT −0.0032 (n.s.), UCE −0.0659 (n.s., negative direction).

**Fragility (25 leave-one-fatal-out refits):** Geneformer−PCA 6/25 lose
significance (pilot: 0/25 — ensemble slightly less bulletproof but still
19/25 robust); STATE−PCA 7/25 lose significance (pilot: 17/25 — **ensemble
is more robust than the pilot suggested**, not less).

**Consensus (5-FM mean) AUC = 0.7678, beats PCA-200 by +0.0255** — reverses
the pilot's finding that consensus underperformed PCA (0.6995, ensembling
"dilutes" signal); the dilution story does not hold in the official version,
reported as-is. Consensus still trails the best single FM (Geneformer 0.8067,
Δ=−0.0389).

Correlation matrix (Pearson, 6x6): 0.707–0.904 (pilot: 0.61–0.88 — range
shifted up, same order of magnitude).

Outputs: `results_patient_v2/paired_inference_full/{winner_test_ensemble,
fragility_ensemble,consensus_ensemble,score_correlation_ensemble}.csv`,
synced to local mirror.
