# Job 2 preparation packet — v1 (2026-10-11)

**SUPERSEDED by `job2_preparation_packet_v2.md` (2026-10-11, same day).** A
review (`codex_output/20261011_124253.md`, R4) found this v1 incomplete: the
v2 wrapper fix had not actually been deployed anywhere yet, resource/sbatch
details were missing, it used the wrong manifest field name
(`source_hashes` instead of the real `source_file_sha256`), and acceptance
criteria did not cover SLURM exit status or per-lineage/embedding
completeness. v2 fixes all of these, including actually deploying
`step2_run_diagnostics_v2.sh` to a new staging path and verifying its
sha256. Kept for the audit trail, not for use.

Read-only inspection, no SLURM submission. Produced in response to
`codex_output/20261011_113512.md` §4 item 5. This packet does not authorize
`sbatch`; it exists so the one remaining authorization decision (§5) can be
made from real, checked facts rather than a prior round's restated plan.

## 1. What this job is

`run_batch_diagnostics.py`, invoked via a wrapper SLURM script, run once on
Job 1's completed output (`provisional_diagnostic_input_v3/input.h5ad`,
24,841 cells × 11,709 genes, 933 MB). Produces provisional cohort-level
batch/outcome diagnostics (metrics, optionally Harmony-corrected embeddings)
for exploratory review — not a final correction method, not a canonical
grouping decision, not a manuscript result.

## 2. Exact deployment identity (checked 2026-10-11, this round)

| Item | Value |
|---|---|
| Deployment root | `/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_aea9a61` |
| Deployment `git rev-parse HEAD` | `aea9a6197539ddd5086e46ddd34561617ad250d7` (matches the directory name) |
| Deployment `git status --porcelain` | empty — clean, no uncommitted changes |
| Deployment HEAD commit | `aea9a61 Fix premature extracted status in preparer; wrapper checks clean tree too` |
| `run_batch_diagnostics.py` sha256 (in deployment) | `2dadf4f394d8fa198ac547d2129b7710e33c8668ba2d0b39ececa659d67722fe` |
| `step2_run_diagnostics.sh` (v1) sha256 (in deployment) | `a1b283ca740d0dfc791d8e18b145384a0ff39b18688adae8707565a6673b29e6` |

This is the SAME frozen deployment Job 1 ran from — not touched by this
round or any prior round. It remains distinct from, and does not affect,
the working-tree repo (`302006_updategit`, now at commit `fde3eea` after
this round's F2/F4/F5 fixes) or any `_testrun_*` scratch deployment.

## 3. The one real wrapper defect found and fixed in a NEW, unDEPLOYED file

`step2_run_diagnostics.sh` (v1, the version actually in the frozen
deployment above) line 60:

```sh
POP_STAGE=$($PYTHON -c "import anndata; print(anndata.read_h5ad('$H5AD', backed='r').uns['population_stage'])")
```

A 2026-10-10 review (R1, `codex_output/20261010_191146.md`) read the
installed `anndata` source directly and confirmed `backed="r"` still
materializes every sparse `layers` entry in full — so this line, used only
to peek at one tiny `uns` string, was silently forcing a full read of the
whole h5ad's layers before `run_batch_diagnostics.py` itself even starts.

**Checked this round**: for the real 933 MB `input.h5ad`, this is a
correctness/consistency defect, not a demonstrated OOM risk under the 1-hour/
16 GB allocation — the file is small. Fixed in a new file,
`scripts/33_batch_effect_audit/code/slurm/step2_run_diagnostics_v2.sh`
(working tree only, **not deployed**, not run), replacing the line above with
a direct `h5py` read of the `uns` scalar (the same pattern
`validate_job1_output.py`'s `read_uns_scalar()` already uses). Verified on
HPCC this round that the h5py read and the `anndata.read_h5ad(...,
backed="r")` read return the IDENTICAL string
(`'pre-Gate-2 PROVISIONAL, not final (see docs/current-status.md)'`) — not
assumed equivalent.

Deploying `step2_run_diagnostics_v2.sh` into a new, separately-reviewed
deployment (its own commit/checksum, its own `EXPECTED_REVISION` update if
the deployment root changes) is a separate step from this packet, not done
here. `run_batch_diagnostics.py` itself has never used `backed="r"` — it
reads the full h5ad into memory at its own line 549, which is the intended,
necessary behavior for the PCA/Harmony computation it performs, not a bug.

## 4. Inputs, identity, and parameters

| Item | Value |
|---|---|
| `--h5ad` | `/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/input.h5ad` |
| `input.h5ad` sha256 | `e2da0e4f44350a1ecb28c1894308812f9234c331420672a52177f9445305d6ba` (from this round's `validate_job1_output.py` run; unchanged since) |
| `--out-dir` | `/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/diagnostic_run` |
| `--out-dir` current state | **does not exist** (checked this round) — first run, no overwrite |
| `--population-stage` | read from the h5ad itself at run time, not retyped (see §3) |
| `--stage-id` | `pre_gate2_decontx_provisional` (matches Job 1's `validate_job1_output.py` `stage_id_matches_expected: true` this round) |
| `--code-revision` | resolved from `DEPLOY_ROOT`'s own `git rev-parse HEAD` at run time, hard-fails if unresolvable or dirty (§2) |
| `--expression-source` | `decontx_counts` (matches Job 1's `expression_source_matches_expected: true` this round) |
| `--batch-col` / `--patient-col` / `--lineage-col` / `--outcome-col` | `cohort` / `patient_id` / `lineage` / `outcome` — grepped directly; **no `cv_group` anywhere in this call path**, confirmed again this round |
| `--run-harmony --seed 0` | exploratory diagnostic embedding only; not a final correction method |

## 5. The CV0201/CV0326 identity conflict's concrete relevance to THIS job

Checked this round with the new `stephenson_disputed_donor_impact.py` tool
(real run, `stephenson_disputed_donor_impact_20261011_123120/`):

- CV0201 and CV0326's cells ARE present in `input.h5ad` (40 and 33 cells
  respectively, Platelet + Monocyte_Myeloid) — **not zero**, so this is
  recorded explicitly rather than assumed away.
- Job 2's call path uses `cohort`/`patient_id`/`lineage`/`outcome` only —
  none of which depend on BGCV06/BGCV08/BGCV13 pool membership. The disputed
  fact (which physical pool a donor's cells came from) is not read by this
  job at all.
- Each disputed donor's own identity label (`CV0201`, `CV0326`) is itself
  independently confirmed consistent across GEX's own `patient_id` field,
  `sample_id`-derived donor, BCR, and (where observed) TCR — 0 disagreements
  found this round and in the prior round. The pool-membership dispute does
  not call either donor's own cell-level identity into question.
- **Conclusion**: this job's result would not change under either proposed
  resolution of the BGCV06/BGCV13 conflict. This is a fact about THIS job's
  input columns, not a general claim that the conflict never matters
  anywhere (e.g. it would matter to a batch/pool-column-based analysis,
  which this job does not perform).

## 6. Acceptance criteria for a run, if authorized

Before accepting any output from a real run:

1. `--out-dir` did not exist beforehand (checked above) and the run's own
   non-empty-dir guard (`run_batch_diagnostics.py` line ~491) was not
   bypassed.
2. The wrapper's own `code_revision`/working-tree-cleanliness checks passed
   (visible in the SLURM log) — not merely that the job exited 0.
3. `manifest.json` (written by the run itself) records the real
   `source_hashes`/`code_revision`/`population_stage` actually used, matched
   against §2 and §4 above.
4. The log path (`06_LOGS/stage33_diag_run_<jobid>.log`) is read in full, not
   sampled, before calling the run successful.
5. Harmony-corrected embeddings (if produced) are treated as exploratory
   diagnostics only — not as a correction decision, not as a basis for any
   manuscript claim, per `AGENTS.md`'s human-in-the-loop requirement for
   scientific decisions.

## 7. What this packet does NOT do

Does not submit `sbatch`. Does not deploy `step2_run_diagnostics_v2.sh`
anywhere. Does not modify `302006_updategit_deploy_aea9a61` or any other
existing deployment. Does not change `canonical_donor_labels.csv`, QC, or
DecontX. Does not resolve the CV0201/CV0326 identity conflict (§5 states its
relevance to this specific job only). The decision to run Job 2 remains the
user's, per `AGENTS.md`.
