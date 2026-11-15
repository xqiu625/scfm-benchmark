# Job 2 preparation packet — v2 (2026-10-11)

Supersedes `job2_preparation_packet_v1.md`, which a review
(`codex_output/20261011_124253.md`, R4) found incomplete: the v2 wrapper fix
was not yet actually deployed, resource/sbatch details were missing, the
wrong manifest field names were used, and acceptance criteria did not cover
SLURM exit status or per-lineage/embedding completeness. All fixed below.
Read-only + one new file copy (the wrapper, into a NEW staging path, not
into the frozen application); no SLURM submission.

## 1. What is new in this version

- `step2_run_diagnostics_v2.sh` is now actually **deployed** to
  `/bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/2ff4859/step2_run_diagnostics_v2.sh`
  (a NEW staging directory — checked it did not exist before copying).
  Post-copy sha256 verified identical to `git show
  2ff4859:scripts/33_batch_effect_audit/code/slurm/step2_run_diagnostics_v2.sh`
  (both `e4a7734470f49409f45bd68330f6502096eac0f734c1e8007d6283bfcf07d3ea`) —
  the file's bytes are unchanged since that commit, so `2ff4859` is the real
  source commit, not a placeholder.
- Verified this round, without executing the full diagnostic script:
  `bash -n` on the deployed file reports no syntax errors; the file's own
  `population_stage` h5py read (its only change from v1) was independently
  re-run against the real `input.h5ad` and still returns the identical
  string anndata's `backed="r"` read returns.
- `302006_updategit_deploy_aea9a61` (the frozen application the wrapper
  `cd`s into) re-checked this round: `git status --porcelain` is still
  empty. The wrapper copy above was placed in `stage33_wrappers/`, NOT
  inside `302006_updategit_deploy_aea9a61` — that deployment was not
  touched.
- Real field names corrected below (`source_file_sha256`,
  `code_revision_verified`, top-level `status` ∈
  `{"complete","complete_with_partial_results","failed"}`) — v1 used the
  wrong name `source_hashes`.

## 2. Exact deployment identity

| Item | Value |
|---|---|
| Application deployment root | `/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_aea9a61` |
| Application `git rev-parse HEAD` | `aea9a6197539ddd5086e46ddd34561617ad250d7` |
| Application `git status --porcelain` | empty (re-checked this round) |
| Wrapper staging path | `/bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/2ff4859/step2_run_diagnostics_v2.sh` |
| Wrapper source commit | `2ff4859` (bytes unchanged since that commit — verified by sha256, not asserted) |
| Wrapper sha256 (deployed AND git blob, identical) | `e4a7734470f49409f45bd68330f6502096eac0f734c1e8007d6283bfcf07d3ea` |
| `run_batch_diagnostics.py` sha256 (in the `aea9a61` application) | `2dadf4f394d8fa198ac547d2129b7710e33c8668ba2d0b39ececa659d67722fe` |

## 3. Resources and exact sbatch command

| Resource | Value |
|---|---|
| Partition | `epyc` |
| Nodes | 1 |
| Tasks | 1 |
| CPUs per task | 4 |
| Memory | 16 GB |
| Walltime | 1:00:00 |
| Log path | `/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stage33_diag_run_%j.log` (stdout+stderr combined) |

Exact command, if authorized:

```sh
sbatch /bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/2ff4859/step2_run_diagnostics_v2.sh
```

The wrapper itself invokes (no flags added or removed from the committed
file):

```sh
python scripts/33_batch_effect_audit/code/run_batch_diagnostics.py \
  --h5ad /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/input.h5ad \
  --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/diagnostic_run \
  --population-stage "<read from the h5ad at run time via h5py>" \
  --stage-id pre_gate2_decontx_provisional \
  --code-revision <resolved from DEPLOY_ROOT's own git rev-parse HEAD at run time> \
  --expression-source decontx_counts \
  --batch-col cohort --patient-col patient_id --lineage-col lineage --outcome-col outcome \
  --run-harmony --seed 0
```

`--run-harmony` is already part of this committed command, run as an
exploratory/transductive diagnostic step by design — not a "possibly
optional" output, and not a claim that Harmony is being selected as a
production correction method or compared against other methods.

## 4. Input identity and re-checked state (this round)

| Item | Value | When checked |
|---|---|---|
| `input.h5ad` sha256 | `e2da0e4f44350a1ecb28c1894308812f9234c331420672a52177f9445305d6ba` | From `validate_job1_output.py`'s `job1_validation_addendum_v3_20261011_123120/report.json`; not rehashed this round |
| `--out-dir` (`diagnostic_run`) exists? | **No** — re-checked this round (`ls` on the parent dir), immediately before writing this packet | 2026-10-11, this round |
| Active SLURM jobs for this user | **None** (`squeue -u xinruqiu` empty) | 2026-10-11, this round |

A submission, if authorized, must re-run both of these checks again
immediately before `sbatch` — this packet's checks are not a substitute for
a check at submission time.

## 5. The CV0201/CV0326 identity conflict's relevance to this job (conditional, not unconditional)

See `stephenson_identity_decision_sheet_v3.md` §3. `run_batch_diagnostics.py`
never reads pool/`cv_group` identity, so this job's result is unaffected ONLY
under the stated condition: the input cells, `patient_id`, `outcome`,
`lineage`, and parameters all remain exactly as they are now. This is not a
claim that the identity question is resolved.

## 6. Acceptance criteria for a real run, if authorized

A run is acceptable only if ALL of the following are true — not merely "job
exited":

1. SLURM accounting shows `State=COMPLETED` and `ExitCode=0:0` for the job
   (check with `sacct -j <jobid> --format=JobID,State,ExitCode`), not just
   that the wrapper script's own `set -euo pipefail` didn't trip.
2. The wrapper's own `code_revision`/working-tree-cleanliness checks (inside
   the SLURM log) show no `ERROR:` line.
3. `manifest.json`'s top-level `status` field is `"complete"` — NOT
   `"complete_with_partial_results"` or `"failed"`. If it is
   `"complete_with_partial_results"`, the exact incomplete items (which
   `manifest["per_lineage"][i]["embeddings"][label]["status"]` is `"failed"`,
   `"partial"`, or `"unsupported"`, and why) must be listed explicitly before
   this round can be called successful — never summarized as "mostly done."
4. `manifest.json`'s `code_revision_verified` is `True` (not merely a
   non-null `code_revision_asserted`) — per `_resolve_code_identity()`'s own
   documented distinction between an asserted string and a verified one. If
   `code_revision_verified` is `False` because the deployment is a loose-file
   scratch copy rather than a git checkout, `source_file_sha256` must be
   checked against the real committed files instead of treating `False` as
   disqualifying by itself.
5. For each requested lineage × embedding, `run-harmony`'s
   `harmony_corrected` scope is present in `metrics_long.csv` with real
   values (not absent/NaN) wherever the uncorrected scope succeeded —
   otherwise Harmony's requested-but-failed state must be visible in
   `manifest.json`, per `run_batch_diagnostics.py`'s own design (a
   requested-but-failed Harmony correction fails the embedding's status, not
   just an inner field).
6. `metrics_long.csv` and every lineage's saved representation/UMAP files
   referenced in `manifest.json` exist and are readable (non-zero size,
   loadable).
7. The log path above is read in full, not sampled.

If a run comes back `"complete_with_partial_results"` or `"failed"`, this
packet's conclusion is NOT "ready" — it is "report the exact gap and let the
owner decide whether to accept, re-run, or investigate," per `AGENTS.md`.

## 7. What this packet does NOT do

Does not submit `sbatch`. Does not modify
`302006_updategit_deploy_aea9a61` (the wrapper copy went to a separate
`stage33_wrappers/` staging path). Does not change
`canonical_donor_labels.csv`, QC, or DecontX. Does not resolve the
`CV0201`/`CV0326` identity conflict. The decision to run Job 2 remains the
user's.
