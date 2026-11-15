# Stephenson monocyte scPRINT-2 technical pilot — execution packet v4

**v1, v2, and v3 are left on disk unchanged as the historical record of what those rounds
produced — v3's output validator and job-script time/revision binding had real gaps, closed
below, per `review:20261016_202903`.** This version is the current one.

**Date:** 2026-10-16. **Authorization basis (unchanged):** `docs/decisions.md` — Cindy's "允许"
(provisional raw-count use for an input→embedding/state-discovery TECHNICAL pilot only; full
visit/donor/sample/cell identity retained; Gate 1 flags and historical DecontX scores attached
as annotations only, never filtering cells or replacing counts; results stay provisional) and
her "ok" (scPRINT-2 `small-v2.ckpt` as the FIRST technical-pilot model, not a confirmed final
manuscript backbone). **Not authorized by this packet or those decisions:** SLURM submission,
dependency installation, push, final model/checkpoint selection, or any scientific DE/masking
conclusion.

## What changed from v3 (do not re-read v3 as current)

v3 closed the deployment and output-identity gaps from an earlier review, but the review's
follow-up pass (`review:20261016_202903`) found 4 real, independently reproduced
problems in v3 itself, all closed this round:

| # | v3 problem | Fix |
|---|---|---|
| R1 | `--time=3:00:00` exceeds `short_gpu`'s real `MaxTime=02:00:00` (`scontrol show partition short_gpu`, queried directly) — the "ready" `sbatch` command could be rejected or killed by the scheduler, independent of whether 3h was scientifically enough | lowered to `--time=2:00:00`. Whether 2h is actually ENOUGH for this larger input remains genuinely **unverified** — not newly measured this round, not assumed sufficient; a timeout still fails closed (`--no-requeue`, no auto-retry) |
| R2 | the validator's `_qc_column_report` set `recovered_from_expected_manifest=True`/`n_recovered` in the JSON report for a missing-but-id-matching column, but nothing anywhere actually wrote that column back — the review's own `persisted_missing_qc` probe called `main()` on a real 3-cell h5ad, reopened the written output, and found both QC columns still absent while the validator had returned exit code 0 (PASS) | removed the fake-recovery branch entirely: a required column missing from the real output is now **always** a hard failure. No sidecar/recovery mechanism was added (scope stayed to the minimal fix the review itself offered) |
| R3 | the validator checked only 2 QC columns (`gate1_flag_would_delete`, `decontx_estimated_contamination`) — no donor/sample identity check at all. The review's `swapped_patient_id`/`swapped_Sample_ID` probes (same ids, same correct primary embedding, same correct QC values, only the identity column's values reversed) both returned `hard_failures=[]` | `patient_id` and `Sample_ID` added to the required-columns set, checked by id with the same (now simplified) logic as the QC columns. `visit` was the review's own synthetic probe field name, not a claim the real object has it — confirmed directly: a bounded h5py read of `stephenson_monocyte_raw.h5ad`'s real obs columns (2026-10-16) returned `['Data_NO', 'Sample_ID', '_index', 'disease', 'monocyte_source_provenance', 'orig_celltype', 'patient_id', 'source_study', 'y']` — there is no `visit` column on this object, so none was invented |
| R4 | the job script's clean-checkout check recorded whatever `HEAD` it found at run time — an updated-but-still-clean deployment (e.g. fast-forwarded while this job sat queued) would still pass | the submitter must now pass the exact pinned commit via `sbatch --export=ALL,EXPECTED_DEPLOY_COMMIT=<hash>` (not baked into the script file itself, which would be self-referential); the job refuses to run if that commit is unset or doesn't match the deployment's actual `HEAD` |

Also corrected: the v3 handoff claimed "49 tests across the 2 touched test files pass" — the
review's own independent run (and this round's) both show the real number is **38** for those
2 files (19 + 19); this packet and the handoff below use the real, re-verified counts.

## Status at a glance

| Component | Status |
|---|---|
| Output identity/value validator (`stephenson_monocyte_scprint2_pilot_validate.py`) | **fixed this round (R2, R3)** — no more fake-recovery claim (fail-closed on any missing required column); now also checks `patient_id`/`Sample_ID` by id, not just the 2 QC columns |
| Job script resource request | **fixed this round (R1)** — `--time` lowered to the partition's real `MaxTime=2:00:00`; sufficiency still unverified, documented honestly |
| Job script revision binding | **fixed this round (R4)** — requires an externally supplied `EXPECTED_DEPLOY_COMMIT`, compared against the deployment's live `HEAD` at run time, not merely "is it clean" |
| Deployment | **re-exported this round** to commit `3de4bca` (was `d0aa764`) — same independent directory, same stale-checkout-untouched guarantee, verified clean and hash-identical to local (§1) |
| Input-prep code, QC source cross-check, gene-symbol mapping, built-in cell filter | unchanged from v3 — not revisited this round, nothing in this review touched them |
| Prep → embed → validate on THIS monocyte input | **ready, not_run** — bundled job script written, not submitted |
| Model/checkpoint | small-v2.ckpt is the confirmed FIRST TECHNICAL-PILOT model (owner "ok"); final manuscript backbone still unconfirmed, UCE not excluded |

## 1. Deployment: re-exported to the fixed commit, stale checkout still untouched

`/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot` (the SAME
directory created in the v3 round) was re-exported this round: a fresh clean `git clone
--no-local` of the local repository at the new fixed commit, rsynced to the existing HPCC
directory (working tree + `.git`, no files deleted — git objects are additive). Verified
directly on HPCC, 2026-10-16:

```
$ git -C 302006_updategit_deploy_scprint2_pilot rev-parse HEAD
3de4bca47ad4ad07d501f557d7b9065540477655
$ git -C 302006_updategit_deploy_scprint2_pilot status --porcelain | wc -l
0
$ sha256sum .../stephenson_monocyte_scprint2_pilot_validate.py   # matches local exactly
c134e6e5ae25d16bd0392413c042772525961b0298ff75c1a006e6a62bc65f1f
$ sha256sum .../run_scprint2_stephenson_monocyte_pilot.sh        # matches local exactly
4bca43007d907754390629e0ac641ebcb237c1d59b763f7fb2306634ba325fd5
$ sha256sum tests/test_stephenson_monocyte_scprint2_pilot_validate.py  # matches local exactly
66f2627f1625e18c0cfdb1717b791ed4cde5ceaa229c883edbeaac6afce69c1e
```

**The stale checkout (`/bigdata/.../302006_updategit`, still at `35131cf` with unrelated
uncommitted changes) was never read, modified, pulled, reset, or cleaned this round either** —
re-confirmed directly: `git -C /bigdata/.../302006_updategit rev-parse HEAD` still returns
`35131cf38da5485c4010676669173ca98f57630f`.

The job script's own R4 fix (§ above) now also requires the **submitter** to assert the exact
expected commit at `sbatch` time, closing the gap the review found in the prior round's
record-whatever-HEAD-is-there check.

## 2. Output identity/value validator — what changed, verified against the review's own probes

Re-ran the review's exact probe methodology (`persisted_missing_qc`, `swapped_patient_id`,
`swapped_Sample_ID`) directly against the fixed code on HPCC:

| Probe | v3 result (the review's own run) | v4 result (this round, re-run against the fix) |
|---|---|---|
| Missing required column, ids fully matching | `exit_code=0` (PASS); reopened output still missing both columns | `exit_code=1` (FAIL); 4 "is missing from the output -- fail closed" hard failures, one per required column |
| `patient_id` values reversed (ids/embedding/QC all otherwise correct) | `hard_failures=[]` | 1 hard failure: `"obs['patient_id'] diverges from the frozen prepared-input manifest for 2 cell id(s)"` |
| `Sample_ID` values reversed | `hard_failures=[]` | 1 hard failure: `"obs['Sample_ID'] diverges from the frozen prepared-input manifest for 2 cell id(s)"` |

Two independent mutation-kills performed on the real file this round (HPCC scratch copy,
hash-matched to the committed file before and after):

1. Reverted the missing-column branch back to the old fake-recovery logic → exactly the 3 tests
   tied to it failed (`fails_hard_on_a_missing_required_column_even_when_ids_fully_match`,
   `..._when_ids_do_not_fully_match`, `fails_hard_on_a_missing_identity_column`).
2. Removed `patient_id`/`Sample_ID` from the required-columns default → exactly the 5 tests tied
   to them failed (the well-formed-pass test plus both identity-column tests × 2 scenarios each).

Both reverted afterward; `sha256sum` matched the committed file exactly both times.

Carried over unchanged from v3 (not revisited, nothing in this review touched them): exact
cell-id coverage, duplicate-id rejection, the real `scprint_emb` primary-key gate, array
dtype/shape/finiteness checks, and the by-id (never by-position) comparison discipline.

## 3. Resource request — now within the real partition limit, sufficiency still open

```
$ scontrol show partition short_gpu
PartitionName=short_gpu
   DefaultTime=02:00:00 ... MaxTime=02:00:00 ...
   State=UP ...
```

`--time=3:00:00` (v2/v3) exceeded this. Lowered to `--time=2:00:00`. **This is a legality fix,
not a sufficiency claim** — whether 2h/128G is actually enough for this pilot's 128,503-cell,
24,737-gene input is still genuinely unverified; no new measurement was taken this round (would
require either an actual run or scPRINT-2's own memory/time-scaling documentation, neither done).
A timeout under `set -euo pipefail` + `--no-requeue` fails closed: it does not corrupt the run
directory, and the validator still requires the output file to exist and pass every identity
check, so a truncated run cannot produce a false PASS.

## 4. Commands (prepared, not executed)

**One bundled job**, pinned to the exact redeployed commit:
```
sbatch --export=ALL,EXPECTED_DEPLOY_COMMIT=3de4bca47ad4ad07d501f557d7b9065540477655 \
  /bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot/scripts/foundation_models/scPRINT/run_scprint2_stephenson_monocyte_pilot.sh
```
Internally, after the `EXPECTED_DEPLOY_COMMIT` check, the deployment-cleanliness+revision check,
the three input/QC sha256 checks, and the checkpoint sha256 check (all fail closed, before any
GPU time is spent):
```
python -c "... stephenson_monocyte_scprint2_pilot_prep.main(out_h5ad=<run_dir>/...)"
cd <run_dir>   # working-directory isolation for the library's relative intermediate writes
scprint2 embed --adata <prepared.h5ad> --ckpt_path <small-v2.ckpt> \
    --species NCBITaxon:9606 --is_symbol true \
    --seed_everything 42 --batch_size 64 --num_workers 8 --how "random expr" --max_len 2000 \
    --output_filename <run_dir>/stephenson_monocyte_scprint2_pilot_embeddings
python -c "... stephenson_monocyte_scprint2_pilot_validate.main(<run_dir>/..._embedded.h5ad, <run_dir>/stephenson_monocyte_scprint2_pilot_input.h5ad, <run_dir>/validation_report.json)"
```
Real exit code of the whole job is the validator's exit code (0 = technical pass, 1 = hard
failure) — not merely whether `scprint2 embed` itself exited 0.

**No prerequisite gap remains for this job's own code path** — §1-§3 above close the ones the
prior round left open. The only remaining item is the SLURM authorization itself (§6).

## 5. Tests and verification, this round (real commands/logs, not reused claims)

- `python -m pytest tests/test_stephenson_monocyte_scprint2_pilot_validate.py
  tests/test_stephenson_monocyte_scprint2_pilot_prep.py -q` on HPCC `302006_env_clean`:
  **41 passed, 1 warning in 95.08s** (41 = 22 in the validate file [19 carried over + 3 net new],
  19 unchanged in the prep file; the prior round's "49" claim was wrong, re-measured as 38 for
  those same 2 files before this round's 3 net-new tests, confirmed by both the review's
  independent run and this one).
- The `EXPECTED_DEPLOY_COMMIT` guard logic was tested twice: first against a disposable stand-in
  git directory (3 cases — correct revision, fast-forwarded/clean-but-wrong revision, correct
  new revision), then against the REAL redeployed `302006_updategit_deploy_scprint2_pilot`
  directory itself (correct commit passes, the stale v3 commit `d0aa764` is correctly rejected
  as a mismatch, an unset variable is correctly rejected by the `:` parameter check) — not only a
  stand-in, closing the production-boundary loop per DEVELOPMENT.md §7.
- `bash -n` syntax check of the full v4 job script passed on the HPCC login node.
- No real cohort prep, model inference, or SLURM job was run this round.

## 6. What remains for the owner

- **Per-job SLURM execution authorization for the bundled job**, pinned to commit `3de4bca` —
  this is the one remaining item blocking this specific run. It is separate from, and not
  granted by, the provisional-input-use or first-pilot-model decisions already given.
- Whether the 2h/128G ceiling (now at least a legal request, per §3) is actually sufficient —
  genuinely unverified; a failure here fails closed, it does not corrupt anything or produce a
  false pass.
- Final manuscript-backbone model/checkpoint confirmation remains separate and not addressed by
  this packet, per Cindy's own framing of the "ok" decision.
