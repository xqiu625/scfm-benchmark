# Stephenson monocyte scPRINT-2 technical pilot — execution packet v3

**v1 and v2 are left on disk unchanged as the historical record of what those rounds produced
— v2's output validator and job-script deployment binding had real gaps, closed below, per
`review:20261016_152556`.** This version is the current one.

**Date:** 2026-10-16. **Authorization basis (unchanged):** `docs/decisions.md` — Cindy's "允许"
(provisional raw-count use for an input→embedding/state-discovery TECHNICAL pilot only; full
visit/donor/sample/cell identity retained; Gate 1 flags and historical DecontX scores attached
as annotations only, never filtering cells or replacing counts; results stay provisional) and
her "ok" (scPRINT-2 `small-v2.ckpt` as the FIRST technical-pilot model, not a confirmed final
manuscript backbone). **Not authorized by this packet or those decisions:** SLURM submission,
dependency installation, push, final model/checkpoint selection, or any scientific DE/masking
conclusion.

## What changed from v2 (do not re-read v2 as current)

v2 fixed R1-R6 from an earlier review but had two remaining real gaps, found by the same
review's follow-up pass and closed this round:

| # | v2 gap | Fix |
|---|---|---|
| Blocker 1 | `stephenson_monocyte_scprint2_pilot_validate.py` took only `expected_n_obs` — a same-count output made of unrelated ids, duplicate ids, or misrouted QC annotations would have passed | now takes the frozen prepared-input `.obs` as `expected_obs` and checks exact cell-id coverage (no missing/extra/duplicate ids), the real primary `scprint_emb` key specifically (not any `X_scprint*`-prefixed key — the installed `scprint2/tasks/cell_emb.py:214,227-251` writes a UMAP to a SEPARATE `X_scprint_umap` key, which the old prefix match wrongly accepted alone), array dtype/shape (rejects string or zero-column arrays), and QC-column correspondence BY ID (recoverable if dropped, a hard failure if present but diverging) |
| Blocker 2 | the job script's `PREP_RESULTS_DIR` pointed at the LIVE `302006_updategit` checkout, confirmed again this round to still be stale (`35131cf`) with unrelated uncommitted changes | job script now points at `302006_updategit_deploy_scprint2_pilot`, an INDEPENDENT clean export of commit `d0aa764` (verified below, §1); the stale checkout itself was never read, modified, or cleaned |

Also fixed, the review's two smaller findings: `validate_prepared_adata()`'s checks used bare
`assert`, which `python -O`/`PYTHONOPTIMIZE` silently strips — rewritten via an explicit
`_require()` helper (verified: reverting to `assert` and re-running the new `-O`-subprocess test
makes it fail, confirming the fix is real, not cosmetic). And the no-op-validator mutation test
previously only checked the no-op's own return value, not that it actually broke the dependent
tests — rewritten to call those tests directly under the monkeypatched no-op and confirm they
fail.

Carried over unchanged from v2 (already fixed, not revisited): R1 (output-suffix handling), R2
(`--is_symbol true`), R3 (cell-filter finding for this object), the checkpoint SHA256 runtime
check, `--no-requeue`, and the fail-stop prep→embed→validate bundling.

## Status at a glance

| Component | Status |
|---|---|
| Input-prep code (GEX subset + QC annotation) | **implemented, tested** (synthetic fixture + an actual h5ad write/read-back round trip); acceptance checks now raise via `_require()`, confirmed to survive `python -O` |
| Output identity/value validator (`stephenson_monocyte_scprint2_pilot_validate.py`) | **rewritten this round** — checks cell-id coverage/duplicates, the real `scprint_emb` key, array dtype/shape, and QC correspondence by id, not just a count; all 9 of the review's own probe scenarios now covered by tests and 2 independent mutation-kills (primary-key gating, id-coverage checks) |
| QC source cross-check (`stephenson_qc_source_crosscheck.py`) | unchanged from v2 — implemented, tested, mutation-killed, independently re-run against real data (128,503/128,503 joined, 0 mismatches) |
| Gene-symbol mapping boundary (`is_symbol=True`) | tested against the real installed gene table (unchanged from v2) — 19,194/24,737 genes map |
| Built-in `min_nnz_genes=200` cell filter, for THIS object | tested against the real object's real sparsity structure (unchanged from v2) — 0/128,503 cells below threshold |
| Deployment binding | **new this round** — `302006_updategit_deploy_scprint2_pilot` is an independent clean git checkout at commit `d0aa764`, verified clean and hash-identical to the authoring commit (§1); the job script refuses to run if that checkout later becomes dirty |
| Working-directory isolation | **new this round** — job script now `cd`s into the per-run directory before calling `scprint2 embed`, so the installed library's relative `data/` intermediate-file writes land inside that run's own directory |
| Input/QC/checkpoint fingerprints | **new this round** — verified BEFORE the run (fail closed), not only recorded after |
| Run parameters (`seed_everything`/`batch_size`/`num_workers`/`how`/`max_len`) | **pinned explicitly this round**, matching the installed `Embedder`/CLI's own current defaults (confirmed via `jsonargparse.ArgumentParser.add_class_arguments`'s `nested_key=None` default, and via reading `Embedder.__init__`'s own defaults directly) |
| Prep → embed → validate on THIS monocyte input | **ready, not_run** — bundled job script written, not submitted |
| Model/checkpoint | small-v2.ckpt is the confirmed FIRST TECHNICAL-PILOT model (owner "ok"); final manuscript backbone still unconfirmed, UCE not excluded |

## 1. Deployment: an independent, verified snapshot — the stale checkout untouched

`/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot` was created this
round as a clean `git clone` of the local repository at the exact commit that fixed both
blockers, then rsynced to HPCC (working tree + full `.git`, never the live `302006_updategit`
working directory). Verified directly on HPCC, 2026-10-16:

```
$ git -C 302006_updategit_deploy_scprint2_pilot rev-parse HEAD
d0aa76416756f53ed71f9ba6c59268ce8a621765
$ git -C 302006_updategit_deploy_scprint2_pilot status --porcelain | wc -l
0
$ sha256sum .../stephenson_monocyte_scprint2_pilot_validate.py   # matches local exactly
07d88e4b5593a1f4ec0b951c5b6af09f3ed2bd0f3dc355c94a14b54dda3c3349
$ sha256sum .../stephenson_monocyte_scprint2_pilot_prep.py       # matches local exactly
43ae5babcf1e2e61f80dd4c3ebba208a954425fdd3c2548281c2b2d05d35d8d4
$ sha256sum .../run_scprint2_stephenson_monocyte_pilot.sh        # matches local exactly
fcde64eaa73407218343fe1e9dabe173bc59ce6e6ab4e3fd61a29696af846d52
```

**The stale checkout (`/bigdata/.../302006_updategit`, still at `35131cf` with unrelated
uncommitted changes) was never read, modified, pulled, reset, or cleaned this round** — it
remains exactly as the prior round found it, untouched, per the review's explicit instruction.
The job script's `PREP_RESULTS_DIR` now resolves inside the new deploy directory and the job
refuses to start (`exit 1`) if that directory is not a git checkout or is not clean at runtime —
so a later direct edit to the deployment itself, not only to the old checkout, is also caught.

Re-deploying a future fix means repeating this same export-and-verify procedure to a NEW commit
inside the SAME deploy directory (or a freshly suffixed one); it is not a code change to the job
script, which reads `DEPLOY_COMMIT` from the live `git rev-parse HEAD` of whatever is there.

## 2. Real verification of the two preprocessing boundaries (unchanged from v2)

Unchanged this round — both probes (gene-symbol mapping, built-in cell filter) were already
verified against the real installed library and real object in the prior round; not re-run,
since nothing here changed. See v2 §2 for the full evidence.

## 3. Output identity/value validator — what it now actually checks

`stephenson_monocyte_scprint2_pilot_validate.py`'s `validate_embedding_output(embedded_adata,
expected_obs, ...)` takes the frozen prepared-input object's real `.obs` (the SAME object given
to `scprint2 embed`), not an integer count, and recomputes every check by cell id:

| Check | Hard failure if... | Verified by |
|---|---|---|
| Duplicate ids | any id repeats in the output | test + mutation-kill |
| Missing ids | any expected id is absent (this round's scope does not accept a partial embedding as a technical PASS) | test + mutation-kill |
| Extra/unrelated ids | any output id is not in the expected manifest | test + mutation-kill |
| Primary embedding key | `scprint_emb` specifically is absent (an `X_scprint_umap`-only file does NOT pass) | test + mutation-kill; confirmed against the installed `cell_emb.py:214,227-251` directly |
| Embedding dtype/shape | non-numeric (e.g. string) dtype, or not 2-D with ≥1 column | test |
| Embedding finiteness | any non-finite value in the primary embedding | test |
| QC-column correspondence | present but diverging from the frozen manifest for a shared id (misrouted annotation) | test |
| QC-column recovery | missing from the output but recoverable by id ONLY when ids fully match; otherwise a hard failure, never a silent pass | test |

Row-order changes (same ids, different order) are recorded (`order_matches_expected`) but never
assumed correct for any other check — every comparison above is by id, via pandas `.loc`
indexing, never by row position.

Two mutation-kills performed and verified this round, on the real file (not a copy with
unverified provenance): (a) reverting the primary-key gate to the old prefix-based check made
exactly the 4 tests tied to it fail; (b) removing the missing/extra/duplicate-id checks made
exactly the 5 tests tied to them fail. Both reverted afterward and hash-verified against the
committed file.

## 4. Commands (prepared, not executed)

**One bundled job**, now pointing at the verified deployment:
```
sbatch /bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot/scripts/foundation_models/scPRINT/run_scprint2_stephenson_monocyte_pilot.sh
```
Internally, after the deployment-cleanliness check, the three input/QC sha256 checks, and the
checkpoint sha256 check (all fail closed, before any GPU time is spent):
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
failure) — not merely whether `scprint2 embed` itself exited 0. A `provenance.json` (installed
package versions, key installed-source-file hashes, deployment commit, pinned run parameters)
is written to the run directory before the prep step runs.

**No prerequisite gap remains for this job's own code path** — §1 above closes the one the prior
round left open. The only remaining item is the SLURM authorization itself (§7).

## 5. Resource ceiling and basis (unchanged from v2 — still genuinely unverified)

`--partition=short_gpu`, 1 GPU, 8 CPUs, 128G RAM, 3h walltime, `--no-requeue`. **Basis**: reused
from the one proven real scPRINT-2 `embed` job in this project (47,977×13,066 input); this
pilot's input is 128,503×24,737 — larger on both axes. **Not claimed sufficient.** No new
measurement was taken this round (would require either an actual run or reading scPRINT-2's own
memory-scaling documentation, neither done); this remains the best available starting point, not
a verified budget. A timeout or OOM fails closed under `set -euo pipefail` and `--no-requeue` —
it does not corrupt the run directory or silently produce a partial result that would pass
validation, since the validator requires the output file to exist and pass identity checks.

## 6. Identity and numeric acceptance criteria — wired into the run path

1. **Input traceability**: every source file's sha256 verified BEFORE the run starts (new this
   round — previously only recorded after), plus recorded again in the prep step's own
   `_summary.json` and in `provenance.json`.
2. **Identity/row-order preservation** (prep step): re-derived directly from the real objects by
   `validate_prepared_adata()`, raised via `_require()` (survives `python -O`), mutation-killed
   against the real production file.
3. **Original counts preserved exactly**: sparse diff against the original GEX-subset `X`,
   never densified, never merely "finite".
4. **Output identity, dimensions, and numeric finiteness** (after `embed`): now cell-id-aware,
   not merely count-aware — see §3's full table. Hard failures produce a real nonzero exit code
   from the bundled job itself.
5. **Model configuration/seed/dependency record**: `seed_everything`/`batch_size`/`num_workers`/
   `how`/`max_len` now PINNED explicitly (not left to the installed library's current defaults,
   which could change silently on a future library upgrade); checkpoint sha256 verified at
   runtime; installed package versions and key source-file hashes recorded in `provenance.json`.
6. **Deployment identity**: the exact deployed git commit is verified clean and recorded at
   runtime (§1), not merely assumed from a path that happens to exist.
7. **Re-readable result and run status**: output opens with `anndata.read_h5ad`; the SLURM
   `.out`/`.err` logs plus `validation_report.json` and `provenance.json` are the run-status
   record.

**Explicitly NOT an acceptance criterion**: whether the embedding shows good-looking clusters or
supports any biological interpretation.

## 7. What remains for the owner

- **Per-job SLURM execution authorization for the bundled job** — this is the one remaining item
  blocking this specific run. It is separate from, and not granted by, the provisional-input-use
  or first-pilot-model decisions already given.
- Whether the 3h/128G ceiling is actually sufficient — genuinely unverified (§5); a failure here
  fails closed, it does not corrupt anything or produce a false pass.
- Final manuscript-backbone model/checkpoint confirmation remains separate and not addressed by
  this packet, per Cindy's own framing of the "ok" decision.
