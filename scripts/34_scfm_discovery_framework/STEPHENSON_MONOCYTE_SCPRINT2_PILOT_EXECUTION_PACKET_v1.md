# Stephenson monocyte scPRINT-2 technical pilot — execution packet v1

**Date:** 2026-10-16. **Authorization basis:** `docs/decisions.md`, "2026-10-16 — Stephenson
monocytes: provisional raw counts accepted for a technical pilot" (Cindy's direct "允许", in
response to the scoped proposal in `review:20261016_130032`). **Scope of that
authorization:** input → embedding/state-discovery technical validation only, full visit/
donor/sample/cell identity retained, Gate 1 flags and historical DecontX scores attached as
annotations only (never used to filter cells or replace counts). Results stay provisional —
NOT evidence for disease/pathway conclusions. **Not authorized by this packet or that
decision:** SLURM submission, dependency installation, push, model/checkpoint final
selection, or any scientific DE/masking conclusion.

**Status at a glance:** `implemented` = code written and present in this repo.
`tested` = exercised by a passing automated test (synthetic fixture unless stated).
`ready, not_run` = a complete, reviewable command exists but was not executed this round.
`not_run` = no execution of any kind this round.

| Component | Status |
|---|---|
| Input-prep code (GEX subset + QC annotation) | **implemented, tested** (synthetic fixture) |
| Input-prep run against the real 5.7 GB object | **not_run** (needs a SLURM allocation; login-node run forbidden by `DEVELOPMENT.md`) |
| QC source cross-check (cell-identity join) | **implemented, tested, and independently re-run against real data** (128,503/128,503 joined, 0 mismatches) |
| `scprint2 embed` CLI itself | **verified present and responsive** (`--help`/`embed --help` captured live, 2026-10-16) |
| scPRINT-2 small-v2.ckpt on a *different* real human cohort | **tested, real artifact exists** (`option_c/scprint2_covid_canonical.h5ad`, 399,538,242 bytes, 2026-05-01) |
| scPRINT-2 embed run on THIS monocyte input | **ready, not_run** — SLURM script written, not submitted |
| Model/checkpoint final selection | **candidate only** — scPRINT-2 small-v2.ckpt is the one with the strongest direct evidence (see §2), NOT a confirmed choice; UCE 4-layer remains a named alternative (`model_registry.py`) |

## 1. Entry points (existing vs. newly written)

**Existing, reused unchanged:**
- `scripts/foundation_models/scPRINT/run_scprint2_human_covid.sh` — the ONLY proven-working
  real-data invocation of this project's `scprint2_env` + `scprint2 embed` CLI. Its exact
  flag set is the basis for §4's command below.
- `scripts/34_scfm_discovery_framework/code/model_registry.py`'s `SCPRINT2_SMALL_V2` entry —
  checkpoint path/hash/species evidence, already verified in a prior round, re-used unchanged.
- `scripts/34_scfm_discovery_framework/results/stephenson_sample_subtype_cell_counts.py`,
  `stephenson_qc_source_crosscheck.py` — supply the already-verified cell-identity join this
  packet's annotation step depends on.

**Newly written this round (the "one small missing piece" the prior review asked to be
named specifically, not glossed over):**
- `results/stephenson_monocyte_scprint2_pilot_prep.py` — the real object mixes 192
  "Antibody Capture" (CITE-seq) features into its 24,929-feature var axis (confirmed via a
  bounded h5py read, 2026-10-16: `var['feature_types']` categories are `Antibody Capture`
  (192) / `Gene Expression` (24,737)). The ONE proven-working scPRINT-2 input in this project
  (`human_platelet_covid_severity_raw_counts.h5ad`) carries ONLY `_index` in `var` — i.e. it
  was ALREADY Gene-Expression-only before scPRINT-2 ever saw it. No existing code in this
  repo performs the equivalent subset for the Stephenson monocyte object, so this script does
  exactly that, plus attaches the already-verified Gate 1/DecontX annotations. **Implemented
  and tested against a synthetic fixture; NOT run against the real object this round**
  (reading/writing the real X requires materializing a 5.7 GB object — forbidden on the login
  node by `DEVELOPMENT.md`, and no SLURM job is authorized this round).
- `scripts/foundation_models/scPRINT/run_scprint2_stephenson_monocyte_pilot.sh` — the actual
  `scprint2 embed` call for this pilot, adapted from the proven script. **Written, not
  submitted.**

## 2. Candidate model/checkpoint (still a candidate — not the owner's confirmed choice)

**scPRINT-2, checkpoint `small-v2.ckpt`** (`model_registry.SCPRINT2_SMALL_V2`):
- Checkpoint file present and hash-verified against the official HuggingFace page
  (`2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822`).
- Species: human + mouse both `checkpoint_metadata_confirmed` (read from the checkpoint's own
  pickled `organisms` list, not just a README) — Stephenson is human, so this is the
  stronger-evidence case for this specific pilot.
- **Real prior successful run on a DIFFERENT real human cohort in THIS project**:
  `02_EMBEDDINGS/option_c/scprint2_covid_canonical.h5ad` (399,538,242 bytes, dated 2026-05-01),
  produced by the exact same env (`scprint2_env`) and CLI call this packet reuses. This is
  concrete runnability evidence for the tool itself — it is NOT evidence that this specific
  monocyte input, its cell-ID alignment, or its current QC version have been validated against
  it (per `MODEL_SELECTION_v1.md`'s own caution, unchanged).
- CLI responsiveness independently re-confirmed 2026-10-16: `scprint2 --help` and `scprint2
  embed --help` both returned real, current flag documentation from the installed env (not
  from memory of an older version) — see §4 for the flags this packet actually uses.

**Named alternative, not excluded**: UCE 4-layer checkpoint (`model_registry.UCE_4LAYER`) —
present on disk, species `documented_unconfirmed` (my own attempt to re-read the primary
source was blocked, per that entry's own note), no real human run logged in this specific
framework yet. Not chosen here because scPRINT-2's species/runnability evidence is directly
stronger for this specific pilot, not because UCE is ruled out.

**This packet does not confirm a final backbone.** That confirmation is the owner's, as
`MODEL_SELECTION_v1.md` already states and this round does not change.

## 3. Input layer, gene axis, and preprocessing contract

| Property | Value | How confirmed |
|---|---|---|
| Source object | `02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad` | existing, already hash-recorded |
| n_obs (cells) | 128,503 | bounded h5py read |
| n_vars (raw) | 24,929 (24,737 Gene Expression + 192 Antibody Capture) | bounded h5py read, 2026-10-16 |
| Gene ID namespace | HGNC gene SYMBOLS (e.g. `MIR1302-2HG`, `SAMD11`, `ISG15`) | bounded h5py read of `var/_index`, 2026-10-16 |
| X layer | raw integer counts stored as `float32` (CSR sparse) | bounded h5py read of `X/data` + `X.attrs['encoding-type']` |
| scPRINT-2's accepted gene-ID namespaces | `ensembl_gene_id` OR `gene_symbol` (official README FAQ, re-verified 2026-10-12) | matches this object's symbol-based var_names directly — no gene-ID conversion needed |
| scPRINT-2's accepted count type | raw counts, accepted directly per README | matches this object's X directly |
| Required preprocessing before `scprint2 embed` | subset `var` to `feature_types == "Gene Expression"` (192 columns dropped) | this round's new finding — see §1 |
| `--is_symbol` flag | left UNSET (CLI default `False`) | the ONE proven-working call used symbol-named input (`human_platelet_covid_severity_raw_counts.h5ad`: `LINC00115`, `FAM41C`, `NOC2L`, ...) with this flag unset and completed successfully — replicating the evidenced behavior rather than guessing a "more correct" value with zero supporting evidence |
| `--use_raw`, `--skip_validate`, `--how`, `--max_len`, etc. | left at CLI defaults, unchanged from the proven call | same reasoning — no evidence motivates deviating |
| Output gene axis | scPRINT-2's OWN internal vocabulary (observed in the proven run: input 13,066 genes → output 85,720 genes, ALL Ensembl IDs) | the model does not preserve the input's gene list; `embed` output var differs in identity and count from input var — must not be conflated |
| Cell retention through `embed` | NOT 1:1 in the one proven run (47,977 in → 47,114 out, 863 fewer) | observed directly in the existing real output file; cause not investigated this round — named as an acceptance-check item, §5 |

## 4. Commands (prepared, not executed)

**Step 1 — prepare input** (requires a SLURM allocation with enough memory to materialize the
5.7 GB object; NOT run on the login node):
```
python stephenson_monocyte_scprint2_pilot_prep.py
```
Reads `RAW_H5AD`/`GATE1_FLAGS_CSV`/`DECONTX_SCORES_CSV` (module-level constants, all existing,
already-verified paths), writes `OUT_H5AD` =
`02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_scprint2_pilot_input_v1.h5ad` (refuses if
it already exists) plus a `_summary.json` report.

**Step 2 — embed** (SLURM, GPU; script written, NOT submitted):
```
sbatch scripts/foundation_models/scPRINT/run_scprint2_stephenson_monocyte_pilot.sh
```
Which runs, inside the allocation:
```
scprint2 embed \
    --adata /bigdata/.../stephenson_monocyte_scprint2_pilot_input_v1.h5ad \
    --ckpt_path /bigdata/.../scprint2_ckpts/small-v2.ckpt \
    --species NCBITaxon:9606 \
    --output_filename /bigdata/.../stephenson_monocyte_scprint2_pilot_embeddings_v1.h5ad
```

## 5. Resource ceiling and basis

Reused, unchanged, from `run_scprint2_human_covid.sh`: `--partition=short_gpu`, 1 GPU, 8 CPUs,
128G RAM, 2h walltime, `--requeue`. **Basis**: this is the ONLY resource spec in this project
that has ever produced a complete real scPRINT-2 `embed` output (the proven run above).
**Caveat, stated plainly, not glossed over**: that proven run's input was 47,977 cells x
13,066 genes; this pilot's prepared input will be 128,503 cells x 24,737 genes — larger on
both axes. Whether the same 2h/128G ceiling is SUFFICIENT for this larger input is
**UNVERIFIED** — this packet reuses it as the only evidenced starting point, not as a
confirmed-adequate value. Fail-stop only: the script exits non-zero on any missing
input/checkpoint/pre-existing-output condition; no automatic retry beyond SLURM's own
`--requeue` (preemption handling, not a logic retry).

## 6. Identity and numeric acceptance criteria (technical validation only — NOT a scientific-quality judgment)

1. **Input traceability**: every source file's sha256 recorded in the prep step's
   `_summary.json` (already implemented, `build_prepared_adata`/`main` in
   `stephenson_monocyte_scprint2_pilot_prep.py`).
2. **Identity/row-order preservation**: prepared input's `obs_names` == raw object's
   `obs_names`, same order, same count (128,503) — enforced by `validate_prepared_adata()`,
   mutation-killed (a cell-count change or a missing QC annotation is each independently
   caught; see `tests/test_stephenson_monocyte_scprint2_pilot_prep.py`).
3. **Output dimensions and numeric finiteness** (after `embed`, not yet run): `obsm` must gain
   at least one `X_scprint*`/`scprint_emb*` key (matching the proven run's own obsm keys);
   every value in that embedding must be finite (no NaN/Inf); `n_obs` of the output must be
   reported and compared against 128,503 EXPLICITLY (per §3's finding that the proven run was
   NOT 1:1) — a difference is not itself a failure, but it must never go unreported.
4. **Model configuration/seed/dependency record**: `--seed_everything` CLI default (42,
   confirmed via `--help`, 2026-10-16) recorded in the job log; checkpoint sha256 and env path
   recorded (already in `model_registry.py`).
5. **Re-readable result and run status**: output must open with `anndata.read_h5ad` without
   error and expose the obsm keys above; the SLURM `.out`/`.err` logs are the run-status
   record (per `DEVELOPMENT.md`: read the log the job actually writes, don't infer completion from
   absence of a `/tmp` file or an empty `pgrep`).

**Explicitly NOT an acceptance criterion this round**: whether the resulting embedding shows
"good-looking" clusters, separates disease/severity groups, or supports any biological
interpretation — all deferred to the real scientific-pilot phase, which this technical
validation does not substitute for.

## 7. What remains for the owner

- Final model/checkpoint confirmation (§2) — a scientific choice, to be brought to Cindy with
  this evidence once she is ready to make it; not resolved by this packet.
- Per-job SLURM execution authorization for `run_scprint2_stephenson_monocyte_pilot.sh` —
  separate from, and not granted by, the provisional-input-use decision this packet executes
  against.
- Whether the §5 resource ceiling needs adjustment for this input's larger size — can only be
  answered by either checking scPRINT-2's own memory-scaling documentation (not done this
  round) or an actual run.
