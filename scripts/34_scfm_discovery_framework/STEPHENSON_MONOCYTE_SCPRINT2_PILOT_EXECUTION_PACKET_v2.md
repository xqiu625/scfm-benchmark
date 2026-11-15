# Stephenson monocyte scPRINT-2 technical pilot — execution packet v2

**v1 (`STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v1.md`) is left on disk unchanged
as the historical record of what that round produced — several of its claims were wrong,
corrected below, per `review:20261016_141242`.** This version is the current one.

**Date:** 2026-10-16. **Authorization basis (unchanged):** `docs/decisions.md` — Cindy's
"允许" (provisional raw-count use for an input→embedding/state-discovery TECHNICAL pilot
only; full visit/donor/sample/cell identity retained; Gate 1 flags and historical DecontX
scores attached as annotations only, never filtering cells or replacing counts; results stay
provisional) and her "ok" (scPRINT-2 `small-v2.ckpt` as the FIRST technical-pilot model, not a
confirmed final manuscript backbone). **Not authorized by this packet or those decisions:**
SLURM submission, dependency installation, push, final model/checkpoint selection, or any
scientific DE/masking conclusion.

## What changed from v1 (do not re-read v1 as current)

v1 made three claims this round's independent review found wrong, each now fixed with real
evidence (not a guess) from reading the ACTUALLY INSTALLED library source directly:

| # | v1 claim | What the installed source actually does | Fix |
|---|---|---|---|
| R1 | `--output_filename` IS the final output path | installed `scprint2/cli.py` (lines ~261-269) writes to **`<value>_embedded.h5ad`** literally | pass a path STEM; check `<stem>_embedded.h5ad` |
| R2 | `--is_symbol` should stay unset (False), replicating an old run | installed `scdataloader/preprocess.py` (lines ~260-280) merges `adata.var.index` against `genesdf`'s **Ensembl-indexed** rows when `False` — our var_names are SYMBOLS, so `False` would map almost nothing | pass `--is_symbol true`; real probe (below) confirms this maps enough genes |
| R3 | prep not filtering cells means `embed` won't either | installed `scdataloader/preprocess.py` (line ~238) ALWAYS runs `sc.pp.filter_cells(min_genes=200)` regardless of `skip_validate` | real, bounded probe (below) confirms 0/128,503 cells of THIS object would be dropped by it anyway |

Also fixed: R4 (two of this packet's own "mutation-kill" tests were self-defeating — `try:
call(); assert False; except AssertionError: pass` catches its own `assert False`, so they
passed regardless of whether the real call raised; `validate_prepared_adata()` trusted a
caller-supplied report instead of re-deriving identity/order from the actual objects) — both
rewritten, re-verified by mutation-killing the real production file. R5 (the job script lacked
a prep step of its own, used `--requeue`, had no runtime checkpoint hash check, and ended with
only `ls`) — now one bundled, fail-stop job with `--no-requeue`, a runtime SHA256 check, and a
real validation step. R6 (the new `stephenson_qc_source_crosscheck.py` reintroduced the exact
`-1`-code defect fixed elsewhere) — fixed, re-verified against real data.

## Status at a glance

| Component | Status |
|---|---|
| Input-prep code (GEX subset + QC annotation) | **implemented, tested** (synthetic fixture + an actual h5ad write/read-back round trip) |
| `validate_prepared_adata()` identity/numeric guard | **implemented, tested, mutation-killed against the real production file** (6 independent bypass scenarios from the review's own probe table now all caught) |
| Embedding-output validator (`stephenson_monocyte_scprint2_pilot_validate.py`) | **implemented, tested** (synthetic fixture) — new this round |
| QC source cross-check (`stephenson_qc_source_crosscheck.py`) | **implemented, tested, mutation-killed, independently re-run against real data** (128,503/128,503 joined, 0 mismatches; `-1`-code and duplicate-ID guards added and verified) |
| Gene-symbol mapping boundary (is_symbol=True path) | **tested against the REAL installed gene table** (no model weights touched) — 19,194/24,737 genes map, clears the library's own 10,000-gene floor |
| Built-in `min_nnz_genes=200` cell filter, for THIS object | **tested against the REAL object's real sparsity structure** (bounded, chunked, login-node-safe) — 0/128,503 cells below the threshold |
| `scprint2 embed` CLI itself | verified present and responsive (live `--help`/`embed --help`, 2026-10-16) |
| scPRINT-2 small-v2.ckpt on a *different* real human cohort | real artifact exists (`option_c/scprint2_covid_canonical.h5ad`, 2026-05-01) — NOT proof the current installed library behaves identically to whatever produced it (version cannot be confirmed) |
| Prep → embed → validate on THIS monocyte input | **ready, not_run** — bundled job script written, not submitted |
| Model/checkpoint | **small-v2.ckpt is the confirmed FIRST TECHNICAL-PILOT model** (owner "ok"); final manuscript backbone still unconfirmed, UCE not excluded |

## 1. Entry points

**Existing, reused:** `model_registry.SCPRINT2_SMALL_V2` (checkpoint path/hash/species
evidence). **`run_scprint2_human_covid.sh` is no longer treated as a behavioral reference for
the installed CLI's flags** — v1's mistake was inferring the current library's behavior from
an old artifact's existence; this version verifies claims against the ACTUAL installed source
(`scprint2/cli.py`, `scdataloader/preprocess.py`) instead.

**Written this round:**
- `results/stephenson_monocyte_scprint2_pilot_prep.py` — GEX-only subset (drops the 192
  Antibody-Capture features) + Gate 1/DecontX annotation. Fixed this round: the overwrite
  guard now also checks for an existing `_summary.json` (an independent fixture proved it did
  not before); `validate_prepared_adata()` now re-derives every check directly from the actual
  objects passed in (cell order/identity, GEX var correspondence, exact original-count
  preservation via a sparse diff, QC-annotation completeness and valid range) instead of
  trusting a caller-supplied report.
- `results/stephenson_monocyte_scprint2_pilot_validate.py` — **new this round**: checks the
  REAL `embed` output for a recognized embedding key, zero cells, non-finite values, and
  reports (without failing on) a cell-count difference from the input and missing QC columns.
- `results/stephenson_qc_source_crosscheck.py` — fixed this round: categorical-code validation
  (the `-1` defect), duplicate-cell_id rejection, strict boolean parsing.
- `scripts/foundation_models/scPRINT/run_scprint2_stephenson_monocyte_pilot.sh` — rewritten to
  bundle prep → embed → validate into one fail-stop job, `--no-requeue`, a runtime checkpoint
  SHA256 check, correct output-path handling, and `--is_symbol true`. **Written, not
  submitted.**

## 2. Real verification of the two preprocessing boundaries (not guesses)

Both probes below were run directly against the actually-installed library and the actual
object — bounded, read-only, no model weights touched, safe on the login node (the second one
is explicitly chunked to avoid materializing the 5.7 GB object):

**Gene-symbol mapping** (`scdataloader.utils.load_genes("NCBITaxon:9606")`, the real function
the installed preprocessor calls): returned a real 85,720-row, Ensembl-indexed gene table with
a `symbol` column. Of this object's 24,737 Gene-Expression var names, **19,194 match a `symbol`
value exactly** (5,543 do not) — comfortably over the installed `Preprocessor`'s own
`min_valid_genes_id=10_000` floor, so `--is_symbol true` will not trigger "Dataset dropped due
to too many genes not mapping to it." The 85,720-row table size also exactly matches the one
proven real scPRINT-2 output's own var count (`option_c/scprint2_covid_canonical.h5ad`: 85,720
genes) — independent confirmation this is the actual model vocabulary, not a coincidence.

**Built-in cell filter** (`sc.pp.filter_cells(min_genes=200)`, unconditional in the installed
preprocessor): read `X`'s `indptr`/`indices` directly via `h5py` in 10,000-row chunks (never
loading the full ~813 MB `indices` array at once, never densifying), and counted each cell's
nonzero count within the Gene-Expression columns only. Result: **min 211, median 1,325, max
5,231 GEX-detected genes per cell — 0 of 128,503 cells fall below 200.** This is a fact about
THIS object, not a general guarantee the library never filters cells elsewhere.

**Not tested this round** (would require actually running the model): whether `--use_raw`'s
exact semantics matter for this raw-count input: left at the CLI default (`false`), unchanged,
since no evidence this round motivates a different value and changing it without evidence
would repeat the same mistake R2 corrected.

## 3. Input layer, gene axis, and preprocessing contract

| Property | Value | How confirmed |
|---|---|---|
| Source object | `02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad` | existing, hash-recorded |
| n_obs (cells) | 128,503 | bounded h5py read |
| n_vars (raw) | 24,929 (24,737 Gene Expression + 192 Antibody Capture) | bounded h5py read |
| Gene ID namespace | HGNC gene SYMBOLS | bounded h5py read of `var/_index` |
| X layer | raw integer counts stored as `float32` (CSR sparse) | bounded h5py read |
| Required preprocessing before `embed` | subset `var` to `feature_types=="Gene Expression"` | implemented, `stephenson_monocyte_scprint2_pilot_prep.py` |
| `--is_symbol` | `true` | real probe against the installed gene table, §2 |
| Built-in cell filter impact | 0/128,503 cells below the library's 200-gene floor | real, bounded, chunked probe, §2 |
| Output gene axis | the model's OWN 85,720-gene Ensembl vocabulary, not the input's | confirmed both via the proven real output file AND via `load_genes()` returning the same row count |
| Cell retention through `embed`, for THIS run | **UNVERIFIED** — not run this round; must be checked explicitly against `n_cells_in`, never assumed 1:1 (the one proven run on a different cohort was NOT 1:1: 47,977→47,114) |
| Output file naming | `<--output_filename value>_embedded.h5ad` | read directly from installed `cli.py` |

## 4. Commands (prepared, not executed)

**One bundled job** (fail-stop; prep → embed → validate; writes to a new, timestamped run
directory under `02_EMBEDDINGS/option_c_monocyte/scprint2_pilot_runs/`):
```
sbatch scripts/foundation_models/scPRINT/run_scprint2_stephenson_monocyte_pilot.sh
```
Internally, after a runtime checkpoint SHA256 check:
```
python -c "... stephenson_monocyte_scprint2_pilot_prep.main(out_h5ad=<run_dir>/...)"
scprint2 embed --adata <prepared.h5ad> --ckpt_path <small-v2.ckpt> \
    --species NCBITaxon:9606 --is_symbol true \
    --output_filename <run_dir>/stephenson_monocyte_scprint2_pilot_embeddings
python -c "... stephenson_monocyte_scprint2_pilot_validate.main(<run_dir>/..._embedded.h5ad, n_cells_in, <run_dir>/validation_report.json)"
```
Real exit code of the whole job is the validator's exit code (0 = technical pass, 1 = hard
failure) — not merely whether `scprint2 embed` itself exited 0.

**Known prerequisite gap, not fixed by this packet**: the HPCC checkout at
`/bigdata/.../302006_updategit` (which the job script's `PREP_RESULTS_DIR` points at) is still
at a stale commit (`35131cf`) with unrelated uncommitted changes, confirmed again 2026-10-16.
Before any real submission, that checkout must be updated to the authorized commit (or the
script's `PREP_RESULTS_DIR` repointed at a freshly-exported snapshot) — this is an operational
step, not a code change, and is explicitly NOT done by this packet.

## 5. Resource ceiling and basis

`--partition=short_gpu`, 1 GPU, 8 CPUs, 128G RAM, **3h** walltime (up from 2h in v1, to
accommodate the newly-bundled prep step — UNVERIFIED margin, not a measured requirement),
`--no-requeue` (R5 fix — v1's `--requeue` conflicted with this round's no-auto-retry scope).
**Basis**: the 128G/1 GPU/8-CPU spec is reused from the one proven real scPRINT-2 `embed`
job in this project; walltime is a judgment addition, not evidence-backed. **Still
UNVERIFIED for this larger input** (128,503×24,737 vs. the proven run's 47,977×13,066) — not
claimed sufficient, only the best available starting point. Fail-stop throughout
(`set -euo pipefail`); the only retry mechanism is none.

## 6. Identity and numeric acceptance criteria — now wired into the run path, not just documented

1. **Input traceability**: every source file's sha256 in the prep step's `_summary.json`.
2. **Identity/row-order preservation**: re-derived directly from the real objects by
   `validate_prepared_adata()` (R4 fix) — mutation-killed against the real production file
   this round (6 independent bypass scenarios, all now caught; see
   `tests/test_stephenson_monocyte_scprint2_pilot_prep.py`).
3. **Original counts preserved exactly** (not merely finite): `validate_prepared_adata()` now
   diffs the prepared `X` against a direct GEX-subset slice of the ORIGINAL `X` (sparse diff,
   never densified) — a changed count value is caught, not just a non-finite one.
4. **Output dimensions and numeric finiteness** (after `embed`): `stephenson_monocyte_
   scprint2_pilot_validate.py`'s `main()` is called as the job's own last step and returns a
   real nonzero exit code on a hard failure (no recognized embedding key, 0 cells, or a
   non-finite embedding value); a cell-count difference from the input and missing QC
   obs-columns are reported in `validation_report.json` but are NOT themselves hard failures
   (matching the already-observed real behavior of a different cohort's proven run).
5. **Model configuration/seed/dependency record**: `--seed_everything` default (42, confirmed
   via `--help`), checkpoint sha256 verified AT RUNTIME (R5 fix) against `model_registry.py`'s
   recorded value, not merely recorded once historically.
6. **Re-readable result and run status**: output opens with `anndata.read_h5ad`; the SLURM
   `.out`/`.err` logs plus `validation_report.json` are the run-status record.

**Explicitly NOT an acceptance criterion**: whether the embedding shows good-looking clusters
or supports any biological interpretation.

## 7. What remains for the owner

- Per-job SLURM execution authorization for the bundled job — separate from, and not granted
  by, the provisional-input-use or first-pilot-model decisions already given.
- The stale-checkout prerequisite (§4) needs resolving before any real submission.
- Whether the 3h/128G ceiling is actually sufficient — genuinely unverified; only an actual run
  (or scPRINT-2's own documented memory-scaling guidance, not checked this round) can confirm
  it, and a job that runs out of time/memory fails closed (no auto-retry), it does not corrupt
  anything.
- Final manuscript-backbone model/checkpoint confirmation remains separate and not addressed
  by this packet, per Cindy's own framing of the "ok" decision.
