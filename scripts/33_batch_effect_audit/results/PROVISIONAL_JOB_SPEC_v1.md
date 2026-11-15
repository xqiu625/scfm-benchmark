# Provisional batch-effect diagnostic job — specification for approval

**Status: specification only. No SLURM job has been submitted.** Per `CLAUDE.md`: "No SLURM
submission without explicit authorization, per job." This document exists so that authorization
can be requested and given against a concrete, already-partially-verified plan, not an abstract
one.

## What this job is, and is not

It runs the stage-33 diagnostic CLI (`run_batch_diagnostics.py`, synthetic-fixture-validated this
session — see `scripts/33_batch_effect_audit/results/synthetic_fixture_v2/`) against a **bounded,
patient-balanced subsample of the real, existing, pre-Gate-2 assembled benchmark input**. It is
explicitly labeled `pre-Gate-2 PROVISIONAL, not final` throughout — Gate 2 (doublet scoring) is
still blocked (`docs/current-status.md`), so this cannot be, and does not claim to be, a
final-benchmark result. It answers one question: does the already-tested diagnostic pipeline run
correctly end-to-end on real (if provisional) data, and what does it show on that data — not
whether a batch effect "is present" in some final, citable sense.

## Step 1 — input preparation (new, bounded, streaming; verified up to the point below)

**Command:**
```
python scripts/33_batch_effect_audit/code/prepare_provisional_diagnostic_input.py \
    --input-root /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2 \
    --out /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v1/input.h5ad \
    --n-per-donor 20 --seed 0
```

**Code identity:** commit `<filled in at the commit that accompanies this spec>` (this file is
committed in the same commit as the script).

**Real, verified input paths and sizes** (not invented — checked directly this session):
- `.../platelet_monocyte_benchmark_input_v2/raw/cells.csv` — 208 MB, 1,236,881 rows.
- `.../platelet_monocyte_benchmark_input_v2/raw/counts.mtx` — **22.2 GB** (MatrixMarket, genes × cells).
- `.../platelet_monocyte_benchmark_input_v2/decontx/counts.mtx` — **50.1 GB**.
- `.../platelet_monocyte_benchmark_input_v2/genes.csv` — 11,709 genes (11,710 lines incl. header).

**Real, verified estimated output size** — from an actual `--dry-run` execution against these
exact files this session (reads only `cells.csv`, never touches the two `counts.mtx` files; ~10s
wall time on the login node):

```
1,236,881 / 1,236,881 cells match lineages [Platelet, Monocyte_Myeloid] (already lineage-restricted)
Patient-balanced cap (<=20 cells/donor, seed=0): 665 donors -> 13,232 cells
```

Per-cohort/lineage breakdown (the real numbers, not estimated):

| Cohort | Monocyte_Myeloid | Platelet |
|---|---:|---:|
| COMBAT 2022 | 2,396 | 43 |
| Chongqing sepsis 2025 | 1,794 | 686 |
| Cillo 2021 | 566 | 14 |
| Liu 2021 | 622 | 38 |
| Qiu 2021 | 106 | 34 |
| Ren 2021 | 3,428 | 481 |
| Stephenson 2021 | 1,682 | 882 |
| Wilk 2021 JEM | 392 | 68 |

Output h5ad: **13,232 cells × 11,709 genes**, two sparse layers (raw + decontx counts) — a few
hundred MB on disk, nowhere near the 22 GB/50 GB of the source files it is streamed from. The
`read_mtx_columns()` streaming reader this script reuses (already used in production by
`assemble_platelet_monocyte_benchmark_input.py` for the same files) never materializes the full
matrices; memory is bounded by the 13,232-column output, not the input file size.

**Untested part, named exactly:** the `--dry-run` execution above verified everything up to
(not including) the two streamed `counts.mtx` reads. The streaming read itself was validated on
this same code path (`read_mtx_columns`) in production by the original assembly script, and
separately smoke-tested in this session only on small synthetic matrices inside
`tests/test_batch_metadata_inventory.py`-adjacent unit tests — not yet on these exact 22 GB/50 GB
files. That is precisely the part this job authorization would exercise for real.

**Resources (Step 1):** `epyc` partition (default, 30-day cap — chosen over `short`'s 2h cap
because this is I/O-bound sequentially scanning 72 GB total from shared storage, and `short`'s
margin is too tight to be safe if the filesystem is under load that day), 4 CPUs, 32 GB memory,
**walltime ceiling 4:00:00** (requested conservatively; the actual compute is a single sequential
scan with no retry logic — if it does not finish in 4h, something is wrong, not merely slow).
Single job, no array, no checkpointing needed (one-shot, not easily resumable mid-scan, but
bounded and non-destructive — reruns are safe, nothing is deleted).

## Step 2 — diagnostic CLI on the prepared input

**Command:**
```
python scripts/33_batch_effect_audit/code/run_batch_diagnostics.py \
    --h5ad /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v1/input.h5ad \
    --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v1/diagnostic_run \
    --population-stage "pre-Gate-2 decontx_input PROVISIONAL, not final (2026-10-09)" \
    --batch-col cohort --patient-col patient_id --lineage-col lineage --outcome-col outcome \
    --run-harmony --seed 0
```

(No `--embeddings` flag: this input has no FM embeddings — 5-FM full-scale extraction has not
started, `docs/current-status.md` — so this runs expression/PCA mode only, exactly the mode this
session's synthetic validation exercised most thoroughly.)

**Resources (Step 2):** same `epyc` partition, 4 CPUs, 16 GB memory (13,232 cells is small; the
synthetic-fixture run of comparable scale — 1,500-2,400 cells/lineage after its own
patient-balanced cap — completed in well under 10 minutes per lineage on a login node), **walltime
ceiling 1:00:00**.

**Acceptance criteria:**
- Exit 0; `manifest.json` `status` is `"complete"` or `"complete_with_failures"` (not a crash).
- `manifest.json`'s `input_sha256` matches a hash taken of Step 1's actual output file (binds the
  two steps together).
- Per-lineage `n_diagnostic_cells` + `n_excluded_by_cap` sums to the real per-lineage cell count
  from the dry-run table above (1,994 Platelet-eligible + 11,238 Monocyte-eligible... i.e. the
  script's own accounting must reconcile with Step 1's known output, not just "run without
  crashing").
- `metrics_long.csv` has rows for both lineages × `expression_pca` × both scopes
  (`uncorrected`/`harmony_corrected`), i.e. 2 × 1 × 2 = 4 embedding/scope combinations minimum.
- No `population_stage` string in any output is missing the "PROVISIONAL" label.

**Failure criteria (stop and report, do not retry automatically):**
- Non-zero exit, or `manifest.json` absent.
- Any `embeddings[...].status == "failed"` for `expression_pca` specifically (an FM-embedding
  failure is tolerable since none are requested here; expression_pca failing is not).
- A cell-ID or shape mismatch between Step 1's output and what Step 2 reads.

## What this job does NOT authorize

- No cell deletion, no Gate 2 scoring, no change to any existing checkpoint or retained-cell set.
- No promotion of this run's numbers into the manuscript or `docs/current-status.md` as a final
  result — it is explicitly provisional, and the per-cohort imbalance above (COMBAT/Ren dominate
  Monocyte_Myeloid; Qiu/Cillo are tiny) means any batch-effect reading from it must account for
  wildly unequal cohort representation, not just interpret the aggregate number.
- No automatic retry if it fails; a failure gets reported, not resubmitted.
- No push and no change to `canonical_donor_labels.csv` or any other canonical artifact.
