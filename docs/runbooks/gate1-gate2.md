# Runbook — Upfront QC (Gate 1 / Gate 2)

Governing document: `scripts/32_rebuild_upfront_qc/qc_checklist.md` (rev.4). The checklist holds
the frozen rules and the reasoning; this runbook holds the commands and the artifact layout.

**Two-gate discipline.** Flags are computed first and nothing is deleted until the gate is
reviewed and frozen:

```
compute basic-QC flags  →  GATE 1 review → freeze → apply
                        →  Scrublet scoring (deletes nothing)
                        →  GATE 2 review → freeze → delete
                        →  DecontX refit on the QC-retained full-lineage set
```

Two standing constraints override any convenience:

- **Outcome blindness.** Thresholds are never tuned using mortality labels or predictive
  performance. Outcome-group retention is disclosure only.
- **Retention is a trigger, not a target.** A surprising retention rate prompts investigation;
  it never justifies reverse-engineering a threshold to reach a desired cell count.

---

## Paths

| What | Where |
|---|---|
| HPCC repo mirror | `/bigdata/godziklab/shared/Xinru/302006/302006_updategit` |
| Stage code | `scripts/32_rebuild_upfront_qc/code/` |
| Results root (`$R`) | `/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc` |
| Logs | `/bigdata/godziklab/shared/Xinru/302006/06_LOGS` |
| Env | `/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean` |

`$R` layout: `decontx_input/<cohort>/` (pre-QC) · `gate1_qc/<cohort>/` · `gate2_input/<cohort>/`
· `gate2_scrublet/<cohort>/`.

**Do not read `gate1_qc/{liu,wilk}/` as authoritative.** Those hold superseded generic-rule
output. Liu and Wilk's frozen results are in `gate2_input/`.

---

## Cohort layouts and rules

Two I/O layouts, because cohorts arrive differently:

| Layout | Cohorts | Script |
|---|---|---|
| `combined` (one `counts.mtx` + `cells.csv`) | Qiu, Cillo, Chongqing, Stephenson, Liu, Wilk | `gate2_run.py --layout combined` |
| `split` (`split_per_channel/<ch>/…`) | Ren, COMBAT | `gate2_run.py --layout split` |

Two cohorts additionally have frozen cohort-specific Gate 1 rules, because their input state
differs (checklist §6e). `qc_gate1_apply.py` implements only the cohort-general rules — it is
**not** the whole of Gate 1:

- **Liu** — `liu_barcode_rank_cellcall.py` → `liu_inflection_cellcall.py`. Liu ships a fixed
  30,000-barcode-per-lane demuxlet whitelist, not a cellranger cell call (verified: all 55 lanes
  have exactly 30,000 barcodes), so barcode-rank cell calling is required. Primary criterion is
  the **inflection** point; knee and valley are diagnostic only.
- **Wilk** — `wilk_seqwell_qc.py`. Seq-Well, median 441 counts/cell, so the 10x-derived
  500/200 floors are a platform-level distribution shift rather than a tail. Those floors become
  diagnostic-only columns; MAD is lower-tail on log1p counts/genes and upper-tail on mito, per
  capture × lineage with a cohort-lineage fallback below 70 cells.

`materialize_liu_wilk_for_gate2.py` writes both cohorts' frozen sets plus `prior_counts.csv`
into `gate2_input/`.

---

## Running Gate 2 — superseded commands removed

The command block that stood here invoked `--counts`, `--cells`, `--platform` and
`--prior-counts`, and called `gate2_scrublet_split.py` directly. **None of that works any
more**: those flags do not exist on `gate2_run.py`, and the split script is retired. The
current invocation is in "The safety layer" above, and it requires `--runtime-config` and
`--base-seed`.

`prior_counts.csv` is no longer produced by `materialize_liu_wilk_for_gate2.py` and must not be
passed to anything. It counted cells *entering Gate 1* under the undefined name `n_pre_qc`,
which is not the loading-level quantity the rate formula is defined over. Existing copies under
`gate2_input/{liu,wilk}/` are kept only as inputs to already-quarantined runs.

`run_gate2.sh` targets `short` (2h cap), 4 CPU, 120G. Chongqing (861K cells / 124 captures),
Ren (284 captures) and COMBAT (70 captures) do not finish in one slot and must resume.

---

## The safety layer (implemented 2026-09-21)

`gate2_run.py` replaces both old runners. It will not start without an authorized runtime config,
and none exists — that is intended, not a bug.

```bash
python3 gate2_run.py --cohort <name> --layout combined|split \
    --input <dir> --out-dir <dir> --runtime-config <csv> --base-seed <int>
```

Exit 2 = configuration refused. Exit 3 = checkpoint refused.

Per-cell `call_status` is one of `predicted_doublet` · `scrublet_noncandidate` · `too_small` ·
`threshold_failed` · `method_unreliable` · `runtime_failed`, kept separate from the per-capture
execution status. **`call_status == "predicted_doublet"` marks a MECHANICAL CANDIDATE, not a cell cleared for
removal.** It is necessary but not sufficient: cohort eligibility, an approved runtime rate, a
capture reliability criterion and a sensitivity analysis are all still required, and no
final-authorization artifact exists. A keep mask is not evidence of singlet status, and candidate calls are tri-state so "unknown" has
a representation. `method_unreliable` is reserved for a policy decision and is never assigned by
code.

Checkpoints bind 18 fingerprint fields (schema version, cohort, raw channel, ordered cell-ID
checksum, cell count, matrix fingerprint, config row checksum, rate, donor proportions and
stage, platform, base and capture seed, Scrublet parameters, software versions, git commit, code
fingerprint, call-status schema). Any mismatch is refused outright. Writes are temp + fsync +
`os.replace`, keys are hash-based so `A-B` and `A_B` no longer collide, and metadata is JSON
rather than a pickled object array.

Tests: `tests/test_gate2_safety_layer.py` — **190 passed, 0 skipped, 0 xfail, 0 xpass**, all
under temporary directories with Scrublet injected as a fake.

**Limitations of that number, stated so it is not over-read:**
- **Real `sc.pp.scrublet` executions in this suite: 0.** Interface compatibility was
  confirmed separately by one independent synthetic smoke test (1 real call, verdict
  COMPATIBLE, 300x200 synthetic integer matrix on HPCC). The two figures are about
  different things and must not be written as contradicting each other.
- One test reflects on the Scanpy signature and **skips** where Scanpy is absent (it runs on
  HPCC, where Scanpy 1.11.3 is installed).
- `adhoc_gate2_param_contract_check.py` is a **read-only diagnostic, not a test**: it compares
  `SCRUBLET_PARAMS` against the installed Scanpy's signature and defaults and does not execute
  Scrublet.

Checkpoint schema is **`gate2-checkpoint/4`** (`/3` superseded when the persisted diagnostics gained the score-distribution fields; a `/3` file would restore without them and a resumed run would disagree with a fresh one). Diagnostics carry a separate `diagnostics_schema_version`, **validated on restore** (missing, older or unknown is refused, never assumed compatible), and `gate2_diagnostics.py` is inside the code fingerprint, so editing a metric definition invalidates stored checkpoints even when the checkpoint schema itself is unchanged. `/1` and `/2` are both rejected outright, with no migration path. `/2` was superseded when the `call_status` value domain changed (`auto_singlet` → `scrublet_noncandidate`); `/1` was superseded because cell IDs were stored
in a fixed-width `U64` array, which truncated anything longer and made the capture permanently
unrestorable while the error blamed "identity/order". IDs are now stored at a width derived from
the longest present, and the save verifies them round-trip before the file is installed.

**Single-writer precondition.** Temp files are unique per call (pid + uuid4, `O_EXCL`), but the
final checkpoint is last-writer-wins and no lock arbitrates between two writers holding
different fingerprints for one capture. Running two writers against one `--out-dir` is a
precondition violation, not a supported mode.

## Historical defects — fixed, kept for the record

These are open as of 2026-09-21 and are the P0 items in `docs/execution-plan.md`.

1. **`--prior-counts` fails open.** Guarded by `os.path.exists()`; a missing file silently
   reverts the expected rate to post-QC cell counts, which the checklist forbids. Only Liu and
   Wilk have priors. The other six cohorts' Gate 2 output is quarantined as
   `gate2_scrublet/<cohort>.INVALID_postqc_prior_20260921/`.

2. **Threshold failure is recorded as success.** `scrublet_one_channel` reads
   `predicted_doublet` unconditionally after its `try/except`; scanpy warns rather than raising
   when auto-thresholding fails, so the capture gets `status='ok'` with every cell called
   singlet. Signature: `status=='ok'` and `threshold` is NaN.

3. **Checkpoints carry no provenance.** `.npz` (combined) and `DONE` (split) are keyed by
   channel name only, with no binding to code version, input checksum, or expected-rate
   configuration. A re-run after fixing (1) will silently reuse the wrong-prior results. Clear
   or version the checkpoint directory before any re-run until sidecars exist.

4. **Split-path summary is incomplete after resume.** `scrublet_filter_split.py` accumulates
   `summary_rows` only for channels processed in the current invocation; channels restored from
   `DONE` never enter it, so a resumed run's summary is not a whole-cohort summary.

5. **Gate 2 is not wired downstream.** `assemble_platelet_monocyte_benchmark_input.py` reads the
   original `*_decontx_input` directories, so Gate 2 output currently reaches nothing.

---

## Expected-rate sensitivity design — ruling of 2026-09-24

The operative decision is
`reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`.
Its follow-up `GATE2_PARAMETER_RULING.md` settles all eleven parameters, including their scopes
and invalidation conditions. The methodological mapping is complete but the production slot
remains inactive. Run only after at least one specific capture passes C1–C3, its exact config
passes audit, and the reviewer records the delegated small-pilot authorization. The first pilot
is one complete capture, never a cell subsample; no full-cohort or deletion permit is implied.
Use the fixed prior stress grid `0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20` only after
the input/parameter conditions and execution gates are satisfied. These are assumed algorithmic
settings, not prevalence estimates, bounds, removal quotas or vendor-derived residual rates.

The former instructions to cap `0.008 * n_recovered / 1000`, substitute observed post-demux
donor shares, and reuse a 10x rate for Wilk are withdrawn. Chemistry-specific vendor mappings
do not supply the unknown selection from recovery through upstream filtering and Gate 1.
Post-demux shares describe survivors, with no identified direction of bias for loading shares.

For explanation only, a two-cell selection model gives
`D = r * (s*a_W + (1-s)*a_X)` and `r_res = D / ((1-r)*a_S + D)`, where the `a` terms are class-
specific survival probabilities. The familiar `r*s / ((1-r)+r*s)` requires full cross-source
removal and no singlet/within-source loss. Those conditions are not measured here; the formula
does not generate runtime parameters.

Qiu/Chongqing are conditional standard-scoring candidates. Cillo/Liu/COMBAT support residual-
input diagnostics, not calibrated same-donor rates; Ren remains residual-scoring-only.
Stephenson and Wilk remain ineligible on present evidence. Wilk cannot borrow a 10x prevalence
model; this does not imply that all Seq-Well expression-based diagnostics are impossible.

Freeze simulations, preprocessing, PCA and neighbours within each capture/seed across priors.
Rank/AUC invariance is expected under the monotone score transformation. Report threshold
failures and call-set sensitivity without selecting a favourable rate. Missing or failed calls
remain unknown; no cell deletion is authorized by this design.

---

## MAD specification

Unscaled (`scale=1.0`), 5 MADs, computed per capture. Low and high tails are flagged separately
and only high tails delete (cohort-general rules; Wilk inverts this deliberately). `MAD == 0`
flags nothing, NaN is never flagged, and captures under 70 cells get no MAD at all.
