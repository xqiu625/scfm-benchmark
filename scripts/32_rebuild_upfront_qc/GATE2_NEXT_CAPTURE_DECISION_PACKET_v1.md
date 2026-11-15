# Gate 2 next-capture decision packet (v1)

**2026-10-12. Read-only. No job submitted, no permit created, no cell deleted, no outcome field
read, `PRODUCTION_PARAMETER_RULING` untouched, no Harmony/correction-method choice made.**

Responds to `review:20261012_001137` next-step items 4–6. Prepares, but does not act on,
a concrete candidate for extending Gate 2 scoring beyond the single accepted `GSM5102900_HC1`
pilot. Selection below is based on provenance and technical readiness only — never on scoring
outcome, "nicer" results, or anything from `GSM5102900_HC1_PILOT_INTERPRETATION_v1.md`.

**Revised in place (2026-10-12) per `review:20261012_124351` (F2, F4).** Fixes: the
literal-name-match count was internally inconsistent (§1 said "6 of 11" while listing 5); C1/C2
for `GSM5102901_HC2` were written as PASS on an inherited, series-level citation rather than
`GSM5102901`'s own sample-level evidence — now replaced with GSM5102901's own GEO SOFT record,
read from the already-cached `family.soft.gz` on HPCC, not inherited from `GSM5102900_HC1`'s
audit or substituted with a downstream paper; §4's HC1 feasibility line used stale
pre-preflight audit wording ("CANNOT BE CONFIRMED") that ignored the later-completed HC1 preflight
job, which did measure and confirm feasibility — corrected below.

## 1. Remaining Qiu captures — identity and readiness

Qiu (GSE167363) has 12 captures total; `GSM5102900_HC1` is the one already scored. Identity
columns (`gsm`, `sample_name`, `donor`, `disease_state`, `draw_time`) are from
`scripts/32_rebuild_upfront_qc/channel_inventories/qiu_gse167363_channel_inventory.csv`.
`n_called_pre_qc`/`n_gate1_retained` are from the already-committed, bug-fixed crosswalk
`reports/manuscript/revision/literature_audit/T7_doublet_prior_evidence/gate2_source_crosswalk.csv`
(commit `52c9ba1`, 2026-09-21; produced by `gate2_source_crosswalk.py`, which superseded
`gate2_prior_inventory.py`'s three join defects) — **not reverified this round**, cited as-is.
The container name-match column is **independently recomputed this round** (method in §3) by
reading the canonical container's own `obs/Sample_ID` values and checking whether the deposit's
GEO filename token appears there literally.

**The `disease_state` column below is outcome-adjacent metadata from the cited inventory, shown
for identity/complexity description only** — see §2's disclosure for exactly what this round did
and did not do with it (no h5ad outcome column was read; no candidate is ranked by outcome label).

| GSM | deposit filename token | sample_name (GEO) | donor | disease_state | draw | deposit barcodes (pre-QC) | Gate 1 retained | container name match |
|---|---|---|---|---|---|---:|---:|---|
| GSM5102900 | `HC1` | HC1 | Tissue donor 1 | healthy control | — | 7,666 | 7,453 | **literal** (`HC1`) — **scored (this pilot)** |
| GSM5102901 | `HC2` | HC2 | Tissue donor 2 | healthy control | — | 8,345 | 8,185 | **literal** (`HC2`) — confirmed §3 |
| GSM5102902 | `P25_T0` | NS LS T0 | Tissue donor 3 | non-survived septic | T0 | 4,106 | 2,991 | **ambiguous** — no `P25` token in container; candidate match `NS1_T0`, unconfirmed |
| GSM5102903 | `P25_T6` | NS LS T6 | Tissue donor 3 | non-survived septic | T6 | 4,622 | 3,814 | **ambiguous** — candidate match `NS1_T6`, unconfirmed |
| GSM5102904 | `P50_T0` | S1 T0 | Tissue donor 4 | survived septic | T0 | 3,248 | 2,641 | **ambiguous** — no `P50` token in container; candidate match `S1_T0`, unconfirmed |
| GSM5102905 | `P50_T6` | S1 T6 | Tissue donor 4 | survived septic | T6 | 4,566 | 4,167 | **ambiguous** — candidate match `S1_T6`, unconfirmed |
| GSM5511351 | `NSES_T0` | NS ES_T0 | donor 5 | non-survived septic | T0 | 4,286 | 3,999 | **ambiguous** — no `NSES` token in container; candidate match `NS2_T0`, unconfirmed |
| GSM5511352 | `NSES_T6` | NS ES_T6 | donor 5 | non-survived septic | T6 | 6,257 | 5,868 | **ambiguous** — candidate match `NS2_T6`, unconfirmed |
| GSM5511353 | `S2_T0` | S2_T0 | donor 6 | survived septic | T0 | 6,963 | 6,632 | **literal** (`S2_T0`) |
| GSM5511354 | `S2_T6` | S2_T6 | donor 6 | survived septic | T6 | 5,682 | 5,412 | **literal** (`S2_T6`) |
| GSM5511355 | `S3_T0` | S3_T0 | donor 7 | survived septic | T0 | 2,853 | 2,816 | **literal** (`S3_T0`) |
| GSM5511356 | `S3_T6` | S3_T6 | donor 7 | survived septic | T6 | 5,650 | 5,500 | **literal** (`S3_T6`) |

Across all **12** Qiu captures, exactly **6** have a **literal, inference-free** name match
between the deposit's GEO filename token and the canonical container's own `Sample_ID` value
(`HC1`, `HC2`, `S2_T0`, `S2_T6`, `S3_T0`, `S3_T6`) — exactly reproducing `CAPTURE_C1_C3_AUDIT.md`'s
"6 of Qiu's 12 GEO names do" statement. `HC1` is one of those 6 and is already scored, so of the
**11 unscored** captures, **5** are literal matches (`HC2`, `S2_T0`, `S2_T6`, `S3_T0`, `S3_T6`) and
**6** are not (`P25_T0/T6`, `P50_T0/T6`, `NSES_T0/T6`), each of which would require an unverified
naming inference (`P25`→`NS1`, `P50`→`S1`, or `NSES`→`NS2`) the original audit explicitly declined
to make for `GSM5102900_HC1`'s own identity close; this packet does not attempt it either.

## 2. Candidate selection

**Proposed: `GSM5102901_HC2`** — same cohort as `GSM5102900_HC1`, another independent donor/sample
(not a "twin" of the same donor; `GSM5102901` is Tissue donor 2, a separate person from
`GSM5102900`'s Tissue donor 1).

**Outcome-boundary disclosure, stated plainly rather than only implied by scope statements
elsewhere in this packet**: the candidate-evidence script (`gate2_hc2_candidate_evidence.py`) did
not read any h5ad outcome column — confirmed, it opens only `obs/_index` and `obs/Sample_ID`. The
channel inventory `§1` cites (`qiu_gse167363_channel_inventory.csv`) **does** carry
outcome-bearing metadata (`disease_state`: healthy control / non-survived septic / survived
septic), because that is the inventory's own identity-mapping content, not something this round
added. Selection below uses **only** verifiable provenance readiness, literal container-name
match, and capture/draw structural complexity (donor count, timepoint pairing) — **not**
survived/non-survived status, and no candidate is ranked or excluded by outcome label. A prior
version of this section used the phrase "septic-outcome proximity" as part of the stated
rationale for preferring `HC2`; that phrasing is retracted here as imprecise — the actual,
outcome-free reason to prefer `HC2` over `S2`/`S3` is the structural one below.

Reasons, all provenance/readiness-based:

- **Literal container name match**, independently confirmed (§3) — no naming inference needed,
  the same property that made `GSM5102900_HC1` "the case most likely to close" per
  `CAPTURE_C1_C3_AUDIT.md`.
- **Same donor/draw structure as the already-accepted pilot**: one donor, one draw, no timepoint
  pairing. This extends the existing technical acceptance to a second independent draw without
  introducing a new structural complexity class (a second capture per donor, paired across two
  timepoints) that `S2_T0`/`S2_T6`/`S3_T0`/`S3_T6` would add — those four are also literal
  matches and are listed as the next-best alternatives, but each is one of a same-donor pair
  across two draws, which is more moving parts for a capture whose sole purpose here is extending
  technical acceptance to one more independent sample, not analysing within-donor timepoint
  structure.
- **Qiu-cohort reuse**: `build_qiu_decontx_input.py` and the already-reviewed C1–C3 methodology
  for Qiu's per-GSM CellRanger triplets apply unchanged; no new cohort-level loader or join logic
  is needed.

## 3. C1–C3 evidence for `GSM5102901_HC2` — independently gathered this round

**Status, corrected: file-structure checks are recorded below; full C1/C2 provenance is now
closed using `GSM5102901`'s own GEO sample record (not a prior-version substitute), C3 passes on
the deposit as with `GSM5102900_HC1`.** A prior version of this section wrote C1/C2 as PASS by
citing `CAPTURE_C1_C3_AUDIT.md`'s HC1-specific, series-level text (depositor processing statement,
chemistry) as if it covered `GSM5102901` directly — it does not, by itself, establish
`GSM5102901`'s own processing history or GEM-well physical-capture unit. This round located
`GSM5102901`'s own sample record in the GEO series SOFT file already cached on HPCC
(`/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/family.soft.gz`, extracted
read-only via `zcat | awk`) rather than relying on the live GEO webpage (attempted; returned a bot
verification page, consistent with `review:20261012_124351`'s own attempt) or a
downstream reanalysis paper.

Evidence gathered read-only on HPCC via a new, unreviewed, ad-hoc script
(`gate2_hc2_candidate_evidence.py`, saved under
`gate2_pilot_interpretation_evidence_20261012/` with this packet — **not** a modification of the
frozen, already-reviewed `verify_capture_evidence.py`, which remains hardcoded to `GSM5102900_HC1`
and untouched). Opens only the three deposit files and the two named container datasets
(`obs/_index`, `obs/Sample_ID`); no outcome column read. Fixed this round (F5): refuses to
overwrite an existing non-empty output directory; decodes the container's `Sample_ID` categorical
with an explicit out-of-range/missing-code check (`c3_sample_id_decode` in the JSON: 0 missing,
0 out-of-range, of 413 categories); checks the matrix header's declared `n_cells`/`n_genes`
against the actual barcode/feature file counts (both matched, `dims_consistency` in the JSON);
checks for duplicate/malformed extracted container tokens (0 of either, of 82).

**Input paths and checksums:**

```
/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted/
  GSM5102901_HC2_barcodes.tsv.gz      42,960 bytes
    sha256 3113aede3fcf401d613d02ef27ab6ab3b05b568b123e21ea57cf19d1950ee644
  GSM5102901_HC2_features.tsv.gz     304,728 bytes
    sha256 0774e088c47a7cd2cb99a22c7d41c50405940ec0f84682c6f1f6a7652cd23587
  GSM5102901_HC2_matrix.mtx.gz    39,369,513 bytes
    sha256 8181ebe77603c8c607ebcf9264f18c43b8de4bf8e63bce2efb5e1768054d8fe5
```

**`GSM5102901`'s own GEO sample record** (from `family.soft.gz`, the `^SAMPLE = GSM5102901` block,
read this round — not HC1's record, not a downstream paper):

```
!Sample_title = HC2
!Sample_source_name_ch1 = PBMCs from non-sepsis donor 2.
!Sample_characteristics_ch1 = donor: Tissue donor 2 | disease state: healthy control
!Sample_characteristics_ch1 = Sex: Female | age range: 45-50
!Sample_extract_protocol_ch1 = ...15,000 cells/sample single-cell suspensions were loaded into
  a Chromium Single Cell Chip (10x Genomics)...  [ONE GEM-well load, stated directly for this
  sample, not inferred from the barcode suffix alone]
!Sample_extract_protocol_ch1 = ...Chromium Next GEM Single Cell 3ʹ Reagent Kits v3.1...
  sequenced at UCSD Genomics center on the NovaSeq platform
!Sample_data_processing = The Cell Ranger software pipeline (v3.1.0, 10x Genomics) was used to
  demultiplex cellular barcodes and map reads to the human reference genome
  (refdata-cellranger-GRCh38-3.0.0) (command cellranger count).
!Sample_data_processing = We obtained the feature-barcode matrix that contains gene expression
  counts in the outputs.
!Sample_library_strategy = RNA-Seq
!Sample_library_source = transcriptomic
!Sample_supplementary_file_1/2/3 = .../GSM5102901_HC2_{barcodes,features,matrix}...  [matches the
  three deposit files' names exactly]
```

**C1 — un-normalised raw-count provenance:**

| evidence | value |
|---|---|
| MatrixMarket header | `%%MatrixMarket matrix coordinate integer general` |
| dimensions / nnz | 33,538 genes × 8,345 cells, 10,024,339 non-zeros (streamed and counted; declared = observed; also cross-checked against the barcode/feature file row counts, `dims_consistency` in the JSON) |
| values | integral, min 1.0, 0 non-integer of 10,024,339 |
| `features.tsv.gz` feature types | `{'Gene Expression'}` only — no Antibody Capture |
| sample-level source statement | `GSM5102901`'s own `!Sample_data_processing`/`!Sample_library_strategy`/`!Sample_library_source` above — genuine for this sample, not inherited |

**PASS** — on `GSM5102901`'s own sample-level GEO record plus the file-level checks, not on an
inherited series-level citation.

**C2 — physical capture and raw-barcode namespace:**

| evidence | value |
|---|---|
| barcodes | 8,345, all unique |
| token length | 16 nt for every barcode |
| GEM-well suffix | `-1` for all 8,345 |
| sample-level loading statement | `GSM5102901`'s own extract protocol: "15,000 cells/sample single-cell suspensions were loaded into a Chromium Single Cell Chip" — one GEM-well load, stated for this sample directly |

**PASS** — the uniform `-1` suffix is now corroborated by `GSM5102901`'s own stated single-chip
loading, not file-count/suffix inference alone.

**Recorded, not resolved** (same caveat `CAPTURE_C1_C3_AUDIT.md` records for `GSM5102900_HC1`,
now independently confirmed to apply to `GSM5102901` too, from its own record above): the deposit
states 3′ v3.1 chemistry; this blocks any implementation depending on a vendor rate, not a pure
algorithmic-stress design.

**C3 — complete declared-stage population and parent-state support:**

First computation attempt used a naive barcode-token split (`token = barcode.split('-')[0]`) on
both sides of the join and found **0/82** container cells matching the deposit — because the
container's cell index is formatted `<Sample_ID>_<16nt-barcode>` (e.g.
`HC2_AACCACAGTTGTGTAC`), not a bare token, so the split-based extraction left the `HC2_` prefix
attached on the container side and never matched. Re-run using the same barcode-extraction regex
`verify_capture_evidence.py` uses for `GSM5102900_HC1` (`re.compile(r'[ACGTN]{8,}')` against the
container's `_index`, matched against the deposit's token set):

```
container cells, Sample_ID="HC2"      82
deposit barcodes                   8,345
container tokens found in deposit   82/82   (100.0%)
container is a subset of the deposit   True
fraction of the deposit retained       0.9826%
```

Same pattern as `GSM5102900_HC1`'s 65/65 (0.85%): the canonical container holds a small,
lineage-selected (platelet) slice of the full CellRanger-called deposit. **C3 PASS for the
deposit, FAIL for the container** — a future scoring experiment on this capture would, as with
`GSM5102900_HC1`, have to take the deposit's 8,345-cell matrix as input, not the 82-cell container
slice. 82/82 container cells matching the deposit establishes the identity join is clean; it does
not, by itself, prove the deposit's full 8,345-cell population satisfies C3's "complete
declared-stage, not lineage-preselected" requirement — that rests on the same reasoning
`CAPTURE_C1_C3_AUDIT.md` gives for `GSM5102900_HC1` (the deposit is the CellRanger-called output
at the declared stage, the container is a lineage-selected slice of it), which applies here by the
same structure, not by inheritance of `GSM5102900_HC1`'s specific numbers.

**Verdict: C1–C3 are closable for `GSM5102901_HC2`**, now on `GSM5102901`'s own sample-level GEO
record plus independently-gathered file-level evidence — not an inherited citation. This is a
candidate-evidence finding only, not a new formal audit sign-off with the same rigor as
`CAPTURE_C1_C3_AUDIT.md` (e.g. no full HVG/PC feasibility measurement yet, §4) — see §5 for what
remains before any permit could exist.

## 4. What is NOT yet established for this candidate

Per `review:20261012_001137`'s explicit instruction, none of `GSM5102900_HC1`'s own
measured values (7,666 cells, `n_neighbors=44`, graph shape `[22998,131]`) are assumed to carry
over. For `GSM5102901_HC2`, all of the following are **unmeasured**:

- **`n_cells_effective`/`n_usable_features` strict-PC feasibility — unmeasured for `HC2`,
  corrected for `HC1`.** A prior version of this section said `GSM5102900_HC1`'s own feasibility
  was still `CANNOT BE CONFIRMED` — that was the state at `CAPTURE_C1_C3_AUDIT.md`'s 2026-09-24
  writing, **before** `GSM5102900_HC1`'s own score-free preflight (job `29073092`,
  `preflight_commit=9aab7bef...`, the same commit the real scoring manifest records as
  `preflight_commit`) completed and measured it. Read directly from that preflight's own manifest
  this round (`CAPTURE_PREFLIGHT_MANIFEST.json`, `/bigdata/.../gate2_capture_preflight/9aab7bef_20260924_225552/`):
  `status=complete`, `n_cells_effective=7666`, `n_usable_features_hvg=1576`,
  `feasible=True` (`30 < min(7666,1576)=1576`), `n_neighbors=44`, `k_adj=132` — all four of these
  numbers match what the real scoring job later actually used. **`GSM5102900_HC1`'s feasibility is
  measured and confirmed, not unresolved**; only `GSM5102901_HC2`'s is unmeasured, and none of
  `HC1`'s specific numbers (7,666 cells, 1,576 HVG, `k=44`) are assumed to carry over to `HC2`.
- **`n_neighbors` (`k`) for `GSM5102901_HC2`.** The approved formula
  `k = round(0.5*sqrt(post-filter N))` is already a cohort-general ruling
  (`GATE2_PARAMETER_RULING.md`), not specific to `GSM5102900_HC1`, so it would carry over
  unchanged — but the `N` it is evaluated against (post `min_genes=3` filter for *this* capture)
  has not been measured, so no numeric `k` value is asserted for `HC2` here.
  `min_genes_per_cell`/`cells_passing_min_genes_3`/`genes_passing_min_cells_3`/`n_usable_features_hvg`,
  all now confirmed measured for `GSM5102900_HC1` (above), are not yet measured for
  `GSM5102901_HC2`.
- **Graph construction shape**, which follows directly from `N` and `k` once both are measured for
  `HC2` — not computable yet.
- **A score-free preflight is needed for `HC2` specifically**, mirroring `GSM5102900_HC1`'s own
  completed preflight (`capture_preflight_GSM5102900_HC1.py`, job `29073092`, now confirmed
  `status=complete` above) before any scoring config could be written. No such preflight exists
  yet for `GSM5102901_HC2` — §5 below prepares a candidate implementation for it.
- **Config and entry-point generalization.** `gate2_scoring_config_GSM5102900_HC1.py` and
  `gate2_score_capture_GSM5102900_HC1.py` are both capture-specific by name and by hardcoded
  identity constants (`CAPTURE`, `INPUT_SHA256`, `config_sha256`, etc.) — they are not
  parameterized over capture. Scoring `GSM5102901_HC2` would need either a new
  `gate2_scoring_config_GSM5102901_HC2.py`/`gate2_score_capture_GSM5102901_HC2.py` pair following
  the same pattern, or a generalized single config/entry-point pair taking the capture as a
  parameter. Either way this is **new, unreviewed implementation** requiring its own scoped
  methodology/code review before any permit could reference it — the existing
  `implementation_commit=07b5662`/`implementation_sha256` pin is specific to the `HC1` entry point
  and does not cover a hypothetical `HC2` one.

## 5. `GSM5102901_HC2` score-free preflight — candidate implementation, NOT submitted

Per `review:20261012_124351` item B (5–7): a reviewable candidate implementation, built
and small-fixture-tested this round, but **not run against the real matrix and not submitted**.

**Files** (new directory, HC1's own preflight files untouched):

```
scripts/32_rebuild_upfront_qc/hc2_preflight/capture_preflight_GSM5102901_HC2.py
scripts/32_rebuild_upfront_qc/hc2_preflight/run_capture_preflight_HC2.sh
scripts/32_rebuild_upfront_qc/hc2_preflight/test_capture_preflight_hc2.py
```

**Candidate commit**: the commit that lands these files on `main` — resolve with
`git log -1 --format=%H -- scripts/32_rebuild_upfront_qc/hc2_preflight/` rather than a hash
hardcoded into this document, so this reference cannot go stale if the commit history around it
changes.

**Structural relationship to the reviewed `GSM5102900_HC1` preflight**
(`reports/manuscript/revision/gate2_capture_audit_2026-09-24/capture_preflight_GSM5102900_HC1.py`):
identical Scrublet pipeline call sequence (filter → normalize → HVG → simulate doublets →
normalize → `zscore`/`pca` manifold, stopping at `set_manifold`), identical score-free boundary,
identical imports from `gate2_parameter_ruling`/`gate2_stress_design`/`gate2_stress_runner`
(unchanged, generic across captures — confirmed this round: `grep` of both modules for
`HC1`/`GSM5102900`/hardcoding returned nothing). Differences: `GSM5102901_HC2`'s own 3 input
hashes (from §3, not `GSM5102900_HC1`'s); the relative path to `scripts/32_rebuild_upfront_qc/code`
is corrected for this file's new location (`../code`, one level, not the original's four — the
exact bug a naive copy-paste would have introduced); no measured value
(`n_cells_effective`/`n_usable_features_hvg`/`n_neighbors`/`k_adj`/`manifold` shape) is pinned in
advance — all are computed fresh from `GSM5102901_HC2`'s own data when the script runs, exactly as
the original computes them fresh for `GSM5102900_HC1`.

**Score-free boundary, tested not just asserted.** The script never calls or imports
`Scrublet.calculate_doublet_scores`/`Scrublet.call_doublets`. 8 small-fixture tests
(`test_capture_preflight_hc2.py`, tens of synthetic cells/genes — never the real 8,345-cell
matrix, never the real 5 GB container) ran on HPCC inside a directory mirroring the real repo's
relative layout (so the corrected import path is exercised for real, not assumed): input-identity
refusal on a missing file, a hash mismatch, a byte-count mismatch, and a container mismatch;
output-directory non-overwrite refusal; the infeasible-stop path (tiny 20-gene fixture, deterministically
below the `n_prin_comps=30` floor regardless of random variance); and the feasible-complete path
with `Scrublet.calculate_doublet_scores`/`call_doublets` monkeypatched to raise if ever called.
**All 8 passed.** Per DEVELOPMENT.md rule 2, a mutation (inserting a call to
`scrub.calculate_doublet_scores()` right before the script marks itself `status='complete'`) was
then injected into a scratch copy on HPCC, and the boundary test was confirmed to **fail** and
report exactly that call — proving the test can catch the violation it claims to guard, not merely
pass by construction. The original file was restored and the full suite re-confirmed clean
(`grep -c "calculate_doublet_scores()"` = 0 in the restored file; 8/8 passing again).

**Resource ceiling and basis**: `short`, `2:00:00`, `48G`, 4 CPUs — copied from
`GSM5102900_HC1`'s own preflight (`run_capture_preflight.sh`), which ran the structurally
identical score-free pipeline to completion well inside this ceiling on a comparable-scale input
(job `29073092`). This is a starting basis stated explicitly (per the review's own instruction not
to let a citation substitute for a stated one), not a guarantee for `GSM5102901_HC2`'s 8,345 cells
(vs. `HC1`'s 7,666) — if this job is ever authorised and run, its own actual resource usage should
be recorded and compared, not assumed to match.

**Acceptance conditions for this candidate preflight, if it is ever run**: `status="complete"`
(not `"stopped_infeasible"`, `"stopped_invalid_manifold"`, or an unhandled exception);
`produced_scores_calls_or_rates=False`; `outcome_columns_read=[]`;
`production_parameter_ruling_written=False`; `inputs` in the manifest match the three
`GSM5102901_HC2`-specific hashes in §3 exactly; `manifold.all_finite=True` and
`manifold.observed_shape_matches_expected=True`/`manifold.simulated_shape_matches_expected=True`
(revised 2026-10-12 per `review:20261012_145004` F2 — the old field name
`rank_supports_n_prin_comps` was renamed because it only ever checked declared shape, never
numerical rank; see `capture_preflight_GSM5102901_HC2.py`'s R1 fix); `neighbor_graph.verified=False`
recorded explicitly (the graph is never built by a score-free preflight — this is expected, not a
defect). The wrapper itself (the same gate `run_capture_preflight.sh` uses for `GSM5102900_HC1`)
enforces, at submission time, that `SCRIPT`/`WRAPPER`/all 4 `IMPORTED` files specifically (not a
whole-working-tree guarantee) are tracked and have no staged or unstaged changes; the manifest's
separate `working_tree_clean` field is a different, strictly stronger whole-repo check and should
not be read as the same guarantee the wrapper enforces.

**Full, un-executed sbatch command** (not run this round):

```bash
mkdir -p logs   # SLURM opens --output at submission time; must exist first
sbatch scripts/32_rebuild_upfront_qc/hc2_preflight/run_capture_preflight_HC2.sh
```

This round ran only the small-fixture test suite above, directly on HPCC (not via `sbatch`, not on
the login node against anything beyond KB-scale synthetic data) — never the real deposit, never
the real container, never `sbatch`.

## 6. Two options and consequences

**Option A — authorize the `GSM5102901_HC2` score-free preflight (§5), then continue toward a
scoring config/entry-point review.** The preflight step itself is read-only, scoring-free, and
produces exactly the measurements §4 lists as missing (`n_cells_effective`, `n_usable_features_hvg`,
`n_neighbors`, `k_adj`, manifold identity) — it commits to nothing beyond that. A further,
separate decision would still be needed to write and review a scoring config/entry point, and
another to permit and run the actual scoring job. Consequence: extends Gate 2's technical
acceptance toward a second capture under the same rigor as `GSM5102900_HC1`, at the cost of
repeating the full preflight→review→permit→dry-run→scoring cycle this pilot just completed.

**Option B — pause and consolidate the threshold-sensitivity evidence first.** Treat
`GSM5102900_HC1_PILOT_INTERPRETATION_v1.md`'s finding (candidate count and threshold both
non-monotonic in the prior; 4 distinct candidate sets, not 8) as something worth absorbing into
the eventual Gate 2 production ruling or prior-selection discussion before investing further in a
second capture. Consequence: no new review burden now, but the single-capture evidence base does
not grow, and whether the non-monotonic pattern is capture-specific or general remains unknown
until a second capture is scored. The candidate preflight implementation in §5 stays ready, unrun,
for whenever this is revisited.

**Recommendation: Option A's first step only** — authorizing just the preflight job in §5, which
is itself read-only and scoring-free, the same low-risk first step `GSM5102900_HC1` took before
any config existed. Writing and reviewing the `HC2` config/entry point and any further job remain
separate, later decisions this packet does not request. **The final choice between A and B, and
whether/when to authorize the preflight job itself, is the owner's.**

This packet creates no permit, submits no job, activates no production ruling, changes no
deletion policy or DecontX unit, and makes no Harmony/correction-method choice.

## 7. Scope and verification boundary

All HPCC reads this round were read-only (deposit files, two named container datasets, the
already-cached GEO SOFT file, the completed `GSM5102900_HC1` preflight manifest) or small-fixture
test execution (tens of synthetic cells/genes, HPCC only). No outcome column read. No SLURM job
submitted. No cell deleted. `gate2_hc2_candidate_evidence.py` is evidence-gathering only — not
part of `IMPLEMENTATION_FILES`, not referenced by any permit, not run against the scoring entry
point. `scripts/32_rebuild_upfront_qc/hc2_preflight/` is a candidate implementation, reviewable
but unauthorized — not deployed via the push-free bundle process real permits use, not run against
real data, not submitted.
