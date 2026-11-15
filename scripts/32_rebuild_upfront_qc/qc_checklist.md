# Phase 1 QC checklist — standards audit (2026-09-21, rev. 4)

Written after the user asked, in three successive rounds, whether Scrublet / mito% + MAD /
ribosomal filtering had been applied — and the answer each time was no. Root cause: the six-phase
architecture named only "per-channel DecontX/Scrublet" as Phase 1, so standard cell-level QC was
never in the plan; per-cohort execution then focused on each cohort's own data-wrangling problem,
treating "DecontX ran clean" as the deliverable. Missing steps throw no error, so a "did it run"
standard cannot catch them — only an outside-in comparison against a domain checklist can.

**Revision history.** Rev. 1 conflated conventional thresholds, model-behaviour inferences and
project choices as if all were universal standards. Rev. 2 separated those categories and withdrew
two invalid inferences. Rev. 3 added donor-level auditing, a flags-first gate, outcome-blindness
and Scrublet/DecontX input constraints. **Rev. 4 corrects the post-demultiplexing doublet-rate
formula (rev. 3's r/n_donors was wrong), specifies MAD completely, and splits the single review
gate into the two gates the actual execution order requires.**

**A standard QC step must be *evaluated*; it does not follow that every common filter must be
*applied to delete cells*.** Whether a filter belongs depends on target-population biology and on
what downstream algorithms actually require.

Status legend: **DO** = implement · **DONE** = already handled (acceptance criteria below) ·
**RECORD** = compute and store, do not exclude on it · **NO** = deliberately not applied ·
**N/A** = not a QC-exclusion step.

---

## 0. Execution controls (rev. 3 — these govern every row below)

### 0.1 Flags-first — TWO gates, not one

Rev. 3 said "all flags are computed first", which did not match the actual execution order:
Scrublet cannot score cells until the basic QC filters have been applied, so its calls cannot be
part of the same freeze. Two separate gates:

```
basic QC flags computed (nothing deleted)
  → GATE 1: review, freeze and APPLY basic QC rules
  → Scrublet scores + candidate calls (nothing deleted)
  → GATE 2: review and freeze doublet calls
  → exclude doublets
  → DecontX
```

After Gate 2, the donor x lineage retention report is **regenerated**, because Scrublet removes a
second tranche of cells and the donor-level picture from Gate 1 is no longer current.

Neither gate deletes on its first pass. Each computes flags/scores, emits the audit tables in
§0.2, and stops. Deletion happens only after the rules are reviewed and **frozen**. Rationale: with several overlapping filters, per-filter flag counts
and cumulative removal counts diverge — a cell failing three filters is counted once in the
cumulative total but three times across per-filter tallies, and only reporting one of the two
hides which filter is actually driving attrition.

Both are therefore required, per filter and cumulatively:
- `n_flagged` per filter (independent of other filters)
- `n_removed_cumulative` after filters are combined
- the overlap structure (how many cells fail ≥2 filters, and which pairs)

### 0.2 Audit granularity — channel, donor AND lineage

Lineage-level retention alone is insufficient. **The endpoint is patient-level mortality**, so what
matters is whether a filter removes *patients*, not only cells. A 10% overall platelet loss
concentrated in a few low-yield donors can eliminate those donors from the analysis entirely while
the aggregate number looks harmless.

Required tables, for each filter and cumulatively:

| Level | Reported |
|---|---|
| Channel | cells in / flagged / retained |
| **Donor** | cells in / retained, **per lineage**; explicit list of donors dropping to **zero retained platelets**, and donors falling **below the prespecified per-donor platelet minimum** |
| Lineage (preliminary labels) | cells in / flagged / retained |

### 0.3 Outcome blindness

**QC thresholds are set blind to mortality labels and to predictive performance.** After the rules
are frozen, retention is audited by cohort and by outcome group — as a *check on whether QC
changed patient composition*, never as a tuning signal. Tuning a threshold until AUC improves, or
until an outcome group's retention looks balanced, converts QC into label leakage.

### 0.4 Retention is a trigger, not a justification

Disproportionate platelet attrition **triggers re-examination**; it does not by itself license
loosening a threshold until the retention number looks acceptable. Retaining noise to protect a
target population is as much a failure as deleting the population. Any adjustment must be argued
from evidence beyond the retention count: expression complexity, **multiple consistent** platelet
markers (a single PF4/PPBP signal is not an exemption), the presence/absence of non-platelet
signal, per-cell contamination estimates, and consistency across replicates/channels.

### 0.5 Grouping unit — which boundary MAD and Scrublet are estimated within

QC runs on the confirmable independent capture unit where one exists. Where capture information is
missing, a sample, pool or technical batch supported by the experimental design is used as a proxy
and the basis recorded. MAD thresholds are estimated within the chosen grouping; Scrublet is run,
wherever possible, on the full cell population of a single capture. Where capture boundaries
cannot be recovered, Scrublet output is treated as an approximate prediction; if the data's
heterogeneity does not support a reasonable joint analysis, it is marked **not reliably assessed**
rather than defaulting every cell to singlet. **Missing capture information is not by itself a
reason to exclude a cohort.**

| Information available | MAD grouping | Scrublet grouping | Decision |
|---|---|---|---|
| Explicit capture / channel ID | by capture | by capture | execute directly |
| Each sample/donor confirmed to be exactly one independent capture | by sample/donor | by sample/donor | execute as a trusted proxy |
| Multiple donors confirmed co-captured in one pool | by pool | by pool, retaining the major cell types | do **not** run separately on post-demultiplexing donor splits |
| Capture unclear, technical batch/library info available | finest trustworthy technical grouping | approximate, if that grouping may span several captures | record basis and uncertainty |
| Only a merged matrix, no recoverable technical grouping | compute within current scope as approximate QC | exploratory only where source and composition are similar; no forced unified calls when they differ | flag the detection limit; **undetected ≠ singlet** |

Two qualifications that govern how the table is applied:

1. **One donor or one pool may still have been loaded more than once.** Only a *confirmed* single
   capture is a trusted proxy.
2. **The expected doublet rate is never derived from a merged matrix or from post-QC cell counts.**
   It is a property of what was loaded. Where loading information is missing, the assumed value and
   its source are stated explicitly (see §2).

#### Per-cohort classification (2026-09-21)

| Cohort | `channel_id` is | Row | Basis / uncertainty |
|---|---|---|---|
| COMBAT | `Channel_ID` from the deposit's own obs (70 GEM wells = 10 pools × 7) | 1 — explicit capture ID | Depositor-provided capture identifier |
| Liu | 10x lane (`B1_10xlane1` …, 55 lanes) | 3 — pool | 14–19 donors multiplexed per lane; grouped by lane, **never** by post-demuxlet donor |
| Cillo | GSM = one GEX library | 3 — pool | 27/48 libraries carry 2 samples, 10 pool 2 distinct donors |
| Ren | GSM = `sampleID` | 2 — sample as trusted proxy | Channel inventory records one library per GSM, 1:1 with GEM well |
| Qiu | GSM = one donor × one timepoint | 2 — sample as trusted proxy | One library per GSM in the GEO record |
| Chongqing | GSM = one donor | 2 — sample as trusted proxy | One library per GSM |
| Stephenson | `sample_id` (143 across 130 patients) | **4 — technical grouping, approximate** | Whether one `sample_id` equals exactly one capture is **not established** from the deposit; used as the finest trustworthy grouping, with the uncertainty recorded. Scrublet output for this cohort is an approximate prediction. |
| Wilk | GSM = one Seq-Well array | 2 — sample as trusted proxy, **but platform differs** | Seq-Well, not 10x. No loading record exists for this cohort, so the 10x 0.8%/1k multiplet figure is reused only as a **stated assumption**, flagged unverified for this platform (§2). |

---

## 1. Cell-level QC

| # | Step | Status | Notes |
|---|---|---|---|
| 1 | Empty-droplet / cell calling | **DONE** (unverified) | All 8 cohorts use CellRanger's own *filtered* (called-cell) matrices. See §4 for what "DONE" must show before final acceptance. |
| 2 | Absolute floor: `min_genes` / `min_counts` per cell | **DO**, provisional | Complements #4: MAD is relative to each channel's own distribution, so a uniformly low-quality channel can pass MAD while below any sane absolute floor. Both **potentially bias against platelets** — which of them actually does, and by how much, is an empirical question the flags report answers; rev. 2's claim that `min_counts` is "the filter most likely" to do so had no data behind it and is withdrawn. |
| 3 | Mitochondrial fraction (`pct_counts_mt`) | **DO**, provisional | High mito% *can* indicate damage; it does not alone prove it. Rev. 3 said low-RNA cells "inflate" mito% — too absolute: a small denominator does not deterministically raise the ratio, it makes the ratio **unstable and high-variance**, so a few counts either way move it far more than in a high-RNA cell. **Platelets are anucleate but do carry mitochondria — the nucleated-cell interpretation does not transfer**; OSCA likewise notes high mito% can reflect genuine biological variation. 20% is a candidate to validate, not a quality truth. |
| 4 | MAD outliers on `log1p_total_counts`, `log1p_n_genes_by_counts` | **DO**, provisional | Per channel, on log1p values. Full specification in §2b — rev. 3's bare "5 MADs" left too much unstated. |
| 5 | Ribosomal fraction (`pct_counts_ribo`) | **RECORD** | Ribosomal content varies genuinely with cell type and translational state, so a hard threshold **may remove real biological variation and cannot specifically distinguish artifact** — it is not that it removes only biology, but that it does not separate the two. This argument stands alone (see §3 for the withdrawn DecontX-based half). |
| 6 | Hemoglobin fraction (`pct_counts_hb`) | **RECORD** | Blood-derived material carries free-hemoglobin background by nature; a hard cut risks deleting real Erythroid cells. Same withdrawal applies. |
| 7 | Doublets — cross-donor (multiplexed pools) | **DONE** (unverified) | Liu: demuxlet SNG-only. COMBAT / Cillo: depositors' HTO calls. Non-multiplexed cohorts have one donor per library, so cross-donor doublets cannot occur by construction. |
| 8 | Doublets — expression-profile detection (Scrublet) | **DO**, see §2 | Catches doublets whose *combined expression profile* is distinguishable — **a portion of** same-donor events, not all. **Same-donor ≠ homotypic**: same-donor T+B is comparatively detectable; homotypic (e.g. same-donor T+T) is much harder for any expression-based method. **A portion may remain undetected, and the residual proportion is unknown** — stated as a limitation, not quantified or implied away. |

## 2b. MAD specification (rev. 4)

"5 MADs" alone is underspecified. The implemented rules:

| Item | Rule |
|---|---|
| Scale | **Unscaled**: `median_abs_deviation(..., scale=1.0)`, stated explicitly rather than left to a default. `scale='normal'` divides by 0.67449, which would silently make the same nominal "5 MADs" ≈1.48× tighter. |
| Tails | Low and high tails are **flagged separately** (`flag_mad_*_low` / `flag_mad_*_high`) so the audit shows which end drives attrition. |
| Which tails delete | Only the **high** tails delete by default. The low tail is computed and reported but not used for removal: it overlaps heavily with the absolute floors (#2) and with genuinely small cells such as platelets, so letting it delete would double-penalise exactly the population at risk, and less legibly than an explicit floor. Revisit only with the audit in hand. |
| `MAD == 0` | **No cells flagged**, status `mad_zero`. The degenerate case must not delete: with MAD 0 the bounds collapse onto the median and every cell differing at all would be flagged. A near-constant metric means the statistic is uninformative here, not that the channel is one big outlier. |
| Missing values | `NaN` is **never flagged**; a missing metric is not evidence of an outlier. Counts of non-finite values are recorded per channel. |
| Small channels | Channels under 70 cells get **no MAD flags at all** (status `too_small`) — a median/MAD pair fitted on a handful of cells is not a distributional summary. Project rule, not a general boundary. |

Per-channel median, MAD, status and low/high counts are written to
`audit_mad_per_channel.csv` so every one of these decisions is inspectable rather than implicit.

## 2. Scrublet execution constraints (rev. 3, rates corrected in rev. 4)

| Constraint | Requirement |
|---|---|
| Input values | **Unnormalized raw counts.** Not log-transformed, not scaled, not DecontX-corrected. |
| Input cell population | **The capture channel's full cell population**, including its major cell types — Scrublet's simulated-doublet null is built by combining the observed cells, so running it on an already-extracted Platelet-only subset would simulate doublets from a population that does not reflect what was actually co-captured. Target-lineage extraction happens *after* Scrublet, never before. |
| Expected doublet rate | Two distinct quantities, both recorded — **rev. 4 corrects rev. 3's formula**. (a) **loading-based** rate `r` for the original capture (10x ≈ 0.8% per 1k cells recovered, scaled by channel size). (b) **post-demultiplexing residual** rate. Rev. 3 used `r / n_donors`, which is wrong: `1/n` is the probability that a doublet *happens to be same-donor*, not the doublet *fraction of what survives* — removing cross-donor doublets also shrinks the denominator. Correct conditional proportion, with `p_d` each donor's **observed** share of that channel's cells (not an assumed even split) and `s = Σ_d p_d²`: `r_residual = r·s / ((1−r) + r·s)`. Worked example (5 donors evenly pooled, r = 10%): s = 0.20, r_residual = 0.02/0.92 ≈ **2.17%**, versus 2.00% from the naive `r/n`; the gap widens with r. **Stated assumptions**: all cross-donor doublets detected and removed, all singlets retained, nothing else filtered — none exactly true (real demultiplexing misses some doublets, may drop ambiguous singlets, and later QC changes composition again), so (b) is a documented estimate used to set Scrublet's prior, not a measured quantity. Scrublet gets (b) for multiplexed channels, (a) for single-donor ones; both values, `s`, and the choice go to the per-channel diagnostics. |
| Threshold | Automatic per channel is the **starting point, not the finish**. The selected threshold, the observed score distribution and the simulated-doublet score distribution are saved per channel so the separation can be inspected before the calls are trusted. |
| Minimum channel size | 70 cells — a **project rule** for skipping threshold fitting on tiny channels, not a general reliability boundary. |
| Status preservation | `skipped` (channel too small), `failed` (Scrublet errored) and `unreliable` (threshold inspected and rejected) are **distinct recorded statuses, never collapsed into "singlet"**. A cell that was never successfully scored is not thereby a confirmed singlet, and must not be silently promoted to one downstream. |

## 3. Gene-level QC and ambient correction

| # | Step | Status | Notes |
|---|---|---|---|
| 9 | `min_cells` per gene | **NO** | Per-channel low-frequency gene removal gives each channel a different gene set, and a gene in 2 cells of a small channel is not thereby biologically void. **Retain the full input gene space**; each FM then handles its own vocabulary, gene ordering and sequence-length limits at its own input. Revisit only at a specific algorithm's input if that algorithm demonstrably requires it — never irreversibly on the shared matrix. |
| 10 | Ambient RNA correction (DecontX) | **DONE** (refit required) | **Fit on the QC-retained, FULL-LINEAGE reference population of each channel, then extract the target lineages — these are not the same cell set.** DecontX estimates contamination using the channel's cell-type composition (and can optionally take an explicit empty-droplet background), so restricting its input to Platelet+Monocyte would change what it treats as ambient. The refit exists to keep the correction consistent with the post-QC cell set, not because the earlier fit was invalid. |

### Inferences withdrawn (rev. 2, retained here)

Both read model *output composition* as evidence of model *correctness* — different questions:

1. ~~"RPL/RPS were the dominant genes DecontX removed, therefore ribosomal contamination is already
   correctly handled."~~ Removal volume is not validation.
2. ~~"Liu's Erythroid 79% contamination matches the expected hemoglobin-ambient pattern, therefore
   it is fine."~~ Consistency with one explanation is not confirmation; mis-annotated clusters and
   over-correction were never ruled out.

Neither withdrawal changes the **RECORD** decisions for ribo/hb — those rest on the
biological-variation argument, which is independent.

## 4. Acceptance criteria for a "DONE" claim (rev. 3)

"DONE (unverified)" is an honest status, but it is a *claim of processing*, not traceable evidence.
Before final acceptance each DONE row must show:

- the **input source** (file path / accession / deposit version),
- the **filter fields** actually applied (column names and values, not a prose description),
- **before/after cell counts**, per channel,
- **barcode-level correspondence** between input and output (so the retained set is reconstructible,
  not merely counted).

## 5. Corrected pipeline order

```
CellRanger filtered matrices (DONE)
  → cross-donor demux / HTO filtering (DONE, multiplexed cohorts only)
  → [NEW] per-cell QC FLAGS computed, nothing deleted  ← flags-first gate (§0.1)
  → [NEW] audit: per-filter + cumulative, by channel / donor / lineage (§0.2)
  → [GATE] rules reviewed and FROZEN, blind to outcome (§0.3)
  → [NEW] apply frozen filters; Scrublet per channel on full-population raw counts (§2)
  → DecontX refit per channel on QC-retained FULL-lineage set (§3 #10)
  → extract target lineages (Platelet + Monocyte_Myeloid)
  → platelet + monocyte assembly (RERUN)
  → 5-FM embedding extraction, full scale (RERUN)
```

## 6. Thresholds — all provisional, frozen only after the §0.2 audit

| Parameter | Candidate | Basis |
|---|---|---|
| `min_genes` per cell | 200 | Convention for 10x; re-examine for anucleate low-RNA platelets. |
| `min_counts` per cell | 500 | Same caveat; potentially biases against platelets — magnitude to be measured, not assumed. |
| `pct_counts_mt` cap | 20% | Candidate only; nucleated-cell interpretation does not transfer to platelets. |
| MAD threshold | 5 MADs | Convention, on log1p values, per channel. |
| Per-donor platelet minimum | to be prespecified | Needed to define "donor lost to QC" in §0.2; must be fixed before the audit, not chosen after seeing it. |
| Scrublet | see §2 | Rate (a) and (b) both recorded; threshold inspected, not assumed. |

## 6b. Gate-1 trial findings (Qiu + Cillo, 2026-09-21)

Basic-QC flags only — **not** final Phase 1 retention; Scrublet removes a second tranche and this
table must be regenerated after Gate 2 (§0.1).

**Donor impact, both denominators required** (Cillo):

| Measure | Value |
|---|---|
| donors in cohort | 29 (only 27 have any platelets; COVID_17 / COVID_23 have none at all) |
| donors ≥20 platelets **before** QC | **11** |
| donors ≥20 platelets **after** QC | 8 |
| QC-newly-failing / all donors | 3/29 = **10.3%** |
| QC-newly-failing / previously qualifying | 3/11 = **27.3%** |

The three: COVID_12 (22→17), COVID_2 (21→17), COVID_3 (24→18). "No donor zeroed" is not the same
as "no patient leaves the analysis" — at the pre-registered ≥20 threshold these three exit.

**Attribution — the mito cap does essentially all of it.** Platelet flag rates: `flag_mito` 14.08%
vs 0.78% in monocytes (18×); `min_counts` 1.68%, `min_genes` 0.60%. In the three affected donors
specifically, 14 of 15 removed platelets were mito-flagged and **`min_genes`/`min_counts` flagged
zero**. Rev. 2's guess that the count floor would be the platelet-hostile filter was wrong in both
directions — this is the empirical answer the flags report exists to give.

**Why**: the same absolute 20% cap sits at very different places in the two distributions —
platelet `pct_counts_mt` p50 8.2 / p90 23.8 / p95 31.3, monocyte p50 4.9 / p90 9.0 / p99 18.7. One
global cap therefore cuts ≈p90 of platelets but beyond p99 of monocytes: the stringency is
lineage-dependent by construction, which is exactly the "nucleated-cell interpretation does not
transfer" concern in §1 #3, now observed rather than predicted.

**Ambiguity, stated rather than resolved**: high-mito platelets do have lower complexity (median
1,493 counts / 608 genes vs 2,545 / 874), so a real quality gradient exists — but **0 of 117 were
flagged by any count/gene/MAD filter**, so no independent filter corroborates that they are
degraded rather than simply platelet-baseline. Per §0.4 this triggers review; it does not by itself
license moving the number.

## 6c. GATE 1 — FROZEN 2026-09-21

**Decision: Scheme A. Global `pct_counts_mt` < 20%, unchanged.** Frozen by the reviewer after the
targeted platelet review and the A-vs-B comparison below. Basic QC is now fixed; execution
proceeds to Scrublet (Gate 2).

**Frozen rule set**

| Rule | Value | Note |
|---|---|---|
| `pct_counts_mt` | < 20%, global (all lineages) | Not a claim that 20% is optimal — see rationale |
| `min_genes` | ≥ 200 | |
| `min_counts` | ≥ 500 | |
| MAD | 5, unscaled (`scale=1.0`), per channel, **high tail deletes only** | §2b |
| `min_cells` per gene | **not applied** | §3 #9, full gene space retained |
| `pct_counts_ribo`, `pct_counts_hb` | recorded, never filtered | §1 #5/#6 |
| per-donor platelet minimum | 20 | Pre-registered 2026-09-15 (benchmark_design_decisions.md Decision 3), with 10 / 50 sensitivity also pre-registered |

**Rationale — this is a decision under acknowledged uncertainty, not a finding that A is correct.**
The evidence supports only this: the high-mito platelets are consistent with platelet identity
(median 6 of 13 canonical markers detected, 66.7% with ≥5, versus 3 / 20.2% in monocytes; no
obvious heterotypic mixed-expression signal) **and** have systematically lower library complexity
(median 1,493 counts / 608 genes / 8.2% of counts in the top gene, versus 2,545 / 874 / 3.6%). It
does **not** establish that mito > 20% must be removed, nor that retaining all of them is without
effect.

**Why B was rejected** (per-lineage MAD on `pct_counts_mt`):

1. **It does not only adjust platelets.** It would move every lineage's threshold — monocyte to
   12.3% (flagging ×5: 0.78% → 4.11%), T cell to 10.0% (×7.5), NK to 12.1% (×5.7). Increased
   filtering is not itself evidence of error, but a platelet-focused review cannot justify those
   other changes.
2. **The computed B is a cross-channel pooled estimate**, not the within-channel per-lineage rule
   it was described as. Platelets average ~24 cells per channel here (831 across 35), far below
   the 70-cell floor the MAD spec requires, so per-channel fitting is infeasible — the pooled
   version mixes in between-channel variation and is a different scheme, not the same one
   approved.
3. **Saving donors is not a valid basis for choosing it.** B retains COVID_2 and COVID_3 (and not
   COVID_12), but COVID_2 lands on exactly 20 — donor qualification here turns on single cells.
   That sensitivity is reported through the pre-registered 10 / 20 / 50 thresholds, not engineered
   away by threshold choice.

**Accepted, documented cost of A**: in Cillo, 3 of 11 previously-qualifying donors (COVID_2,
COVID_3, COVID_12) fall below the 20-platelet minimum; platelet mito flag rate is 14.08% versus
0.78% in monocytes; 692 of 831 platelets retained. Qiu: no donor-level loss.

**Correction recorded**: an earlier note in this review described median + 5·MAD as landing at
"≈p95" for every lineage. That is wrong — where median + k·MAD falls depends on each
distribution's shape, and it does not delete an equal proportion across lineages. The B column
itself shows the spread: 8.18% of platelets versus 3.76% of T cells.

**Scope**: these are basic-QC numbers only. Scrublet removes a further tranche, so the donor ×
lineage retention table is regenerated after Gate 2 before any of it is treated as final Phase 1
retention.

## 6d. Input-state diagnostics — Q1–Q4 (2026-09-21)

**Method principle (governs every threshold below).** *We harmonize QC objectives and audit fields
across cohorts, while calibrating operational thresholds to platform and upstream input state. The
aim is comparable retained-cell quality, not comparable deletion rates.*

### Q1 — is Gate 1 computed on raw integer counts? **PASS, all 8 cohorts**

| Check | Result |
|---|---|
| X is non-negative integers | Pass — all 8 `all_integer=True`, min value 1. Liu/Stephenson/COMBAT are float64 *dtype* but hold integer values (MatrixMarket "real" storage), not normalized data |
| Not `normalize_total` / `log1p` / `scale` / DecontX-corrected | Pass — non-integer values would appear otherwise; input read from `*_decontx_input/counts.mtx`, the PRE-DecontX raw input (DecontX output is a separate `*_decontx_counts.mtx`) |
| QC metrics recomputed from that matrix, not read from archived fields | Pass — `qc_gate1_apply.py:compute_flags()` derives total/n_genes/pct_mt from the matrix directly |
| Q2 percentiles from the pre-filter full set, not the retained object | Pass — `qc_flags_per_cell.csv` holds ALL cells with flags; deletion is applied separately to `cells.csv`. **Caveat**: Ren/COMBAT figures come from their first 40 channels (full cells within those channels), not all 284/70 |

No cohort's Gate 1 result is void on these grounds.

### Q2 — which cohorts show hard truncation from upstream filtering?

**Clear evidence: COMBAT, Ren, Stephenson.** Possible count floor: Qiu, Chongqing.

| Cohort | min_counts | min_genes | max_mito | Reading |
|---|---|---|---|---|
| COMBAT | 336 | **300** | **10.0** | mito capped at exactly 10.0%, genes floor at exactly 300 |
| Ren | 982 | 378 | **10.6** | high floors on both, mito capped ≈10% |
| Stephenson | 305 | **200** | 11.4 | genes floor at exactly 200 |
| Qiu | **500** | 14 | 99.8 | counts floor at exactly 500; no mito cap |
| Chongqing | **500** | 9 | 99.8 | same |
| Cillo / Liu / Wilk | 203 / 100 / 48 | 174 / 11 / 43 | 68.2 / 92.0 / 82.3 | no comparable hard boundary |

These discrete boundaries **strongly support upstream threshold filtering**, but final confirmation
still rests on the raw-matrix verification above plus the deposits' own processing records.
**For already-truncated cohorts, this project's QC can only state that no ADDITIONAL cells violate
our rules — not that those cohorts share the same original quality distribution as the others.**

### Q3 — is Wilk's 57.5% a low-quality tail? **No — a platform-level distribution shift.**

Wilk counts percentiles: p10 126, p25 195, **p50 433**, p75 1,074, p90 2,504. Genes: p10 96,
p25 126, **p50 244**, p75 535. The median already sits below the 500-count candidate floor, and
54.4% / 43.1% of cells fall below the count / gene floors — against 10x cohort count medians of
≈3,000–4,500. Applying counts ≥500 / genes ≥200 here is not trimming an aberrant tail; it is
cutting the body of a Seq-Well distribution with a 10x ruler — **platform-mismatched filtering that
manufactures artificial cohort shift**.

**Result (implemented 2026-09-21)**: *Platform-adaptive filtering eliminated the severe
lineage-dependent attrition produced by applying 10x-derived absolute thresholds to Seq-Well data,
while retaining removal of within-platform extreme outliers.* Cohort removal 4.0% (was 57.5%);
Platelet 2.46% (was 76.7%); per-lineage removal now spans 2.5–6.2%. Thresholds came from
capture × lineage for 245,143 cells and a borrowed cohort-level same-lineage threshold for 864
cells in strata too small to fit their own MAD. Donor-level platelet retention: ≥20 platelets
22 → 21 donors; ≥50 platelets 20 → 20. **This shows the artificial platform bias was corrected; it
does not by itself establish comparable biological quality** — that rests on the subsequent marker,
contamination and donor audits.

**Decision for Wilk**: do not apply the 500/200 10x floors. Keep the shared QC objectives and audit
fields; calibrate with Seq-Well-appropriate rules (within-platform distribution / MAD, an explicit
very-low-complexity floor, marker-consistency checks). 500/200 are retained as **diagnostic columns
only**. Thresholds frozen before any mortality outcome or model performance is seen.

### Q4 — what are Liu's flagged barcodes? **Background-region barcodes, not low-RNA cells.**

| | n | counts median (p25–p75) | genes median | mito median |
|---|---|---|---|---|
| Removed | 161,527 | **234** (180–315) | 165 | 1.6% |
| Retained | 629,740 | **4,701** (3,263–6,236) | 1,543 | 3.0% |

**Separation**: removed p90 = 530 counts versus retained p10 = 1,694 — the two populations barely
overlap, a ~3× gap. That is a bimodal split at the barcode-rank knee, not a continuous low-RNA cell
population. The removed set's *lower* mito (1.6% vs 3.0%) points the same way: ambient/background
barcodes are dominated by cytoplasmic mRNA rather than mitochondrial, whereas damaged cells run
high.

This is consistent with the documented fact that Liu's per-lane "filtered" matrices are a fixed
**30,000-barcode demuxlet candidate whitelist**, not an independent cell call. Per the decision
rule: **Liu needs proper cell calling first, then cell-level QC on the called cells** — `min_counts`
/ `min_genes` must not stand in for cell calling. Reporting stays layered:
`candidate barcodes → demuxlet SNG (791,267) → QC-retained (629,740, 79.6%)`, so "not a credible
singlet" stays distinguishable from "credible singlet, low complexity".

**Not yet directly tested** for Liu: coordinated lineage-marker coherence in the flagged set, and
explicit barcode-rank knee placement. The count/gene/mito separation is strong but these two checks
were listed and have not been run.

Separately, 14,173 removed Liu cells were flagged by mito ALONE (otherwise passing) — a distinct
group of genuine cells with high mito, not background barcodes.

## 6e. FROZEN 2026-09-21 — Wilk and Liu handling

*We harmonize QC objectives and audit fields across cohorts, while calibrating operational
thresholds to platform and upstream input state. The aim is comparable retained-cell quality, not
comparable deletion rates.*

| Cohort | Decision | Core reason |
|---|---|---|
| Wilk | Drop the 10x-style `min_counts=500` / `min_genes=200`; use Seq-Well-internal adaptive QC | Whole distribution shifted left; not an aberrant tail |
| Liu | Re-do cell calling first, then cell-level QC | Current input carries many background barcodes; QC cannot substitute for cell calling |

### Wilk — Seq-Well-adapted rules (frozen)

1. No `counts ≥ 500`, no `genes ≥ 200`. Both retained as **diagnostic columns only**.
2. `log1p_total_counts` and `log1p_n_genes` get a **lower-tail 5-MAD cut only** — not two-sided
   deletion.
3. Thresholds estimated per **Seq-Well × broad lineage**, so naturally low-RNA platelets are not
   declared outliers by other cell types' distributions.
4. Estimated **within capture/sample** where a reliable ID exists and that stratum has enough
   cells; small strata borrow Wilk's cohort-level same-lineage threshold — a tiny sample never
   fits its own MAD.
5. **mito: no global 20%.** Seq-Well × lineage **upper-tail 5 MAD**, with per-cell flags kept for
   audit.
6. Donor-level platelet retention keeps being reported, especially donors near 20 / 50 — but
   thresholds are never reverse-tuned on whether a donor is preserved.

5 MAD is chosen to remove only the extreme tail, not to redefine Seq-Well cells as needing to reach
10x complexity. Batches with systematically different sequencing depth warrant separately estimated
adaptive thresholds, while watching that cell types with genuinely different RNA content are not
deleted for it.

### Liu — cell calling, FALLBACK path (confirmed by deposit contents)

**Verified 2026-09-21**: GSE161918 deposits 114 files — 56 RNA + 56 CSP
`*_filtered_feature_bc_matrix.h5` — and **zero raw/unfiltered droplet matrices**. The preferred
EmptyDrops-on-raw-droplets path is therefore unavailable: the lower ambient tail was already
truncated upstream, so EmptyDrops here could not be described as standard raw-droplet cell calling.

Fallback, as frozen:
- per-lane **barcode-rank curve**, with **knee** and **inflection** computed and saved;
- a **data-derived per-lane boundary** from the empirical valley between the two count components
  (the bimodal split already evident: removed p90 = 530 vs retained p10 = 1,694);
- **EmptyDrops only as a sensitivity analysis**, never the primary caller;
- **forbidden**: repackaging fixed 500-count / 200-gene floors as cell calling.

Pipeline of record for Liu:

```
fixed 30,000 candidate barcodes per lane
  → barcode-rank data-derived cell calling  (knee / inflection / valley saved per lane)
  → intersect with demuxlet SNG   (noting demuxlet was originally run on the candidate set)
  → cell-level QC
  → Scrublet
  → DecontX
```

### Liu cell-calling rule — FROZEN (inflection)

**`called cell = total_counts ≥ that lane's own inflection`**, computed independently per lane.

| Metric | Role |
|---|---|
| **Inflection** | **primary cell-calling threshold** |
| Knee | diagnostic only — check it roughly agrees with inflection |
| Valley | diagnostic only — never deletes; peak detection can jump to the wrong mode |

- The threshold is each lane's **own** inflection, never the global 529. The 319–1,043 spread
  across lanes reflects real sequencing-depth and experimental differences and is not flattened.
- Inflection comes directly from the barcode-rank curve's geometry, making it more stable than a
  valley that depends on histogram binning and bimodal peak location, and than a knee that already
  produced a spurious extremum (B1_10xlane4: knee 113 at rank 13,667/13,708).

**Automatic flagging, never automatic threshold replacement.** Saved per lane: inflection, knee,
valley, their ratios, barcode counts either side of the threshold, the barcode-rank curve,
`pct_above_inflection`, and the count/QC profile of barcodes within ±20% of the threshold. A lane
is flagged for human inspection if its inflection is itself an extreme outlier across the 55 lanes,
or if the curve shows no clear turn — but disagreement from knee or valley never auto-substitutes
for inflection. The four known failures are all *auxiliary-detector* failures; inflection stayed in
range in each, so no lane threshold is changed on their account.

**Naming in the methods**: *"Barcode-rank–based cell calling using a lane-specific inflection
threshold."* **Not** EmptyDrops — this is not an ambient-profile hypothesis test on a complete
raw-droplet matrix.

**On 500 vs 529**: *the median lane-specific inflection was approximately 529 UMIs; its proximity
to the conventional 500-UMI threshold was coincidental and served only as external concordance.
Cell calls were determined independently from each lane's barcode-rank curve.* The two may discard
a similar number of barcodes, but they mean different things: 500 treats every input barcode as an
already-called cell and then filters on quality; inflection first decides which barcodes look like
cell-containing droplets.

### Why Liu's low Erythroid / Unclear called-rate is the expected pattern

**Pre-call lineage label ≠ confirmed cell identity.** These labels were assigned to *candidate
barcodes* before cell calling, from profiles that include ambient RNA. Free HBA/HBB is highly
abundant in blood, so empty droplets are readily labelled Erythroid, and background barcodes
lacking clear lineage signal fall into Unclear. Those two categories are therefore exactly where
background barcodes concentrate — **Erythroid called 41.5% does not mean 58.5% of real erythroid
cells were deleted**, and Unclear's 37.3% fits the same explanation most directly.

The rise in mito to 15.2% among barcodes within ±20% of the threshold shows the transition zone is
enriched for damaged or marginal-quality barcodes, supporting the overall placement of inflection —
but it shows only that quality degrades there, and cannot by itself establish that every uncalled
barcode is an empty droplet.

A **light Erythroid diagnostic** is recorded (called vs uncalled counts/genes/mito, erythroid
marker breadth, HBA/HBB dominance, DecontX contamination) purely to document why the pre-call label
cannot be read as true cell identity. It is **not** used to rescue barcodes or to re-optimise any
threshold, and it does not block Gate 2.

### Liu — ACCEPTED 2026-09-21, Gate 2 open

**Methods statement of record:** *Lane-specific barcode-rank cell calling retained 644,679 of
791,267 candidate barcodes (81.5%). No donor was completely depleted of platelet cells; however,
seven donors no longer met the prespecified minimum of 20 called platelet cells, reducing the
eligible platelet donor set from 49 to 42. This attrition was distributed across donors rather than
driven by a small number of anomalous samples.*

**Framing that matters**: not "cell calling caused the loss of seven donors", but *after candidate
barcodes were restricted to confidently called cells, seven donors no longer met the prespecified
minimum*. Cell calling **revealed** that those donors' effective platelet counts were insufficient;
suspected background barcodes are not re-counted as cells in order to preserve a donor total.

Why the seven are not grounds to relax inflection:
- the primary detector was confirmed working in the one outlier lane (B2_10xlane06: slope −9.79 →
  **−79.45** → −60.09, a sharp local minimum; knee 1,048 vs inflection 1,043, a 0.5% difference);
- Platelet called-rate 72.7% sits alongside Monocyte 69.1% and B_cell 73.0% — no platelet-specific
  suppression;
- no donor was zeroed, and loss is distributed (top-5 donors hold 36% of total platelet loss);
- pre-calling "platelet counts" came from **pre-call labels on candidate barcodes**, which can
  include background barcodes labelled Platelet by ambient PF4/PPBP.

**Frozen:**
- lane-specific inflection unchanged;
- the 20-cell threshold remains the prespecified primary inclusion criterion;
- the 50-cell threshold continues as the stricter sensitivity analysis;
- the 7 donors dropping below 20 are excluded from the **platelet-restricted** patient-level
  analysis, but stay in the **Monocyte** analysis where their monocyte counts qualify;
- complete pre/post donor count tables and the 7-donor list are retained
  (`liu_donor_Platelet_callcounts.csv`, `liu_donor_Monocyte_Myeloid_callcounts.csv`);
- **no donor is restored and no threshold is adjusted on the basis of mortality outcome.**

Donors below 20 after calling: AA220014, HGR0000074, HGR0000083, HGR0000089, HGR0000091,
HGR0000136, IM000040. Below 50 (sensitivity arm): CHI014, HGR0000078, HGR0000093, HGR0000098,
HGR0000141, HGR0000143.

### Two outstanding checks

- **Barcode-rank knee plots: required.** Both a cell-calling diagnostic and part of the final
  methods record. In progress.
- **Lineage-marker coherence: recommended, not blocking.** The ~3× non-overlapping gap between
  removed and retained already establishes that a separate low-count component is present; the
  marker analysis serves to confirm whether genuine low-RNA cells sit inside that component.

### The 14,173 mito-only Liu cells

Not to be merged with background barcodes. They pass through **cell calling first**, then are
judged independently by the new lineage-aware mito rule.

### Ren / COMBAT truncation — now a whole-cohort finding (all channels, 2026-09-21)

*Across all 284 Ren and 70 COMBAT capture channels, no cell exceeded 12% mitochondrial counts.
COMBAT showed an exact 10.00% upper boundary in every channel, while Ren showed a median channel
maximum of 9.98%. Together with the corresponding lower boundaries in detected genes and total
counts, these distributions strongly support upstream threshold-based QC rather than naturally
occurring truncation.*

Supporting figures: Ren — 1,462,702 cells, cohort min_counts 878, min_genes 378, max_mito 10.60,
per-channel max_mito median 9.98 / p95 10.03, **0/284 channels above 12%**. COMBAT — 836,148 cells,
min_counts 336, min_genes 300, max_mito 10.00, per-channel max_mito median = p95 = max = 10.00,
per-channel min_genes median 301 / max 306, **0/70 channels above 12%**. No longer a sampled
inference.

## 7. Plan of record

Per the reviewer's judgment: no rewrite needed. Add the donor-level audit, the flags-first gate and
the Scrublet/MAD constraints (this revision), then **a small-scale trial run** on a limited set of
cohorts. Only after the trial passes review does the full rebuild start.
