# GSM5102900_HC1 scoring pilot — result interpretation (v1)

**2026-10-12. Read-only. No job submitted, no cell deleted, no outcome field read,
`PRODUCTION_PARAMETER_RULING` untouched.**

Responds to `review:20261012_001137` §2–§3 and its next-step instruction item 1.
Interprets the single accepted technical run — SLURM job `29882644`, output directory
`GSM5102900_HC1__cfg0f02df7c7e2a__impl30831273d4d3` (see `implementation-note:20261012_000549`
for the run's own technical acceptance). This document adds nothing to that acceptance; it only
explains what the eight frozen priors actually produced.

**Revised in place (2026-10-12) per `review:20261012_124351` (F1, F2, F3).** The first
version of this document (a) wrongly described the automatic threshold as rising monotonically
with the prior "as the algorithm requires" and asserted an unverified "flat-region crossing"
mechanism for the 0.16→0.20 drop, (b) miscounted 4 distinct candidate sets as 3, and (c) claimed
the diagnostic plots were reproducible from tracked files when the plotting script actually read a
live, uncommitted HPCC path. All three are fixed below; the underlying data-level facts (nesting,
argsort identity, moved-in/out counts, Jaccard matrix) were correct and are unchanged.

## 1. What this corrects relative to the submission handoff

`implementation-note:20261012_000549` described `n_called` as "rising with the prior as
expected for an algorithmic stress grid." That is not what the data shows. Recomputed directly
from the raw sidecars (not copied from the manifest's `grid_stability` block or from the review),
using
[`gate2_pilot_interpretation_evidence_20261012/gate2_pilot_interpretation_analysis.py`](gate2_pilot_interpretation_evidence_20261012/gate2_pilot_interpretation_analysis.py),
whose full output is
[`gate2_pilot_interpretation_evidence_20261012/gate2_pilot_interpretation_recomputed.json`](gate2_pilot_interpretation_evidence_20261012/gate2_pilot_interpretation_recomputed.json):

- `n_called` is **not monotonic** in the prior: it falls from 291 (prior 0.16) to 251
  (prior 0.20).
- The eight priors collapse into **four distinct candidate sets** (machine-counted by value —
  `n_distinct_candidate_sets` in the recomputation JSON — not read off a table by eye, which is
  exactly how a prior version of this document miscounted 4 as 3): `{0.005, 0.01, 0.02}` → the
  identical 8-cell set; `{0.04, 0.08}` → the identical 225-cell set; `{0.12, 0.2}` → the identical
  251-cell set; `{0.16}` alone → its own distinct 291-cell set, which is itself the 4th group, not
  an unlabeled remainder. ("Identical" within a group is proved, not inferred from a rounded
  Jaccard: for each pair in a group the recomputation independently confirmed `nested=True` *and*
  equal cardinality, which together force set equality — see §3.)
- The candidate *count* changes because each prior's automatic threshold is fit **independently**
  (§3) — Scrublet does not fit one joint threshold across priors, and nothing in the algorithm
  requires the resulting cut to retain a non-increasing or non-decreasing share of cells as the
  prior changes.

## 2. Per-prior thresholds, candidate counts, and proportions

Recomputed from `PER_PRIOR_CALLS.json`/`PER_PRIOR_SCORES.json` (`n_called` independently summed
from the raw boolean vectors, not read from the manifest) and cross-checked cell-by-cell against
`SCORING_MANIFEST.json`'s per-prior `threshold`: for every one of the 7,666 cells and all 8
priors, `call == (score > threshold)` holds exactly — the mechanical scoring rule was applied
consistently everywhere, confirmed below the counting level, not just at it.

| prior | auto-threshold | candidates | % of 7,666 cells | identical candidate set to |
|---:|---:|---:|---:|---|
| 0.005 | 0.088932 | 8 | 0.1044% | 0.01, 0.02 |
| 0.010 | 0.163174 | 8 | 0.1044% | 0.005, 0.02 |
| 0.020 | 0.283172 | 8 | 0.1044% | 0.005, 0.01 |
| 0.040 | 0.188922 | 225 | 2.9350% | 0.08 |
| 0.080 | 0.311002 | 225 | 2.9350% | 0.04 |
| 0.120 | 0.409099 | 251 | 3.2742% | 0.20 |
| 0.160 | 0.448088 | 291 | 3.7960% | (none — distinct) |
| 0.200 | 0.540506 | 251 | 3.2742% | 0.12 |

**The threshold is not monotonic either**, confirmed by machine check
(`threshold_monotonic_nondecreasing_across_prior_grid_order: false` in the recomputation JSON):
it falls from 0.283172 (prior 0.02) to 0.188922 (prior 0.04). A prior version of this document
claimed the threshold "rises monotonically... as the algorithm requires" — false on both counts:
the recorded values themselves contradict it, and (§3) nothing in Scrublet's automatic-threshold
procedure requires monotonicity across independently-fit priors in the first place. Candidate
count is flat within each of the four groups in §1, then **drops** 291 → 251 going from prior 0.16
to 0.20. 40 cells leave the candidate set between prior 0.16 and prior 0.20; zero cells enter
(`prior_0.16_to_0.2_moved_out_count: 40`, `prior_0.16_to_0.2_moved_in_count: 0`, recomputed).

## 3. Why this happens: the real algorithm, and what is and is not verified

**The algorithm, read from the installed Scanpy 1.11.3 source itself** (`inspect.getsource` on the
live HPCC environment's `scanpy.preprocessing._scrublet.core.Scrublet`, not assumed or quoted from
memory):

- `call_doublets` sets `threshold = skimage.filters.threshold_minimum(self.doublet_scores_sim_)`
  — a histogram-valley finder applied to the **simulated** doublet scores **only** — then sets
  `predicted_doublets_ = self.doublet_scores_obs_ > threshold`. This is **not** a joint fit over an
  "observed/simulated mixture" (a prior version of this document described it that way; that
  description is wrong). Each prior gets its own, independently-fit threshold, with no constraint
  tying one prior's threshold to another's.
- `calculate_doublet_scores` passes `exp_doub_rate=self.expected_doublet_rate` (the prior) into
  `_nearest_neighbor_classifier`, which computes `doublet_scores_obs_`/`doublet_scores_sim_`
  themselves. **The prior is an input to the per-cell score formula, not only to the threshold.**
  Consequently, the eight priors' scores are not eight thresholds cut into one shared numeric
  scale — each prior has its own score distribution and its own threshold search, visible directly
  in the manifest (e.g. prior 0.005's `observed_score_max=0.142` vs. prior 0.02's `=0.40`). A
  prior version of this document treated 0.4481 (prior 0.16) and 0.5405 (prior 0.20) as two cut
  points on one fixed axis and inferred they "cross into an adjacent flat region" of a shared
  histogram; that specific mechanism was **not verified** (no step-by-step trace of
  `threshold_minimum`'s valley-finding was done) and conflated two different priors' separate
  score scales. It is retracted here, not merely softened.

**What IS verified, independently, from the raw per-cell data** (this part is unaffected by the
retraction above):

- **`all_argsorts_identical: true`.** The full rank ordering of all 7,666 cells **within each
  prior's own score vector** is byte-identical across all 8 priors (confirmed with a stable,
  index-tie-broken sort). This is a fact about relative rank order only — it does **not** imply
  the eight priors' score *values* share a common scale (they do not, per the algorithm above).
- **`all_pairs_nested: true`** for all 28 pairs of candidate sets, and **4, not 8, distinct
  candidate sets** (§1). Given stable rank order, every prior's candidate set is some top segment
  of that one fixed ranking, so a changing prior can only move the cut up or down that same
  ordered list — consistent with nesting — but **why** a given prior's independently-fit threshold
  lands where it does, in terms of `threshold_minimum`'s bimodal valley search on that prior's own
  simulated-score histogram, has not been traced step by step and is labeled **UNVERIFIED** rather
  than given a specific causal story. What can be said without further verification: since each
  prior's threshold is fit independently and the score formula itself changes with the prior,
  there is no structural reason for the resulting cut (expressed as a count or fraction of the
  fixed ranking) to move monotonically with the prior, and the data shows it does not.

**Global intersection across all 8 priors: 8 cells** (zero-indexed positions `724, 734, 1281,
1326, 2204, 2810, 3878, 7659` in the shared cell order — the same 8 cells called at every prior,
since the 8-cell sets at 0.005/0.01/0.02 are literally a subset of every larger set by nesting).
**Global union across all 8 priors: 291 cells** (exactly the prior-0.16 set, since it is the
largest and every other set nests inside it).

## 4. Full pairwise Jaccard matrix (recomputed, all 28 pairs)

| | 0.005 | 0.01 | 0.02 | 0.04 | 0.08 | 0.12 | 0.16 | 0.2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **0.005** | 1.0000 | 1.0000 | 1.0000 | 0.0356 | 0.0356 | 0.0319 | 0.0275 | 0.0319 |
| **0.01**  | 1.0000 | 1.0000 | 1.0000 | 0.0356 | 0.0356 | 0.0319 | 0.0275 | 0.0319 |
| **0.02**  | 1.0000 | 1.0000 | 1.0000 | 0.0356 | 0.0356 | 0.0319 | 0.0275 | 0.0319 |
| **0.04**  | 0.0356 | 0.0356 | 0.0356 | 1.0000 | 1.0000 | 0.8964 | 0.7732 | 0.8964 |
| **0.08**  | 0.0356 | 0.0356 | 0.0356 | 1.0000 | 1.0000 | 0.8964 | 0.7732 | 0.8964 |
| **0.12**  | 0.0319 | 0.0319 | 0.0319 | 0.8964 | 0.8964 | 1.0000 | 0.8625 | 1.0000 |
| **0.16**  | 0.0275 | 0.0275 | 0.0275 | 0.7732 | 0.7732 | 0.8625 | 1.0000 | 0.8625 |
| **0.2**   | 0.0319 | 0.0319 | 0.0319 | 0.8964 | 0.8964 | 1.0000 | 0.8625 | 1.0000 |

Minimum: **0.027491408934707903** (any of the three `{0.005,0.01,0.02}` priors vs. prior 0.16),
matching `SCORING_MANIFEST.json`'s `grid_stability.call_set_jaccard_min` exactly and matching the
review's figure. Maximum off-diagonal: 1.0, attained by the three within-block pairs noted above
plus `0.12|0.2`.

As the review states, the minimum Jaccard (`8/291 ≈ 0.0275`) is the ratio of the smallest nested
set to the largest nested set, not a comparison of two unrelated cell populations. The correct
reading is:

> In this capture, with this seed and this prior grid, the per-cell doublet score ranking is
> completely stable. What is unstable is **how much of that fixed ranking a given prior's
> automatically-fit threshold happens to retain** — and that quantity does not move monotonically
> with the prior, because the threshold can land in different flat regions of the score
> distribution.

This is **not** evidence that two mostly-disjoint candidate populations exist, is **not** a
statement about which 8, 225, 251, or 291 cells are "real" doublets, and does **not** by itself
support intersection-deletion, union-deletion, fixed-top-N selection, or picking whichever prior
looks best. `seed_stability` remains `not_computed_pending_secondary_seed_design`; nothing in this
document changes that.

## 5. Diagnostic plots

`gate2_pilot_interpretation_evidence_20261012/gate2_pilot_per_prior_histograms_overview.png`
(all 8 priors) and one PNG per prior in the same directory. **PNGs are gitignored project-wide**
(`.gitignore`: `*.png`, "optionally ignore all figures folders" — the same convention every other
stage's figures follow) and are therefore not committed.

**Reproducibility, corrected.** A prior version of this section claimed the plots were
"regenerated deterministically from the tracked, committed JSON manifest" — false: the plotting
script actually read a hardcoded *live, uncommitted HPCC path*
(`.../GSM5102900_HC1__cfg.../SCORING_MANIFEST.json`), and the committed
`gate2_pilot_interpretation_recomputed.json` did not even carry histogram data. Fixed: the
analysis script (`gate2_pilot_interpretation_analysis.py`) now also writes
`gate2_pilot_histogram_snapshot.json` — an immutable, **tracked** copy of exactly the per-prior
`threshold`/histogram fields needed for plotting, bound to its source manifest by sha256
(`source_manifest_sha256`). The plotting script now takes an explicit `--snapshot` argument and
reads only that tracked file; it verifies (not assumes) that each row's observed and simulated
histograms share identical bin edges before plotting them on one shared x-axis, and refuses to run
if they do not. Script:
[`gate2_pilot_interpretation_evidence_20261012/gate2_pilot_interpretation_plots.py`](gate2_pilot_interpretation_evidence_20261012/gate2_pilot_interpretation_plots.py),
snapshot:
[`gate2_pilot_interpretation_evidence_20261012/gate2_pilot_histogram_snapshot.json`](gate2_pilot_interpretation_evidence_20261012/gate2_pilot_histogram_snapshot.json).
No histogram or threshold was refit or recomputed for plotting in either version.

**Reading the plots.** Observed (blue) and simulated (orange) use independent log-scaled y-axes —
bar heights are not comparable across colors by eye, only each color's own trend across the
x-axis. Each prior's panel has **its own x-axis scale** (§3: the prior is an input to the score
formula itself), so a threshold line's x-position is not comparable across panels — the retracted
claim in the prior version of this section (that the 0.16 and 0.20 panels show the threshold
"crossing an adjacent flat region" of one shared histogram) is exactly the kind of cross-panel
comparison this caveat rules out, which is why both plots now print this caveat directly in their
own captions rather than leaving it to the reader to infer.

## 6. What this document does not establish

- Not a ranking of which prior is "correct"; the grid remains an algorithmic stress axis
  (`prior_label: "assumed_algorithmic_stress_setting"` in the manifest), not a prevalence
  estimate.
- Not a deletion rule. The 8-cell global intersection and the 291-cell global union are reported
  as diagnostics only, not as candidate removal lists.
- Not a statement about doublet calling on any other capture, cohort, or seed.
- Not a basis for choosing DecontX units, a batch-correction method, or a final retained set.

## 7. Scope and verification boundary

Read `SCORING_MANIFEST.json`, `PER_PRIOR_CALLS.json`, `PER_PRIOR_SCORES.json` from the existing
accepted run output directory only, plus (for §3's algorithm description) `inspect.getsource` on
the live, installed `scanpy==1.11.3` in the project's own conda environment — not fetched from an
external source, not quoted from memory. No outcome column read, no rescoring, no refit, no change
to the four original output files, no new SLURM job, no new permit. All three scripts in
`gate2_pilot_interpretation_evidence_20261012/` ran read-only on HPCC (`conda activate
302006_env_clean`), not on the login node for anything beyond trivial JSON/plot I/O, and not
locally. `gate2_pilot_interpretation_analysis.py` now takes an explicit `--manifest-dir` (no
hardcoded path) and both it and the plotting script refuse to write into an existing non-empty
output directory, so a rerun cannot silently overwrite a prior result.
