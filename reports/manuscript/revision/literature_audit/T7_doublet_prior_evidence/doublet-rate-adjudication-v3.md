# Adjudication v3 (final): stage alignment for `expected_doublet_rate` across eight cohorts

Supersedes v2, which superseded v1; neither earlier version was committed to this
repository. Companion tables, under the names they carry here:
`rate-sensitivity-scenarios-v3.csv`, `pooled-same-donor-pairing-scenarios-v3.csv`.
A third companion referenced by the source, `v3_cohort_label_corrections.csv`, was not
supplied and is not present.

Markings: **[direct]** (verbatim in a primary source), **[derived]** (computed here, assumptions
named), **[scenario]** (depends on an assumption the sources do not establish), **[not
identifiable]**.

Gate numbering follows the project's own pipeline naming: Gate 1 is the upstream correction step
applied before doublet scoring, Gate 2 the Scrublet step. This document does not redefine either.

---

## A. Executive verdict — the stage-alignment result

**Every rate in this adjudication is defined at a stage other than the matrix Scrublet will
receive, and no cohort provides evidence on how multiplet prevalence propagates between those
stages. `runtime_eligible = false` for all 21 grid rows.**

That is the v3 finding, and it is stronger than v2's. v2 established that no cohort has a
*validated* prior; v3 establishes that even the scenario values are stage-mismatched. A rate
computed from an author's called-cell matrix describes that matrix. The matrix Scrublet is given
has passed study QC and Gate 1, and neither the direction nor the magnitude of the prevalence
change across those steps is measured anywhere in these eight studies. Numerator and denominator
both move — cells are removed, and the removed cells are not a random sample of either singlets or
doublets.

The practical consequence: `expected_doublet_rate` cannot be set from evidence in this dataset. It
can only be *declared*, with the declaration carrying its source stage, its target stage, and the
explicit statement that the propagation assumption is unknown. The grid is a menu of declarations
with their provenance, not a set of estimates.

The four stage fields are now mandatory on every rate row:

| field | meaning |
|---|---|
| `source_stage` | the matrix or measurement the value actually describes |
| `target_stage` | the matrix it would be applied to |
| `propagation_assumption` | what must hold for the value to transfer between them |
| `runtime_eligible` | whether it may be consumed by a run without an explicit declared assumption |

---

## B. Per-cohort stage alignment

### Qiu / Chongqing / Cillo — `rho_pre_gate1_scenario`

- `source_stage`: author called-cell matrix
- `target_stage`: post-Gate-1 Scrublet input
- `propagation_assumption`: **unknown**
- `runtime_eligible`: **false**

Three candidates (0.4 / 0.8 / 1.2 % per 1,000 called cells) with the v2 labels intact: 0.4 comes
from a different chemistry (GEM-X v4), 1.2 has no vendor basis, and called cells are not recovered
cells, so even the matched-chemistry 0.8 is evaluated on a substituted denominator. Renaming the
quantity to `rho_pre_gate1_scenario` makes the mismatch explicit rather than implicit: these
describe the deposited matrix, not the matrix Scrublet sees.

Qiu additionally carries a diagnostic row: do not enter the vendor table through its "cells loaded"
column, because realized called cells are a median 34.2% of nominal load against the guides'
assumed 56.6–62.5% **[derived]**.

### Wilk — `rho_capture_pre_filter_scenario`

- `source_stage`: Seq-Well loading
- `target_stage`: post-author-heuristic, post-Gate-1 input
- `propagation_assumption`: **unknown**
- `runtime_eligible`: **false**

Scenario A 6.8% and scenario B 8.5%, both on the single matrix-level estimand
P(X ≥ 2 | X ≥ 1), differing only in whether the protocol's "~80%" is read as a per-cell efficiency
or as well-level dropout **[scenario]**. No central value. Two filters sit between the loading
density and the target matrix — the authors' own >75-genes-per-100-UMI heuristic and Gate 1 — and
neither their yield nor their doublet selectivity is reported. The published 2.37% stays a
reference row: Hughes give no denominator, and they applied it in DoubletFinder as a per-cell
percentile **[direct]**.

### COMBAT — totals and the theoretical post-demux figure

- `rho_capture_pre_filter_scenario` = 8 / 16 / 24%, `source_stage` nominal channel loading
  (50,000 cells, vendor table extrapolated), `target_stage` post-Gate-1 Scrublet input,
  propagation unknown, not runtime-eligible.
- `same_donor_pair_probability` = 1/14 = 0.0714, an **equal-mixture scenario** — per-channel donor
  proportions are not public.
- `post_perfect_demux_theoretical_scenario` = 0.57 / 1.14 / 1.71%, whose `target_stage` is a
  **hypothetical matrix immediately after perfect cross-donor demultiplexing — not the Gate 2
  input**. It assumes perfect cross-donor detection *and* no change in prevalence through study QC
  or Gate 1. **It must not be described as an actual Gate 2 residual prior.**

### Ren — residual not identifiable

Unchanged from v2 and restated: the original 0.08 is not reusable on the authors' own filtered
output; residual prevalence is **[not identifiable]** because the upstream filter's precision is
unmeasured and the denominator changes after filtering. The optional {1%, 2%, 4%} band is an
`arbitrary_robustness_scenario` — not an estimate, not derived from the authors' realized 5.1%
median removal, and to be labelled arbitrary wherever it appears.

### Liu, Stephenson

Liu: no identifiable source-stage rate at all (no loading, no recovered cells, fixed
30,000-barcode matrices). Its Σpᵢ² is a measured composition share whose target stage is the
doublet composition at demultiplexing, not a prevalence in any matrix. Stephenson: not identifiable
at capture level, since no public field maps cells to GEM well.

---

## C. Σpᵢ² — target stage stated explicitly

`pooled-same-donor-pairing-scenarios-v3.csv` now carries the four stage fields on all 58 rows. The
`target_stage` is "doublet composition at demultiplexing — NOT a prevalence in any matrix, and NOT
the Gate 2 input", and the propagation assumption records that no relation is established between
this composition share and the doublet prevalence of the post-Gate-1 matrix. Using it as a scaling
factor additionally assumes perfect cross-donor or cross-sample removal. Only Liu's 55 lanes are
measured (median 0.0704, range 0.0555–0.0864); COMBAT, Stephenson Cambridge and Cillo are
equal-mixture scenarios.

---

## D. Scrublet mechanics (unchanged from v2, restated for completeness)

`expected_doublet_rate` is the total multiplet prevalence among barcodes of the matrix passed in,
including undetectable homotypic doublets; φ_D is the detectable share; and
`overall_doublet_rate_ = detected_doublet_rate_ / detectable_doublet_fraction_` "should agree
(roughly) with `expected_doublet_rate`" **[direct]**. scanpy's wrapper defaults to 0.05 where
upstream defaults to 0.1 **[direct]**, so the value must always be set explicitly.

ρ is not a deletion quota in Scrublet: it enters through
`ld = q·ρ/r / (1 − ρ − q(1 − ρ − ρ/r))`, identical in both implementations **[direct]**.
**Each ρ changes the Bayesian score calibration and therefore the automatic threshold, which is
re-estimated by `threshold_minimum` on the transformed simulated-score distribution** — so every
rate requires a full re-score and re-threshold, not a new cut on a fixed score vector. The
simulation may be held fixed across rates via scanpy's `adata_sim` argument **[direct]**. For a
fixed simulation and neighbour graph the cell *ranking* is ρ-invariant (Spearman = 1.000000 across
ρ ∈ [0.01, 0.25]) **[derived]** — useful as a unit test, not as a shortcut. DoubletFinder, by
contrast, *is* quota-like **[direct]**; claims about one tool do not transfer to the other.

---

## E. Platelets (unchanged; numbers illustrative only)

Homotypic platelet doublets are not reliably identifiable by neighbour-based callers or by
library-size filters — Scrublet's own paper reports same-cluster doublets were "virtually
indistinguishable from singlets" **[direct]**. The single-capture panel (Qiu HC1, 50 platelet-like
cells: global 3-MAD flagged 0.0% of simulated homotypic platelet doublets, within-platelet 5.9%,
with the cut at 6,140 UMI against a largest observed platelet of 6,334) is **illustrative only**
**[derived]** and must not set any project-wide threshold. A MAD rule is admissible as an
independent QC flag and sensitivity analysis, estimated within capture × lineage with a
pre-registered minimum stratum size, never as a primary deletion rule, never tuned on outcome.

---

## F. Manuscript wording (v3)

**Prospective manuscript template — do not use until the sensitivity analysis has been
executed.** No grid has been run. Every claim below describes intended work, not
completed work.

> A common Gate 2 scoring implementation will be evaluated where compatible, while
> cohort-specific upstream demultiplexing, doublet filtering, input stages, and rate
> assumptions will be retained and explicitly recorded. Estimated recovered cells are not reported for any capture in
> these cohorts — Cell Ranger `metrics_summary.csv` files are not deposited for any of the eight
> studies — and the deposited matrices sit at heterogeneous processing stages, so the Scrublet
> `expected_doublet_rate` could not be estimated for the matrix each capture actually contributes.
> We therefore treated it as a declared sensitivity parameter rather than an estimate. Candidate
> rates were specified per platform and chemistry from the vendors' published expected-multiplet
> tables and, for the Seq-Well cohort, from Poisson well occupancy at the reported loading density;
> each candidate is reported with the processing stage it describes, the stage it was applied to,
> and the explicit statement that the relationship between the two is unknown (Supplementary Table
> S*). All results will be reported across the full grid. Because the expected rate changes the
> Bayesian score calibration and hence the automatically estimated threshold, scores and thresholds
> will be recomputed for every candidate rate; the rate will not be used as a removal quota, no
> fixed fraction of cells will be deleted, and no rate or threshold will be selected using outcome
> data.
>
> For pooled captures we report the modelled same-donor pairing probability Σpᵢ² per capture
> (measured donor composition where available, Liu median 0.070 per lane; equal-mixture assumption
> otherwise). This is a composition share under an independent-pairing assumption, not a measured
> residual rate, and the post-demultiplexing figures derived from it are theoretical values for a
> hypothetical perfectly demultiplexed matrix rather than estimates for the analysed data.
>
> Ren et al. applied Scrublet before deposit (median 5.1% of cells removed per sample). The
> residual doublet prevalence of that matrix is not identifiable from public data and the authors'
> original expected rate does not apply to their filtered output; results for this cohort are
> reported with the upstream filter recorded as a covariate and, where a rate was required, under
> an explicitly arbitrary robustness scenario.
>
> Homotypic doublets, including platelet–platelet doublets, are **not reliably identifiable or
> separately quantified** by the methods used here; their residual contribution is reported as a
> limitation rather than corrected.

Forbidden phrasings, extended in v3: "uniformly across captures" (the captures enter at different
stages, so uniform *pipeline* is defensible while uniform *correction* is not); "we removed the
expected X% of doublets"; "Scrublet under-detected relative to the expected rate"; "the filtered
matrix contains the recovered cells"; "Liu's matrices contain called cells"; "pooled-corrected"
for an uncorrected total; "residual prior" for the post-perfect-demux theoretical figures; and any
description of Σpᵢ² as a measured residual rate.

---

## G. Still requires author contact or controlled access

1. Cell Ranger `metrics_summary.csv` per capture (Qiu, Chongqing, Cillo, Stephenson, COMBAT).
2. Liu: confirmation of the forced cell count, per-lane loading, retained raw droplet matrices.
3. COMBAT: per-channel recovered/called cells pre-QC; measured per-channel donor proportions.
4. Stephenson: sample → GEM-well assignment for the 15 Cambridge pools; per-lane loading.
5. Cillo: cells loaded per lane.
6. Chongqing: any manuscript; Cell Ranger version; loading protocol.
7. Wilk: which array well count applies (85,000 vs 86,000); per-array loading verification; bead
   occupancy; what layer the "~80% capture efficiency" describes.
8. Ren: the pre-Scrublet matrices.
9. **New in v3, and the one that would change the grid's status column**: any measurement relating
   doublet prevalence before and after an ambient-correction step. Without it, `propagation_
   assumption = unknown` is permanent for every cohort, whatever the upstream evidence improves to.

## Division of labour

This document adjudicates what is scientifically defensible and sets no rate in any run script.
Since `runtime_eligible = false` on all 21 rows, the repository-side contract has a clear default:
refuse any run whose rate references a grid row without an accompanying explicit, human-declared
propagation assumption recorded in the run metadata, and record per capture the `source_stage`,
`target_stage`, `propagation_assumption`, quantity name, grid row id, and status label alongside
the value used.