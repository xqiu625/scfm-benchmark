# Sensitivity-design review packet — Gate 2 doublet rate

**For independent methodological review. Prepared 2026-09-22. Read-only assembly: no rate was
selected, no runtime configuration exists, Gate 2 was not run, no SLURM job was submitted, no
production checkpoint was read, and no outcome/mortality field was consulted.**

This packet exists because an earlier internal proposal applied the 10x "0.8% multiplet per
1,000 cells" model to **called cells**. That was withdrawn: the model's independent variable is
instrument-level recovered/loading count, not called cells, post-demultiplexing cells, or
Gate 1 retained cells. Rather than substitute a different number, the five questions below are
put to review.

**Nothing in this packet is a recommendation.** Candidate structures are offered so the review
has something concrete to accept, modify or reject; the conclusions are the reviewer's.

---

## What is being decided, and what is not

| | |
|---|---|
| **In scope** | how to parameterize a scoring-only sensitivity analysis over doublet rate |
| **Not in scope** | any reliability criterion; any deletion authorization; whether cells get removed |

Gate 2 is a **scoring** pass. `call_status == "predicted_doublet"` marks a *mechanical
candidate*; a final deletion authorization additionally requires cohort eligibility, an approved
rate, a capture reliability criterion and a completed sensitivity analysis. **No
final-authorization artifact exists and none will be created here.**

## Evidence grades used throughout

| Grade | Meaning |
|---|---|
| `measured` | observed in this project's data, or in a per-capture table the depositors published |
| `published` | stated in a paper or vendor document, not re-measured here |
| `derived` | computed from other fields by a stated formula |
| `proxy` | a measured quantity standing in for a different one, with a known stage mismatch |
| `assumed` | chosen for illustration; no evidence behind it |
| `unknown` | not obtainable from anything on hand |

A quantity graded `assumed` or `unknown` must never be presented as an estimate of a real
prevalence.

---

## Source files

Living documentation: `docs/current-status.md`, `docs/execution-plan.md`, `docs/decisions.md`,
`docs/data-contracts.md`, `docs/runbooks/gate1-gate2.md`.

Evidence layer, `reports/manuscript/revision/literature_audit/T7_doublet_prior_evidence/`:
`doublet_prior_capture_evidence.csv` (605 rows), `doublet_prior_source_audit.csv` (23 sources),
`prior_evidence_verdicts.csv`, `capture_evidence_join_key.csv` (746 captures),
`gate2_source_crosswalk.csv` (790 artifacts), `tenx_multiplet_rate_tables.csv` (4 vendor
documents), `liu_demuxlet_per_lane.csv`, `ren_gse158055_capture_qc.csv`,
`rate-sensitivity-scenarios-v3.csv`, `pooled-same-donor-pairing-scenarios-v3.csv`,
`doublet-rate-adjudication-v3.md`, `memo.md`.

Every T7 row carries `runtime_eligible=false` (21/21 rate rows, 58/58 pairing rows) and none was
changed.

---

## Cohort evidence-gap matrix

Machine-readable: `cohort_evidence_gap_matrix.csv`, produced by
`scripts/32_rebuild_upfront_qc/code/build_sensitivity_review_packet.py`.

The narrative fields are a **manual adjudication** — a human read the sources and wrote down what
they say; no script can derive "the deposit does not state which side of the filter it sits on"
from a CSV.

The script additionally **revalidates an enumerated list of invariants** against the T7
artifacts, and fails rather than letting those drift. It does **not** validate everything in this
packet. What is checked:

- capture coverage **and mapping status** per cohort, summing to 746 (so Stephenson's 143
  `unmapped_to_gem_well` cannot silently become `mapped`);
- source-list row counts: 605 capture-evidence, 23 sources, 790 crosswalk artifacts;
- 21/21 rate and 58/58 pairing scenarios `runtime_eligible=false`;
- Liu as the only cohort with a measured donor composition, its 55 lanes, and its stage being
  post-demultiplex — hence graded `proxy`;
- exactly one `assumed` equal-mixture row each for COMBAT, Cillo and Stephenson;
- Ren's 284 samples with `Ncells_manual_removal` summing to 1,462,702;
- the Ren/Stephenson/Wilk classifications, matched **on each cohort's own table row** in
  `docs/current-status.md`;
- **no cohort-level recovered count in either canonical field** — `n_input_barcodes` empty
  across 790 crosswalk rows, and `estimated_or_recovered_cells` holding only agreed sentinels
  across 605 evidence rows. Scope stated deliberately: `tenx_multiplet_rate_tables.csv` also has
  a `cells_recovered` column, populated and numeric in all 54 rows, but that is the **vendor's
  published loaded → recovered reference mapping**, not a measurement of any cohort's capture,
  and it is excluded by design. The check is sentinel-based rather than `float()`-based, so
  `~5000`, `5000 cells`, `estimated 5,000` and `5k` trip it for human review instead of being
  waved through as prose;
- **Liu's 55 scenario lanes match the demuxlet source lane IDs one-for-one**, with no duplicates
  on either side — agreeing on the number 55 is not the same as agreeing on *which* 55 lanes;
- COMBAT's 50,000 loaded exceeding the maximum in every vendor table (29,000 across 4 documents).

What is **not** machine-checked: every narrative judgement, the correspondence between this
prose and the regenerated CSV, and whether `evidence_source` points at the right rows — those
fields name a file and a paper, which improves traceability but is **not row-level provenance**.

| Cohort | Loading (granularity) | Recovered | Upstream doublet handling | Donor proportions | Mapping | Class |
|---|---|---|---|---|---|---|
| **Qiu** | 15,000/sample `published`, **capture-addressable protocol value** | `unknown` | none reported `published` | n/a | 12/12 | standard candidate |
| **Chongqing** | `unknown` | `unknown` | none reported `unknown` | n/a | 124/124 | standard candidate |
| **Cillo** | `unknown` | `unknown` | HTO cross-**sample** only `published` | equal mixture `assumed` | 35/35 | upstream filter applied |
| **Liu** | `unknown` | `unknown` | demuxlet cross-**donor** only `published` | per-lane SNG composition `proxy` (post-demux stage) | 55/55 | upstream filter applied |
| **COMBAT** | 50,000/channel `published`, **capture-addressable protocol value**, beyond every vendor table | `unknown` | demuxlet+Vireo cross-**donor** `published` | equal mixture `assumed` | 70/70 | upstream filter applied |
| **Ren** | 300–14,000/reaction `published` **cohort-level only** | `unknown` | **authors' Scrublet + manual removal** `measured` | n/a | 284/284 | **residual_scoring_only** |
| **Stephenson** | 20–30k/sample `published` **protocol-level only** | `unknown` | **site-dependent** `published` | equal mixture `assumed` | **0/143 unmapped** | `not_eligible_without_additional_evidence` — **structurally not runnable** |
| **Wilk** | 15,000/array `published` **protocol-level only** | `unknown` | heuristic; **deposit side unstated** `published` | n/a | 23/23 inherited | `not_eligible_without_additional_evidence`; **non-10x platform** |

**No cohort has an obtainable recovered count.** This is the single fact that makes the vendor
model inapplicable as-is everywhere, including Qiu.

**Only Liu has measured donor proportions — but measured at the wrong stage.** Its 55 per-lane
rows come from demuxlet SNG composition, i.e. from barcodes remaining *after* DBL and AMB were
excluded. What `s = Σ p_d²` requires is the composition **at loading**. So Liu's values are
`measured` as a quantity and `proxy` as an input to `s`; they must never be described as
"measured loading proportions". COMBAT, Cillo and Stephenson have neither — a single `assumed`
equal-mixture row each, never measured at any stage.

**These are not comparable evidence.** One measured-but-stage-mismatched cohort and three
assumption-only cohorts sit in the same column of the matrix and must be read by their grade,
not by their presence.

---

## The five questions

### Q1 — Defining a sensitivity range when loading is unknown

Chongqing, Cillo and Liu have no loading figure at all. Any grid chosen for them is a statement
about *method sensitivity*, not about prevalence, and the packet must not let the two be
confused.

**For the reviewer to decide:** what principle sets the endpoints and spacing; how the grid is
labelled so a reader cannot mistake it for an estimate; whether a grid may be shared across
cohorts or must be cohort-specific; and whether "no defensible grid exists, do not score" is an
acceptable answer for any of them.

**Candidate structures, offered without preference:** (a) a fixed wide absolute grid applied
identically to every loading-unknown cohort, labelled `arbitrary`; (b) a grid anchored to the
range spanned by cohorts that *do* have loading, carried over as an explicit extrapolation;
(c) no grid — these cohorts get diagnostics only, with rate held at a single declared reference
value whose sole purpose is to make captures comparable.

**Rejection conditions to consider:** any structure that yields a single "best" rate; any
structure whose endpoints were chosen after seeing call counts.

### Q2 — Converting a recovered-stage rate to post-Gate-1 residual uncertainty

The vendor tables are indexed on cells **loaded** and report a multiplet rate among cells
**recovered**. The matrix Scrublet receives is downstream of recovery, of the depositors' own QC,
and of this project's Gate 1. Qiu is the cleanest case and still has this gap.

The chain, with each factor's status:

**Qiu as a worked example.** Published loading evidence exists for five cohorts, at two
different granularities, and the distinction matters more than its presence:

| Granularity | Cohorts | Usable as a per-capture `r`? |
|---|---|---|
| **capture-specific measured** — the actual load for this capture was recorded | **none** | — |
| **capture-addressable protocol value** — a uniform protocol figure that can be attached to each capture, but was not measured per capture | Qiu (15,000/sample), COMBAT (50,000/channel) | attachable, but it is a nominal target, not a measurement |
| **cohort/protocol-level unassignable** | Ren (300–14,000 per reaction), Stephenson (20–30k per sample; 2×25k per pool), Wilk (15,000/array) | **no** — cannot be attached to a specific capture |
| none | Chongqing, Cillo, Liu | no |

**No cohort has capture-specific measured loading.** Qiu and COMBAT are addressable because the
protocol states one figure and every capture inherits it — the between-capture variation in
actual load is unrecorded and is itself a source of uncertainty the propagation must carry.

Qiu is used below because it is the cleanest capture-addressable 10x case. The chain is the same for
COMBAT, with the extra complication that 50,000 loaded is outside the range of every vendor
table on file.

```
loading            capture-ADDRESSABLE protocol value for Qiu and COMBAT
                   (one figure every capture inherits, not a per-capture
                   measurement); cohort- or protocol-level and unassignable for
                   Ren, Stephenson and Wilk; absent for Chongqing, Cillo, Liu
  -> recovered                               UNKNOWN for every cohort
  -> vendor multiplet rate among recovered   published, but chemistry-dependent
  -> after depositor QC                      depletion factor UNKNOWN
  -> after Gate 1                            depletion factor UNKNOWN
  -> residual prevalence in the scored matrix   NOT IDENTIFIABLE
```

**For the reviewer to decide:** whether an interval can be propagated through the unknown
depletion factors at all; if so, what bounds are defensible and on what grounds; whether the
result should be reported as an interval, a grid, or explicitly as non-identifiable; and how the
Qiu chemistry conflict (GEO `3' v3.1` vs paper text `5'`) should be handled — resolved before
scoring, carried as parallel branches, or treated as a further uncertainty band.

**No numbers are filled in here deliberately.** The vendor tables are in
`tenx_multiplet_rate_tables.csv` with document, page, loaded and recovered columns; selecting
values from them is part of what is being reviewed.

### Q3 — Parameterizing upstream removal efficiency

Five cohorts had some doublet-directed filtering before our input. The formula in the existing
scenario table is the **special case `e = 1`** — complete removal of the targeted class:

```
r_residual = r · s / ((1 − r) + r · s),        s = Σ p_d²          [e = 1]
```

The general form, which the review should work from rather than re-derive:

```
                    r · [ s + (1 − e)(1 − s) ]
r_residual = ─────────────────────────────────────────
             (1 − r) + r · [ s + (1 − e)(1 − s) ]
```

with `r` the base total multiplet rate, `s = Σ p_d²` the same-source pairing probability, and
`e ∈ [0, 1]` the fraction of cross-source doublets actually removed upstream. At `e = 1` it
reduces to the expression above; at `e = 0` it returns `r` unchanged.

**Assumptions the general form carries, all of which the reviewer should rule on explicitly:**

1. doublets form by independent random pairing, with propensity independent of donor and cell type;
2. `p_d` is the composition **at loading** — not what survived demultiplexing;
3. removal efficiency `e` is **uniform across all cross-source doublets**, i.e. no donor pair is
   systematically easier or harder to resolve;
4. **no singlets were removed** by the upstream step;
5. upstream QC did not change doublet composition (it depleted classes proportionally);
6. unresolved/ambiguous barcodes are excluded from the denominator, and reassigning them would
   change `s`;
7. the denominator adjustment is the only effect of removal on the surviving population.

**Scope of the "lower bound" claim.** With `r` and `s` held fixed and assumptions 4–5 holding,
`e = 1` minimizes `r_residual`, so the existing table's values are a lower bound **within this
model**. That is not an unconditional lower bound on the real residual prevalence: Ren's manual
removal, Cillo's HTO exclusion and COMBAT's multi-stage processing may each violate assumption 4
(removing singlets) or 5 (enriching another doublet class), and either violation can move the
true residual in the opposite direction. The bound is a property of the parameterization, not of
the data.

`e` has **no measured value in any cohort.**

**For the reviewer to decide:** whether `e` enters as a scalar, an interval, or a separate grid
axis; what range is defensible given that genotype demultiplexing, HTO exclusion and an
expression-based heuristic are different mechanisms with different failure modes; what evidence
would be needed to narrow it; and whether `e = 1` may be used as a labelled bound rather than a
point value.

**Cohort-specific notes.** Liu and COMBAT: genotype demultiplexing, which also produces
unresolved/ambiguous barcodes whose reassignment changes the denominator. Cillo: HTO exclusion
removes cross-*sample* doublets only; within-sample doublets were never addressed. Ren: the
authors applied an expression-based detector on the same axis we would use, plus an undefined
manual removal — this is the case where re-running is not an orthogonal measurement.

### Q4 — Bias from donor/sample proportion proxies

`s = Σ p_d²` requires the donor composition **at loading**. What is available:

| Cohort | Available | Grade | Stage mismatch |
|---|---|---|---|
| Liu | per-lane demuxlet SNG composition, 55 lanes | `proxy` | measured on post-demultiplex singlets, after DBL and AMB were excluded |
| COMBAT | equal mixture, 14 donors | `assumed` | never measured |
| Cillo | equal mixture, 2-donor libraries | `assumed` | never measured |
| Stephenson | equal mixture, pools of four | `assumed` | never measured; pool also spans two lanes |

**For the reviewer to decide:** the direction and plausible magnitude of the bias when
post-demultiplex composition substitutes for loading composition; whether excluding AMB
barcodes biases `s` upward or downward; whether an equal-mixture assumption is admissible at all
when composition was never measured, or whether those cohorts must instead carry a widened
interval on `s`; and whether `s` should be treated as uncertain jointly with `r` rather than
plugged in as a constant.

The bound `1/n_donors ≤ s ≤ 1` holds regardless of composition, with equality on the left only
for exact equal mixture. It is available everywhere but is **too wide to be useful on its own**:
for COMBAT it admits everything from 0.071 to 1.0, a fourteen-fold range in `s` and therefore in
the same-donor residual.

**Further questions for the reviewer, beyond the bare bound:**

- Should `s` carry an interpretable scenario interval derived from plausible **donor imbalance**
  (e.g. the most/least even mixture consistent with a stated protocol) rather than the
  mathematical extremes?
- For Liu, should the AMB/DBL barcodes excluded from the SNG composition be **reassigned under
  explicit scenarios** to bound how far the post-demultiplex composition can differ from loading
  composition? `liu_demuxlet_per_lane.csv` gives per-lane **aggregate** `n_sng`, `n_dbl` and
  `n_amb` — it does **not** record which donors the AMB barcodes would have been assigned to,
  and a DBL barcode is not a singlet that can simply be handed to one donor. These aggregates
  therefore permit **bounded or explicitly assumed reassignment scenarios**, and quantify how
  large the unattributed fraction of the denominator is, but they do **not identify** the donor
  allocation and cannot on their own resolve the composition bias.
- Should `s`, `r` and `e` be varied **jointly**, given that they enter one expression and are not
  separately identifiable?
- If varied jointly, how should **implausible corner combinations** be excluded — taking the
  extreme of all three at once produces scenarios no real capture could exhibit, and reporting
  them as a range would overstate the uncertainty as much as a point estimate understates it.

### Q5 — Two-dimensional `base rate × removal efficiency`, or a wide residual-space grid

Two structures, with different identifiability and interpretability:

| | 2-D `r × e` | Direct residual grid |
|---|---|---|
| Parameters | base rate `r`, efficiency `e`, plus `s` | residual rate directly |
| Interpretability | each axis has a mechanistic meaning | one number, no mechanism |
| Identifiability | **neither `r` nor `e` is identifiable from our data**; only their combination reaches Scrublet | the quantity varied is the quantity used |
| Risk | implies the factors are separately known | discards the mechanism, so cohorts with different upstream histories look alike |
| Cost | product grid, larger | smaller |

**For the reviewer to decide:** which structure to adopt, whether it should differ between the
upstream-filtered cohorts and the untouched ones, and whether the 2-D form should be used for
*reporting* provenance while the *scored* values come from the residual axis.

---

## Cohorts whose classification this review cannot change on its own

- **Ren** stays `residual_scoring_only`. Its input is downstream of the authors' own Scrublet
  (131,962 cells removed, 7.6%) and a further undefined manual removal (136,006, 8.5%).
  Per-capture counts match `Ncells_manual_removal` for **284/284** samples — this is
  **counts-level agreement, not barcode-level**; no identity comparison was performed. The
  authors' `expected_doublet_rate = 0.08` was their *input prior*, not a measured removed
  fraction, and is not reusable as our residual.
- **Stephenson** is **structurally not runnable**, which is stronger than "needs more evidence".
  It is not that a rate is missing: 0 of 143 captures resolve to a GEM well, and upstream doublet
  handling differs by site (Scrublet at Newcastle, emptyDrops + genotype tagging at Cambridge,
  souporcell at London) with no way to attribute a capture to a site. **There is no legitimate
  unit of analysis**, so a rate sensitivity design has nothing to be defined over. Resolving
  identity and processing provenance is a precondition for the question, not a refinement of the
  answer.

- **Wilk** stays `not_eligible_without_additional_evidence`, for two independent reasons.
  (i) The deposit does not state which side of the authors' `>75 genes per 100 UMIs` multiplet
  filter it sits on. (ii) **Seq-Well is not a droplet platform**, so the 10x vendor tables do not
  apply in principle — not merely out of range as with COMBAT's 50,000, but a different physical
  process. Even if (i) were resolved, Wilk would still require either a Seq-Well-specific
  occupancy model (Poisson well loading, e.g. Hughes 2020) with its own validation, or an
  explicitly `arbitrary` sensitivity analysis. It cannot borrow a 10x rate under any label.

An existing scenario row is not evidence; only a new primary source would change either.

---

## What a reviewer is asked to return

1. For each of Q1–Q5: a decision, or an explicit "insufficient evidence, do not proceed".
2. For each cohort: whether a scoring experiment may be designed at all, and under what label.
3. Rejection conditions — what would invalidate the design after the fact.
4. Confirmation that nothing in the design depends on outcome information, on call counts
   already seen, or on a preferred retention.

Only after that does anything become a runtime configuration, and that step is a separate,
explicitly recorded decision in `docs/decisions.md`.
