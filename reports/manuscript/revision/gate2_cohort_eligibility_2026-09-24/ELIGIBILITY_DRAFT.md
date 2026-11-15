# Gate 2: per-cohort C1–C3 eligibility draft — six conditional cohorts

**Draft. No cohort is eligible. Nothing here approves a run.**

The 2026-09-24 ruling permits a scoring *design* for six cohorts, conditional on C1–C5. "Conditional"
means the design is allowed, **not** that any input qualifies. This draft records, per cohort, what
C1–C3 would require and what is actually established today.

Stephenson and Wilk are `not_run` and are not in scope here.

C4 (fixed implementation) and C5 (interpretation) are project-wide. The runner has guards and
unknown/NA semantics, but complete enforcement is not established: direct PCA/graph identity
and the exact production configuration still need audit. The eleven parameters are now approved
in `../methodology_review_2026-09-24/GATE2_PARAMETER_RULING.md`; production remains inactive.

## Status vocabulary

| status | meaning |
|---|---|
| `established` | evidence exists, is located, and says what the condition needs |
| `partial` | some evidence exists; a named gap remains |
| `not_established` | no evidence located; the condition is open |
| `contradicted` | located evidence says the condition does **not** hold |

**Nothing below is `established`.** That is the finding, not an omission.

---

## Qiu — `standard_scoring`

| | |
|---|---|
| **C1** raw counts | `not_established`. Deposited filtered matrix; no evidence located that the scored matrix is the un-normalised RNA counts from the deposit. Integrality alone is not sufficient. |
| **C2** capture unit | `not_established`. `channel_id = GSM` per the builder. Whether one GSM is one GEM well is asserted by the inventory (`GSM sample = one library`), not demonstrated. |
| **C3** parent support | `not_established`. The deposit supplies a filtered/called-cell matrix; what was removed before it, and whether the parent states survive, is not recorded. |
| **extra** | The ruling adds a chemistry conflict (GEO 3′ v3.1 vs paper 5′). It blocks any implementation that depends on a vendor rate; it does not by itself block a pure algorithmic-stress design. |
| **would satisfy C2** | the deposit or paper stating one library per GSM and one GEM well per library |

## Chongqing — `standard_scoring`

| | |
|---|---|
| **C1** raw counts | `partial`. The deposit supplies a filtered/called matrix (`barcodes.tsv` of the filtered triplet). Counts are present; un-normalised provenance is not shown, and there is **no unfiltered droplet matrix** in the deposit. |
| **C2** capture unit | `not_established`. `channel_id = GSM`; 133 GSMs, 133 donors. 1:1 is consistent with one capture per GSM but does not demonstrate it. |
| **C3** parent support | `not_established`. Upstream handling is recorded as "none reported" — which the ruling explicitly forbids reading as confirmed absence of upstream removal. |
| **would satisfy C3** | a methods statement describing the cell-calling and any filtering applied before deposit |

## Cillo — `residual_diagnostic`

| | |
|---|---|
| **C1** raw counts | `not_established`. Post-HTO input. RNA counts are present per GSM, but the scored object is reached through an h5ad join, and which matrix is scored is not pinned to an un-normalised source. |
| **C2** capture unit | `partial`, and weakest here. `channel_id = GSM`, but the join to the authors' h5ad runs through an **empirically mapped Seurat merge suffix** (`"1_1"`, `"2_1"`) recorded per GSM — not derivable. Of 48 libraries: 40 exact, 6 `PARTIAL — h5ad grouping differs from GEO library`, 2 `no h5ad cells`. **8 libraries have no suffix at all.** |
| **C3** parent support | `not_established`. Complete post-HTO stage coverage and surviving parent-state support need evidence. Removal of HTO doublets alone does not prove singlet parent states absent. Separately, unrestricted pairing does not identify the residual doublet mixture; retain `residual_null_unvalidated`. |
| **must not** | conflate HTO sample with donor; read unrestricted pairing as a calibrated within-sample residual |

## Liu — `residual_diagnostic`

| | |
|---|---|
| **C1** raw counts | `partial`. Per-lane CellRanger H5s exist (55 lanes). But the files are **CITE-seq**: RNA and ADT share the container, and the ruling requires the separation at the scored stage to be shown. Also a fixed 30,000-barcode boundary per lane of unknown mechanism. |
| **C2** capture unit | `not_established`. The lane token (`B1_10xlane1`) comes from the deposit's own GSM titles and is the unit demuxlet ran per. Reading it as one 10x chip well is **this repository's annotation**. The repo's own `channel_unit` says "56 wells" while **55** lanes exist. |
| **C3** parent support | `not_established`. Input is post-demultiplex; a donor loaded but not recovered as a singlet leaves no trace, so the parent population is not the loaded one. |
| **must not** | present any output as a same-donor rate estimate; treat SNG composition as loading composition beyond the `proxy` label |

## COMBAT — `residual_diagnostic`

| | |
|---|---|
| **C1** raw counts | `not_established`. The deposited h5ad's matrix provenance is not pinned to un-normalised RNA; `var['feature_types']` mixes Gene Expression and Antibody Capture. |
| **C2** capture unit | `partial`. `obs/Channel_ID` gives 70 channels (10 `Pool_ID` × 7 sorted regions), and published per-channel loading (50,000 cells/channel) supports the channel as *a* loading unit. That one channel is one GEM well is the builder's reading. |
| **C3** parent support | `not_established`. Multi-stage selection (demuxlet + Vireo + QC) before deposit; the surviving population is several filters downstream of loading. |
| **must not** | assign a prior by vendor extrapolation or truncation; read a score as validated residual prevalence |

## Ren — `residual_scoring_only`, never pooled

| | |
|---|---|
| **C1** raw counts | `not_established`. Provenance is **counts-level agreement only** — totals match, barcodes were never compared. |
| **C2** capture unit | `partial`. `sampleID` ↔ `gsm` is 1:1 across all 284 samples, and the inventory states "one library per GSM / sampleID". 284 libraries over 196 donors, so donor ≠ capture. Whether a library is one GEM well is asserted, not shown. |
| **C3** parent support | `not_established`. The input follows the authors' Scrublet pass and manual removal; stage completeness and parent-state support need evidence. This proves prior selection, not that all requisite singlet parent states are missing. A second pass is residual-scoring-only, not independent validation. |
| **must not** | present a second Scrublet pass as independent validation; transfer anything at barcode level without separate identity evidence |

---

## Summary

| cohort | C1 | C2 | C3 | eligible |
|---|---|---|---|---|
| Qiu | not_established | not_established | not_established | **no** |
| Chongqing | partial | not_established | not_established | **no** |
| Cillo | not_established | partial | not_established | **no** |
| Liu | partial | not_established | not_established | **no** |
| COMBAT | not_established | partial | not_established | **no** |
| Ren | not_established | partial | not_established | **no** |

**0 of 6 cohorts have an established C1, C2 and C3.** Cillo and Ren have documented upstream
selection and remain restricted in interpretation. That does not by itself establish missing
parent states: the completeness/support question and residual-calibration question are distinct.

This draft is evidence-status bookkeeping. It does not rule on whether a `partial` can be closed,
and it does not propose how to close one — those are for the reviewer.
