# Phase 0 — COMBAT and GSE290679 labels, plus channel inventory (2026-09-14)

## COMBAT: labels found, open access, 16 events

The per-donor mortality flag is public and does not need EGA. Zenodo record
[10.5281/zenodo.6120249](https://doi.org/10.5281/zenodo.6120249) ships `CBD-KEY-CLINVAR.tar.gz`
(30 KB) containing `COMBAT_CLINVAR_for_processed.txt` — 611 sample rows × 31 columns with
**`Death28`** (0/1), `Outcome`, `Hospitalstay`, `Source`, age/sex/comorbidities, and — critically —
a `scRNASeq_sample_ID` column that joins straight onto the CITE-seq object.

Collapsed to donors with CITE-seq data: **124 donors, 16 dead by day 28** — 8 COVID-19
(`COVID_SEV` 4, `COVID_CRIT` 3, `COVID_LDN` 1), **6 sepsis**, 2 influenza. Donor IDs are in
`combat_citeseq_donor_mortality.csv`.

Note the gap versus the cohort-level table (27/116 COVID, 17/58 sepsis): most COMBAT deaths are
in donors profiled only by whole-blood RNA-seq or CyTOF, not CITE-seq. So COMBAT is worth 16
events here, not ~44.

**The counts are also open.** `COMBAT-CITESeq-DATA.h5ad` (6.4 GB) on the same Zenodo record holds
raw + normalised GEX and ADT with cluster annotations; EGA is needed only for FASTQ, which this
track does not use. This removes the DUA wait I assumed last session — COMBAT can enter Phase 1
on the same footing as the GEO cohorts.

Endpoint compatibility is good: `Death28` is exactly a 28-day all-cause mortality flag, so COMBAT
can serve either the follow-up-death primary or the 28-day sensitivity analysis without
reinterpretation.

## Revised event ledger

| Scoring | Total |
|---|---|
| Previously verified set | 66 |
| + COMBAT `Death28` | **82** |
| + GSE279452 scored at 90-day instead of 28-day (19 → 30) | **93** |

So ≥90 is reachable, but only on the 90-day scoring of GSE279452 — and only if influenza (2
events) and non-COVID ARDS (1 event) are admissible in the pooled model. Per stratum after
COMBAT: **COVID 52, sepsis 27, influenza 2, non-COVID ARDS 1**. The sepsis arm also stops being a
single-cohort result: GSE279452's share drops from 19/21 (90%) to 19/27 (70%), with 6 independent
events from a different country and care system.

## GSE290679: blocked behind dbGaP, not an email away

The CD8 companion preprint (bioRxiv 2026.05.28.728555) states patient-level metadata is in dbGaP
under **controlled access** — ICU participants `phs004258.v1.p1`, healthy controls
`phs004377.v1.p1` — and it fits 14/28/90-day mortality models on those labels without publishing
any per donor. This is a Data Access Request with institutional signatures, not an author
contact. Recommend dropping it unless a dbGaP DAR is already in scope; the matrices are excellent
(unfiltered per-sample h5) but 38 donors of unknown outcome are worth nothing to the endpoint.

Incidental finding from the same preprint, relevant to our own counts: they excluded two
GSE279452 samples — **Abd-S8 and Abd-S80 (GSM8571120, GSM8571121)** — for "ambiguous
metadata-to-matrix pairing." Both are non-fatal in our donor table, so dropping them costs 0
events (110 → 108 sepsis donors, 19 fatal at d28 either way). Independent corroboration that the
GSE279452 title-to-donor join has known soft spots.

## Channel inventory

**GSE161918 (Liu) — resolved, 56 GEM wells.** The GEM well is Batch × 10x Lane, named directly in
the GSM titles (`Batch 1 PBMC 10x Lane 1`), with paired RNA and CSP/ADT `filtered_feature_bc_matrix.h5`
files for each of 56 wells. Do **not** use the `Pool` column in
`GSE161918_allbatches.HTOandTimepointMetadata.csv.gz` as the well identity — 110 of its 213 rows
have `Pool = "misassigned"`. Instead, `GSE161918_demuxletResults.tar.gz` gives per-barcode donor
assignment for 55 of the 56 lanes (B1 lane 6 is absent), which is the authoritative
donor→well map: **14–19 donors per lane** (B1 lanes carry 19 each, the same donor set across the
batch's lanes). Inventory in `liu_gse161918_channel_inventory.csv`.

The `.best` files are **not** singlet-filtered — each holds 30,000 barcodes classified by
`DROPLET.TYPE`, and the cross-donor doublet load is heavy: median **42.5% DBL** per lane (range
27–61%), median 48% SNG and 8% AMB. So demuxlet already identifies the cross-donor doublets for
you, and the sensible split is to drop `DBL`/`AMB` barcodes on genotype evidence first and let
Scrublet handle only the same-donor doublets it can actually detect — running Scrublet per lane
on the raw barcode set and expecting it to find a ~40% doublet rate would be asking it to do a job
demuxlet does better. DecontX `batch` = lane either way, with donor identity carried as an
annotation rather than as the batch.

**COMBAT — resolved, 70 GEM wells.** `obs` in `COMBAT-CITESeq-DATA.h5ad` carries `Channel_ID`,
`Pool_ID`, `Death28`, `scRNASeq_sample_ID` and `QC_scrub_doublet_scores` directly: 10 pools
(`gPlexA`–`gPlexK`) × 7 channels = **70 wells, 14 donors each**. Inventory in
`combat_channel_inventory.csv`. Each channel carries 16–18 of the 18 annotated lineages, so
per-channel DecontX has ample label diversity — no need to fall back to pool-level batching. Like
Liu, donor is a within-well annotation and cross-donor doublets are in-pool. One naming trap:
`GEX_region` (values `A: CD4/naive/reg. T cells` … `F: pDC/HSC/platelets etc.`) is *not* a
physical lane or sort gate — all six regions appear in all 70 channels; it is the paper's
iterative-clustering region label.

### COMBAT's platelet yield is the real constraint, not its labels

Platelets are annotated (`PLT` major subset) but almost absent: **1,563 PLT cells out of 836,148**,
median **8 per donor**. By threshold:

| PLT cells per donor | donors | of which fatal |
|---|---|---|
| ≥1 | 105/124 | 15/16 |
| ≥10 | 46/124 | 10/16 |
| ≥20 | 22/124 | 7/16 |
| ≥50 | 8/124 | 3/16 |

The 16 fatal donors individually have 0, 1, 1, 1, 7, 8, 10, 13, 15, 24, 27, 30, 48, 55, 91, 109
platelet cells. COMBAT sorted PBMCs for viability before CITE-seq, and the deposited object holds
counts only for cells that passed their QC — so unlike the GEO cohorts, platelet droplets lost to
their filtering cannot be recovered by re-calling from raw matrices. Whatever per-donor platelet
threshold the primary analysis uses will decide how many of these 16 events actually carry a
platelet observation: at ≥20 cells it is 7, not 16. Worth checking against the per-donor platelet
counts in the existing cohorts before COMBAT is written into the power statement — it may be a
cohort that contributes to the cross-lineage and monocyte analyses while contributing little to
the platelet endpoint itself.

## Liu label conflict — adjudicated, keep 4

Resolved against the project's own audit trail (`deliverable4_fatal_label_audit.csv`), which holds
verbatim rows from **Liu et al. Table S6B** (`mmc6.xlsx`, PBMC samples sheet): that sheet lists 54
samples over 33 subjects and carries exactly four `Deceased` rows, one per subject —
HGR0000071 (BS-EQ-0005), HGR0000110 (BS-EQ-0013), HGR0000141 (BS-EQ-0022), HGR0000430
(BS-EQ-0044), all Critical with WHO severity 5.

The deposit's HTO sheet is a demultiplexing aid, not an outcome record: its `Status` field is a
free-text column that reads `na` for 73 of 213 rows, and it marks HGR0000110 and HGR0000430
`Alive` with no conflicting rows. Nothing in it indicates a follow-up window. **Keep 4 events for
Liu**, sourced to Table S6B, and treat the HTO sheet as authoritative only for lane/donor
structure. Ledger unchanged at 82.

## Canonical donor-label crosswalk (Phase 2 input)

`canonical_donor_labels.csv` merges every cohort with verified per-donor labels into one table:
**721 donors, 86 events, 65 donors with no usable label** (mostly healthy controls, plus 22
Stephenson "unknown" outcomes). Columns: `cohort`, `accession`, `donor_id` (deposit-facing),
`disease`, `event`, `event_definition`, `platelet_capable`, `label_source`, `join_basis`,
`channel_unit`, `notes`.

Events by stratum: COVID-19 53, sepsis 25, sepsis (paediatric) 2, non-septic ICU control 3,
influenza 2, non-COVID ARDS 1. **75 of the 86 are platelet-capable**; the 11 that are not are
9 COMBAT donors under 20 platelet cells and Ren's 2 deceased donors whose only libraries are
FACS-sorted T/B cells (P-S073, P-S074). Ren has 31 sorted-only donors in total — 2 deceased,
24 discharged, 5 controls — all flagged in `platelet_capable`.

Why this table reads 86 where the ledger reads 82: it counts the 3 non-septic ICU-control deaths
and Ren's 2 sorted-only donors, which the ledger excludes, and it omits Wilk 2020 (GSE150728)
entirely because that cohort's single fatal donor has no verified per-donor ID yet. `event` is
NA rather than 0 for controls and for donors whose outcome the source records as unknown, so
"denominator" decisions stay explicit rather than being baked in.

Three fields deliberately encode uncertainty rather than hiding it: `join_basis` marks Wilk 2021
as `INFERRED record-ID suffix convention`; `event_definition` keeps each cohort's native window
(28-day, 90-day, in-hospital, unspecified) instead of forcing a single label; and
`platelet_capable` carries the per-donor platelet constraint so a threshold change re-derives the
usable event count without touching the labels.

## Staged-file integrity

For verifying what landed on HPCC after the killed monitor — Zenodo's published md5s:

| file | bytes | md5 |
|---|---|---|
| `COMBAT-CITESeq-DATA.h5ad` | 6,409,089,483 | `87c6b1a733ea1adc37c808d9a357de74` |
| `CBD-KEY-CLINVAR.tar.gz` | 39,725 | `c5516db86341dcb1e74000e5a5d218a0` |

My local CLINVAR copy matches, so the label table above came from an intact file. `md5sum` the
staged h5ad against the value above before Phase 1 touches it — a truncated 6 GB HDF5 will open
and read `obs` fine and only fail later in the matrix blocks.
