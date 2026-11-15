# Ren and Cillo channel inventories (2026-09-14)

Both gaps closed. Files: `ren_gse158055_channel_inventory.csv` (284 rows),
`cillo_gse180578_channel_inventory.csv` (48 rows). Each carries the GEM-well identity, its member
samples/donors, and the key needed to split cells into wells.

## Ren (GSE158055) — 284 wells, clean, but two chemistries

One library per GSM, no multiplexing, and the merged matrix splits cleanly: every one of the 284
GSM titles is a `sampleID` value in `GSE158055_cell_annotation.csv.gz`, which also carries
`celltype`/`majorType` per cell for all 1,462,702 cells. So `batch` = `sampleID` and nothing needs
demultiplexing.

Two things to carry into Phase 1:

- **Mixed chemistry within the cohort: 223 libraries are 10x 5′, 61 are 10x 3′** (GEO
  characteristic `single cell sequencing platform_id`). Ambient composition and gene-body coverage
  differ between the two, so platform belongs in the model as a stratum/covariate, and any
  platelet-marker or depth threshold calibrated on 5′ libraries should be re-checked on the 3′
  subset rather than assumed to transfer.
- **Platelet/megakaryocyte cells are already annotated** as `majorType == "Mega"`: 13,056 cells
  total, 106 of 284 libraries with ≥20. All 12 deceased donors have Mega cells — 581, 710, 278,
  190, 119, 115, 96, 72, 28, 24, 21, 20 — so at a ≥20-cell threshold 12/12 Ren events qualify.

This last point revises the earlier "2 deceased donors are sorted T/B only, therefore no
platelets" note: P-S073 and P-S074 do carry 24 and 28 Mega-annotated cells despite having only
FACS/MACS-sorted T/B libraries. That is almost certainly sort carry-over rather than a genuine
platelet population, so `canonical_donor_labels.csv` now flags them `caution (N Mega cells,
FACS/MACS-sorted T/B libraries only — likely sort carry-over)` and reports measured counts for
every other Ren donor instead of a categorical yes/no.

## Cillo (GSE180578) — 48 wells, and it is multiplexed after all

My earlier note called this cohort non-multiplexed. That was wrong, and the structure matters:

- 48 GEX libraries. **27 carry two samples**; 17 carry one.
- **10 libraries pool two different donors** (e.g. `COVID_11_day_10_COVID_12_day_1`,
  `COVID_2_day_10_ARDS_2_day_1`), so cross-donor doublets are real in this cohort too.
- **6 libraries mix PBMC with endotracheal aspirate in the same well** (`*_PBMC_ETA_GEX`; the
  h5ad splits them as `..._PBMC` and `..._ETA` patient_ids). Those wells' ambient soup is a
  blood + airway mixture, which is a different soup from every other well in the pooled dataset —
  worth a sensitivity check, or excluding the ETA fraction outright given the endpoint is blood.

Barcode-level donor assignment: 40 of 48 libraries match an h5ad suffix group by exact member-set
equality (the suffix on `obs._index`, e.g. `AACACGTCAATAGAGT_1_1` → group `1_1`, is the library
key). 6 are partial matches — all of them healthy-donor/ARDS control pairings where the h5ad's
grouping disagrees with GEO's library naming by one donor (GEO `Healthy_3_Healthy_4` vs h5ad
`{Healthy_4, Healthy_5}`) — and 2 libraries have no h5ad cells at all (`COVID_17_day_10_PBMC`,
`COVID_4_day_10_PBMC`, dropped in the authors' QC). No event donor is affected by the partials.

**One caveat that does touch events.** Two of the ten multi-donor libraries have no matching CITE
library in the deposit, so their hashtags are unavailable — and both pair a deceased donor with a
survivor:

| library | donors | HTO in deposit |
|---|---|---|
| GSM5465145 `COVID_23_day_1_COVID_20_day_1` | COVID_23 survived + **COVID_20 deceased** | no |
| GSM5465158 `COVID_6_day_10_COVID_18_day_1` | COVID_6 survived + **COVID_18 deceased** | no |

For those two wells the donor label — hence the event label — rests entirely on the authors'
per-barcode assignment in the h5ad, which cannot be re-derived from the deposit. Both donors in
each pair are genetically distinct, so `souporcell`-style genotype clustering on the library would
independently reproduce the split and confirm the assignment; worth doing before those cells enter
a mortality model, since a mis-assignment moves cells between a fatal and a surviving donor.

## Remaining follow-up

The per-sample GEX tars use internal study IDs (`ALIR_5091_1_5_PBMC` inside
`GSM5465125_COVID_1_day_1_day_5_PBMC_GEX.tar.gz`), and two h5ad patient_ids — `ALIR_5310_10` and
`ALIR_5097_10_PBMC_TEST` — appear under no GEO title at all. Reading the tar member names for all
48 libraries would resolve both the ALIR↔GEO mapping and the six control-pairing mismatches; it
needs ~1.4 GB of downloads, so I left it unless you want it.