# Phase 1 QC caveats

Open, unresolved data-quality notes from the upfront-decontamination rebuild that were
investigated but not fully explained, with an explicit decision on how each was handled.
See `project_302006_rebuild_architecture` memory for the full per-cohort technical writeup;
this file is the manuscript-facing summary of what's still open.

## COMBAT: B-cell contamination is a large, consistent, unexplained outlier

**Finding.** DecontX (per-channel, `z`=harmonized lineage, verified 2026-09-16 to give results
mathematically identical to celda's own multi-batch `batch=` splitting via a controlled
comparison on Ren) estimates COMBAT's B_cell lineage as ~80% contaminated on average (mean
79.9%), consistently across all 70 independently-processed channels (range 52.4-88.7%, channel-
level median 81.2%, 825-2,133 B cells/channel). This is far above every other lineage in COMBAT
(next highest: Erythroid 46.5%, Platelet 37.6%; lowest: T_cell 4.1%) and is not replicated by the
analogous smaller lineage in Liu or Ren (e.g. Liu's DC, also a minority lineage, was mid-range at
31%, not an outlier).

**Two candidate mechanisms investigated and ruled out:**
1. *Small-channel/small-N EM instability* - ruled out for Ren generally (cell-count vs.
   contamination Spearman rho=+0.223, wrong sign; flat above ~500 cells/channel) and B_cell
   channels here have 825-2,133 cells each, not small.
2. *Low native RNA content inflating the apparent contamination fraction* - ruled out directly:
   COMBAT B_cell has median 3,615 / mean 6,143 raw UMI/cell, the second-highest mean of all 9
   lineages (above NK_cell, Monocyte_Myeloid, T_cell, DC, Platelet), not lower.

**Gene-level check (the decisive evidence).** In 5 channels spanning the full contamination
range, the genes accounting for the most removed counts from B cells are universal housekeeping/
ribosomal genes (ACTB, B2M, EEF1A1, the RPL/RPS family, FTL/FTH1, mitochondrial genes) - not
canonical B-cell identity markers. Total removal attributed to a full canonical B-cell marker
panel (CD79A, CD79B, MS4A1, CD19, IGHM, IGHD, BANK1, IGKC, etc.) is only 0.6-1.5% of all removed
counts per channel. **One partial exception: CD74** (MHC-II invariant chain - real, but an
antigen-presentation accessory gene, not a core lineage-defining marker) shows substantial
removal, up to ~90% in some channels.

**Decision (user, 2026-09-16): proceed with the decontaminated data.** Canonical B-cell identity
genes are largely preserved, satisfying the pre-agreed criterion for continuing rather than
pausing (the criterion being: pause only if canonical identity genes are the main thing removed).
The overall 80% magnitude and why it's specific to B cells - not a confirmed low-RNA-content or
small-N effect - remains unresolved. Leading unconfirmed hypothesis: a model-identifiability
limitation where B_cell (a numerically smaller lineage relative to the monocyte/T-cell-dominated
ambient estimate in most channels) has a housekeeping-gene profile not well-separated from the
ambient pool those majority lineages drive. **Not pursued further per explicit decision** -
downstream robustness reruns (before vs. after decontamination, on the pre-specified PCA-200 +
one-FM validation pass) take priority over root-causing this mechanism.

**Practical implication for downstream use:** treat COMBAT B-cell-level results (any B-cell-
specific subgroup analysis, not the whole-cohort platelet/mortality endpoint) with explicit
caution until/unless revisited. This does not block the cohort's overall inclusion.

Diagnostic scripts: `code/diagnostics/combat_bcell_ambient_check.py`,
`code/diagnostics/combat_lineage_umi_check.py`.

## Cillo: 1 of 12 fatal donors (ALIR_5310_10) not recoverable from the public GSM tar files

**Finding.** `build_cillo_decontx_input.py` (2026-09-17) joins each of the 42 PBMC-only raw
per-GSM 10x tar.gz files against the deposited annotated h5ad (patient_id/outcome) via a
barcode+h5ad_suffix key, and recovers 29 of Cillo's 40 canonical donors, **11 of 12 fatal events
(92%)** - every COVID- and ARDS-named fatal donor is present. The one fatal donor NOT recovered,
`ALIR_5310_10`, together with a second (non-fatal) donor `ALIR_5097_10_PBMC_TEST`, use the
depositors' internal ALIR sample-naming scheme rather than the COVID_n/ARDS_n scheme every GSM
title and tar filename uses - project_302006_rebuild_architecture already flagged (2026-09-15)
that these two patient_ids have no obvious matching GEO title, but had incorrectly noted "no event
donor is known to be affected" - corrected 2026-09-17, `ALIR_5310_10` IS event=1.0 in
canonical_donor_labels.csv, so this is a real 1/12-event coverage gap, not a cosmetic one.

**Not pursued further for now.** Resolving it requires reading all 48 tar files' internal
member/sample directory names (not just the 42 already processed) to find where these two
patient_ids' raw counts actually live - flagged as open, deliberately deferred so as not to block
Cillo's DecontX run and the cohort's entry into the benchmark; the other 11 fatal donors and 6
excluded ETA-mixed libraries are unaffected by this gap.

## Cillo DecontX (2026-09-17) - run, no anomaly requiring the same investigation as COMBAT

35 channels (of 42 PBMC-only GSMs; 7 skipped for having no matching h5ad cells - all non-fatal
Healthy/COVID_24 donors, see build_summary.csv), 71,150 cells, 33,538-gene common panel
(Ensembl-ID intersection across channels - raw per-GSM feature files differ by up to 1 gene due to
minor 10x reference-version drift across the deposit's processing batches). Contamination: mean
11.7%, median 1.1% (much lower than Liu 29.4%/Ren 52.7%/COMBAT's overall profile - plausible given
Cillo's smaller, more homogeneous per-channel matrices). By lineage: Erythroid highest (48.9%,
same expected hemoglobin-ambient pattern seen in Liu at 79%), B_cell second (39.6%) and Platelet
third (35.7%, the project's core thesis quantity) - B_cell being elevated here too is a pattern
worth watching across cohorts, but at 39.6% it is nowhere near COMBAT's 79.9% outlier magnitude and
does not on its own warrant the same gene-level investigation; T_cell/NK_cell/Monocyte_Myeloid
lowest (2.9-4.4%). No single channel or lineage stands out as an unexplained anomaly the way
COMBAT's B_cell did.

## Platelet + Monocyte_Myeloid benchmark input assembled across all 4 Phase-1 cohorts (2026-09-17)

`assemble_platelet_monocyte_benchmark_input.py` combines the decontaminated Platelet and
Monocyte_Myeloid cells from Liu/Ren/COMBAT/Cillo (ALL canonical donors, not restricted to
within-COVID the way the validation subset was), labeled with donor_id, disease, event,
event_definition, and a per-cell draw_time whose PRECISION differs by cohort - never presented as
uniform. **Result: 789,227 cells x 16,603 genes (gene intersection across the 4 cohorts), 381
donors, 43 fatal events, 748,528 Monocyte_Myeloid + 40,699 Platelet cells.** Per cohort: Liu 33
donors/4 events, Ren 196/12, COMBAT 123/16, Cillo 29/11 (out of the previously-flagged 12, still
missing `ALIR_5310_10`). Staged at
`05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input/{raw,decontx}/` on HPCC
(counts.mtx: 14.8G raw / 33.5G decontx - too large for git; cells.csv: 107M, also HPCC-only,
summarized here instead).

**draw_time precision by cohort, verified directly against each cohort's own metadata - not
assumed uniform:**
- **COMBAT (282,704 cells, exact_day)**: per-cell `TimeSinceOnset` (numeric days since symptom
  onset) already in the deposited h5ad's obs.
- **Cillo (28,995 exact_day / 526 unresolved)**: exact integer day re-parsed from the h5ad's
  `patient_id` string (e.g. "day_5"); unresolved only for the 2 ALIR-named patient_ids and Healthy
  donors, which carry no day token.
- **Ren (315,591 cells, category)**: categorical `sample_time` (progression / convalescence /
  control) from the channel inventory - not a numeric day, but per-channel exact (Ren has no
  multi-timepoint pooling ambiguity, one library per GSM/sampleID).
- **Liu (138,657 partial_suffix / 22,754 unresolved_mixed_pool)**: GEO's own "timepoint"
  characteristic for Liu's pooled 10x lanes is literally "mixed", and demuxlet resolves donor
  GENOTYPE only (identical across a donor's repeat draws) - so per-cell timepoint is not
  reliably resolvable from first principles. In practice, however, the majority of cells (86%)
  DO carry a demuxlet-resolved donor_id with a "_T<n>" suffix (e.g. "HGR0000430_T1"), meaning
  whoever ran the depositors' own demuxlet pipeline used draw-specific reference candidates for
  most lanes/donors and only fell back to a bare donor ID (no suffix, `unresolved_mixed_pool`
  here) where that could not be resolved. Used where present, flagged separately where absent -
  trusted at the same reliability as the depositors' own demuxlet resolution, not independently
  re-verified by us.

## Phase 1 extended to Qiu, Chongqing, Stephenson, Wilk (2026-09-18)

Per explicit user request to decontaminate the remaining datasets. Schulte-Schrepping (EGA
controlled-access, confirmed directly against the EGA site) and Bernardes (accession genuinely
unresolved in the paper itself, processed data behind a FastGenomics account) are NOT attempted -
both need a user action (Data Access Request; portal registration) this session cannot take.

**Qiu/GSE167363**: 7/7 canonical donors recovered (100%), 2 fatal events, both present. DecontX
mean contamination 7.6%, median 0.7% - low, no outlier lineage or channel.

**Chongqing sepsis/GSE279452**: the originally-downloaded tar had only 108/133 canonical donors;
the missing 25 turned out to be a separate platelet-sorted (PS/PHC/PIC-prefixed) GSM batch not in
that tar, including 2 of the cohort's 22 fatal events - re-downloaded the full
`GSE279452_RAW.tar` (6.1G vs the original 4.7G) and recovered 124/133 donors, **all 22/22 fatal
events**. The remaining 9 missing donors are non-fatal (1 skin-tissue library, correctly excluded
by the PBMC-only policy; 8 platelet-sorted controls with an apparent filename-casing mismatch --
not chased further since no event is affected). DecontX: by lineage, DC highest (38.9%), B_cell
30.5%, Platelet 22.7%, T_cell lowest (6.9%) - no isolated outlier. First run hit the 8h walltime
during the silent decontXcounts()/write phase (decontX() itself had finished at 5.05h); resubmitted
at 20h partition=batch, COMPLETED in 11h17m (contamination numbers identical to the first attempt,
confirming decontX's `set.seed` determinism survived the restart). **User feedback captured for
future large jobs on this project (see feedback_hpcc_checkpoint_short_partition memory): this
single-call, non-checkpointed approach on batch/highmem should be replaced with the checkpointed
per-channel-split + short-partition array pattern already validated on Ren, rather than resubmitting
a long unchunked job when one attempt times out.**

**Stephenson/E-MTAB-10026**: 130/130 canonical donors recovered (100%), 8 fatal events. Public via
ArrayExpress/BioStudies (confirmed NOT controlled-access, unlike Schulte-Schrepping). Real
structural quirk caught before it silently dropped every cell: obs's own row index is stored under
the column name "covid_index", not the usual "_index". DecontX complete: mean contamination 13.0%,
median 2.3%. B_cell is the highest lineage at 40.1% (Erythroid 29.9% second) - elevated similarly
to COMBAT's B_cell pattern but at half the magnitude (COMBAT was 79.9%) and Stephenson's is not an
isolated single-lineage spike the way COMBAT's was (Erythroid is close behind); not investigated
further, but worth noting as a second cohort where B_cell reads relatively higher than other
lineages, which may be worth a cross-cohort look if the pattern recurs in the remaining cohorts.

**Wilk/GSE174072**: raw data is double-gzip-compressed Seq-Well RDS (exon/intron/spanning layers,
gene SYMBOL only, no Ensembl ID) -- used exon-only counts (standard non-nucleus scRNA convention).
Donor/sample resolution required parsing canonical_donor_labels.csv's `notes` column rather than a
clean GSM map; PBMC-only policy applied (ACK-only donors excluded). 23/25 donors recovered, **all
5/5 fatal events present** (both exclusions are non-fatal, ACK-only). **DecontX contamination is
much higher across the board than every other cohort** (mean 42.3%, median 38.8%, vs 7-53% mean
across Liu/Ren/COMBAT/Cillo/Qiu) - elevated fairly uniformly across lineages (Erythroid 71.4%,
B_cell 52.8% highest; T_cell 11.9% lowest), not a single isolated outlier the way COMBAT's B_cell
was. Plausible explanation already at hand, not independently confirmed: Seq-Well is a bead-based
capture platform with different (generally noisier) ambient-RNA characteristics than 10x droplet
chemistry, which is what every other cohort in this benchmark uses - a platform effect, not
necessarily a defect specific to this dataset. Not investigated further given the ready
explanation and that a systematic cohort-wide shift (unlike a single-lineage anomaly) does not
raise the same identifiability-artifact concern COMBAT's finding did.

## Final 8-cohort Platelet + Monocyte_Myeloid benchmark input assembled (2026-09-19)

Extended `assemble_platelet_monocyte_benchmark_input.py` to all 8 Phase-1-complete cohorts
(added Qiu/Chongqing/Stephenson/Wilk to the original Liu/Ren/COMBAT/Cillo) and refactored it into
checkpointed per-cohort extraction (`--cohort`) + a separate `--combine` step, run as 8 independent
short-partition jobs rather than one monolithic job -- direct application of the HPCC job-design
correction from the Chongqing DecontX timeout (see feedback_hpcc_checkpoint_short_partition
memory). Also added `read_mtx_columns()`, a chunked MatrixMarket reader that filters to the wanted
(Platelet/Monocyte) columns while streaming, since Chongqing/Stephenson's decontx matrices (44G/26G)
are too large to mmread whole before subsetting.

**Result: 1,236,881 cells x 11,709 genes, 665 donors, 80 fatal events**, all 8 cohorts. By lineage:
Monocyte_Myeloid 1,123,914, Platelet 112,967. Per cohort (donors/events): Liu 33/4, Ren 196/12,
COMBAT 123/16, Cillo 29/11, Qiu 7/2, Chongqing 124/22, Stephenson 130/8, Wilk 23/5.

**This is the FINAL scope** per Cindy's 2026-09-19 decision (see
project_302006_decisions_2026-09-19 memory and benchmark_design_decisions.md Decision 4):
Schulte-Schrepping and Bernardes are excluded from the benchmark, not deferred -- 8 cohorts / 80
events is the fixed scope going forward, not a placeholder pending future access.

Staged at `05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2/{raw,decontx}/`
(counts.mtx ~22G raw, cells.csv 180M -- both HPCC-only; small summaries committed to
channel_inventories/platelet_monocyte_v2_*.csv). One real bug caught during the refactor: the
per-cohort checkpoint writer initially forgot to write genes.csv in the success branch (only in
the "no cells found" branch) -- caught immediately when the combine step failed with
FileNotFoundError on the first affected checkpoint, fixed by writing genes.csv once per checkpoint
unconditionally rather than duplicating it in each branch.

**IMPORTANT for any future consumer of this dataset: `event` is NaN for 50,588 / 1,236,881 cells
(53 donors across Ren, Stephenson, Cillo, Qiu, Chongqing)** - donors with a genuinely undefined
outcome (e.g. Stephenson's "unknown"-outcome patients), not a missing-data artifact. Caught while
building `pca200_pooled_secondary_cv.py` (2026-09-19): a naive `.astype(int)` cast on `event`
silently turns NaN into a garbage integer instead of raising, corrupting the classification
labels with no visible error - the job had to be cancelled mid-run and restarted after adding an
explicit `event.notna()` filter before the cast. **Every future CV/classification script reading
this dataset must filter to `event.notna()` first** - do not assume every cell/donor has a
resolved label.

**Not yet done**: extraction of actual FM embeddings from this dataset (Geneformer/scGPT/UCE/
scPRINT-2/STATE, per Cindy's 2026-09-19 decision using the SUBSAMPLED strategy, not full-scale
over all 1.2M cells) - this is the assembled, labeled INPUT to that step, not the embeddings
themselves.

## cv_group wired into the secondary pooled CV protocol (2026-09-19, priority (a) done)

`pca200_pooled_secondary_cv.py` is the first script in the rebuild to pass `groups=cv_group`
(not `patient_id`) to `repeated_cv.py`'s `run_repeated_kfold`, for the grouped 10-seed x 5-fold
StratifiedGroupKFold secondary evaluation on the pooled set (LOCO, the primary evaluation, needs
no such wiring - cohort-level splits can never break a well/pool apart). Backfilled `cv_group`
into the final assembled dataset's cells.csv (joined from canonical_donor_labels.csv, 0 missing
across all 1,236,881 cells; 503 unique cv_groups across 665 donors). Directly verified (not just
assumed from sklearn's contract) that no cv_group is ever split across a fold, for 3 seeds.
Smoke-tested on the Platelet-only subset (106,911 cells after excluding undefined-event cells,
605 donors, 446 cv_groups) with the same TruncatedSVD-200+LR baseline used throughout Phase 3.
**The cv_group-purity verification itself passed cleanly** (no group split across folds, checked
independently for 3 seeds, on the corrected labels) - this is the actual property priority (a)
exists to guarantee, and it is confirmed. The illustrative 10-seed x 5-fold AUC run that follows
it timed out at 2h on `short` (fitting TruncatedSVD-200 fresh 50 times on ~85k-cell x 11,709-gene
raw sparse counts per fold is simply slow) - not resubmitted with a longer walltime, since this
was a smoke test of the wiring, not a production run: real FM embeddings will already be dense,
low-dimensional (a few hundred to ~2k dims) per-cell vectors, not raw high-dim sparse counts
needing on-the-fly SVD, so this specific slowness is an artifact of the smoke test's choice of
baseline on raw counts and is not expected to recur once real embeddings are plugged in as X.
Priority (a) is considered done on the strength of the purity verification.
