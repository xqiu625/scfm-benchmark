# Monocyte severity pilot — one-page readiness checklist v1 (2026-10-18)

Per `review:20261018_000122` Task 6. Pure synthesis of
`MONOCYTE_PILOT_PREREG_v1.md`, `MODEL_DECISION_v1.md`, `POSITIVE_CONTROL_EVIDENCE_v1.md`, and
`results/stephenson_monocyte_pilot_design_readiness_donor_visit_subtype_v2.tsv` — no new
computation, no literature search, no re-opening of anything those documents already mark
`confirmed`. Every number below is restated from those sources with its own origin named; none is
recomputed here. This does not authorize, and nothing below should be read as authorizing, any
real DE/masking/GSEA run — see "Current hard blocker" at the top, which is the single most
important line in this file.

## Current hard blocker (read this first)

**No admissible embedding exists yet for this pilot.** Job 30044888 (the only real scPRINT-2
small-v2 technical-pilot run so far) is `VALIDATION_FAILED` on cell identity (`
implementation-note:20261017_164353`). An identity-recovery adapter has since relabeled `obs_names` correctly
(`recovered_v3/`, real, applied, 128,503/128,503, 0 QC mismatches by id) — but that only answers
`identity_validation`. Whether `obsm["scprint_emb"]`'s rows are numerically bound to the correct
cell (`producer_row_alignment`) is still open pending Packet B, which is prepared and dry-run
validated but not yet submitted (see `claude_output/mono_<2026-10-18 handoff timestamp>.md`). **Do
not treat the recovered artifact as DE/masking/GSEA-ready until both checks have a real PASS.**

## A. Confirmed, locked — do not re-ask

| Item | Confirmed value | Source |
|---|---|---|
| Question | Within Stephenson (or Ren), do monocyte states/programs differ severe vs. mild/moderate; donor is the biological replicate | PREREG §1 |
| Main comparison (LOCKED 2026-10-15) | comparison = {Severe, Critical} AND confirmed ICU location at sampling; reference = {Mild, Moderate}; direction = comparison − reference | PREREG header, §7; `docs/decisions.md` 2026-10-15 |
| True donor count | 86 (SDRF-confirmed, not 90 patient_id strings) | PREREG §2 |
| Resolved donor-level arm counts | comparison=25, reference=52, excluded-other-severity=9 (86 total; `sdrf::C-8821` resolved into comparison per Cindy's 2026-10-15 "fewer-members" tie-break) | PREREG §6, `results/stephenson_icu_severity_grouping_resolve_mixed_visits.py` output |
| Cell-level join | 128,503/128,503 Stephenson monocyte-canonical cells matched 1:1 to the raw-file extraction; 0 unmatched/duplicate | PREREG §2 |
| Disease scope | COVID-19 only; Schulte-Schrepping/Bernardes permanently excluded as DATA (literature citation only, never reinstated as a cohort) | PREREG §3; `docs/decisions.md` 2026-09-19 |
| Severity field | `Status_on_day_collection_summary` (+ `Status_on_day_collection` for ICU), never `Worst_Clinical_Status`; outcome/fatality never used to design/select/tune | PREREG §3/§4 |
| Design type | `independent_donors_unpaired` (intercept+condition only) — a donor fixed effect would saturate this between-donor contrast and must NOT be added even where repeat visits exist. **Scope note:** this locks the MVP between-donor comparison's own model; it does NOT by itself settle the formal cross-sectional-vs-longitudinal analysis design or how donor correlation across a donor's multiple retained visits should be handled there (see B.8 for the related, still-open visit-selection question) — the same donor's different visits must never be treated as independent donors in any analysis | PREREG §6 |
| Validation unit | the donor, never the cell — any held-out split for model-selection/masking must be donor-level | PREREG §8 |
| Model direction | single-model direction (one main model + PCA/expression baseline), not a five-model benchmark — already decided | MODEL_DECISION §0/`docs/decisions.md` |
| First-technical-pilot model | scPRINT-2 `small-v2.ckpt` selected for the first technical pilot — do not ask this again. Final paper backbone is NOT locked by this; UCE is not scientifically excluded | `docs/decisions.md` 2026-10-16 "scPRINT-2 small-v2 selected for the first technical pilot only" |
| Provisional-input use for THIS technical pilot | approved — Stephenson provisional raw counts, all visits retained, Gate 1/historical DecontX scores attached as annotations only, no disease/pathway conclusions drawn from it. This is a narrower, already-closed question than final scientific QC admission (see B.2 below — do not conflate the two) | `docs/decisions.md` 2026-10-16 "Stephenson monocytes: provisional raw counts accepted for a technical pilot" |

## B. Still `pending_owner` — genuinely open, not yet decided (do not infer a default)

1. **Final cohort**: Stephenson (proposed primary — larger 3-tier severity distribution, clean
   cell-level join, real donor-identity crosswalk) vs. Ren (proposed alternate — severity values
   `UNVERIFIED`, no severity field in any locally-present Ren raw file at all). PREREG §2.
2. **QC-version acceptability for FINAL scientific admission** (distinct from the technical-pilot
   input-use approval now in table A above — do not conflate the two): this pilot's input is
   `provisional` because the platelet-level QC rebuild is unfinished. Whether a provisional-QC
   result may ever be admitted into a disease/pathway scientific conclusion, or must wait for the
   QC rebuild, is still open. PREREG §3.
3. **Time-window restriction** on `Days_from_onset` (97/102 samples numeric; 5 non-numeric are
   `Healthy`/`nan`/`Not_known`, never cast to 0) — field identified, filter not yet decided.
   PREREG §4.
4. **Whether Stephenson's own "severity" field matches the paper's own Methods definition of
   severity** — a literature read, not a metadata pull. PREREG §4.
5. **Monocyte subtype scope**: all monocytes vs. a classical/CD14-restricted pool. Both cohorts
   have real per-cohort `orig_celltype` subtype labels available (Stephenson: `CD14_mono`,
   `CD83_CD14_mono`, `CD16_mono`, `C1_CD16_mono`, `Mono_prolif`); not yet re-joined from the
   per-cohort raw files for whichever scope is chosen. PREREG §5.
6. **Site/center confound handling**: 3 real sites (Cambridge=36, Newcastle=43, Sanger=11,
   patient_id-level approximation, not re-verified at the true-donor level) — model, stratify, or
   report-only. PREREG §6b.
7. **`CV0201`/`CV0326` pool-membership conflict**: does not affect their donor-level severity
   labels (0 mismatches found across every independently cross-checked source) and does not block
   the current MVP design, but needs resolving before any pool-level dependency grouping is built.
   PREREG §6b.
8. **Within-arm visit-selection rule**: "always take the SDRF-true-donor's `Resample=='Initial'`
   sample" is PREREG §6's proposed, result-independent MVP rule — proposed, not yet confirmed by
   the owner, and not yet applied to any DE-ready artifact (only to the severity-distribution
   recount).
9. **Positive-control direction**: three literature papers disagree — Schulte-Schrepping and
   Combes report an ISG/IFN signature HIGH in mild, LOW/absent in severe; Lee reports type-I IFN
   CO-EXISTING with a TNF/IL-1β program specifically in severe classical monocytes (opposite
   direction, different cohort/sampling/subset). Stephenson's OWN published Table 7 values
   (`CD14_mono`/`CD83_CD14_mono` significantly DOWN in severe, `CD16_mono`/`C1_CD16_mono` not
   significant at the authors' threshold) support the Schulte-Schrepping/Combes direction but are
   CANDIDATE status only — the 6-group→3-group collapse script behind Table 7's own Δ/p-adj
   numbers is not in the authors' public repository, so that collapse rule is unconfirmed from the
   authors directly. Fallback already named if direction can't be frozen: reframe as "detection-
   power calibration" (does an ISG program show up at all) rather than a directional pass/fail.
   POSITIVE_CONTROL §"What this means", items 1/4.
10. **scPRINT-2 vs. UCE for the FINAL paper backbone** (narrower than it looked in the prior draft
    of this file — the first-technical-pilot model choice itself is already decided, see table A
    above, and should not be re-asked): what remains open is only the eventual PAPER backbone.
    UCE's trigger is an explicit runnability/applicability gate failure — **note, not yet
    adjudicated by the owner**: job 30044888's VALIDATION_FAILED is currently understood as an
    identity-pipeline INTEGRATION defect (engineering-remediable, see the hard blocker above), not
    a model runnability/applicability failure per se — whether this distinction holds up, and
    whether it has any bearing on the UCE trigger, is for the owner to judge once
    `producer_row_alignment` is actually closed, not assumed either way here.
11. **Masking-score method** (`MODEL_DECISION` §3's 10-step contract: donor-level dev/held-out
    split, mean-difference or L2-logistic severity direction, real masking mechanism TBD from
    scPRINT-2's own interface, signed vs. `|Δ(g)|` ranking, dual gene-universe coverage reporting)
    — a method PROPOSAL, nothing implemented, no split constructed.
12. **Train/validation/target data separation**: target cohort = `pending_owner` (item 1 above);
    model pretraining-corpus overlap with this project's own cohorts = `unknown`, never asserted
    as "no overlap"; validation split = donor-level, within target cohort, not yet built.
    MODEL_DECISION §4.

## C. `pending_evidence` — checkable facts, not yet checked (agents can investigate; not a Cindy decision by itself)

- Ren's severity vocabulary cannot be checked until Ren's own severity source file (not present
  locally) is located — PREREG §3.
- Model-vocabulary coverage of (a) the eventual final DE gene universe and (b) the genes the model
  can actually attribute a `Δ(g)` for — both distinct from the existing 96/97-symbol raw-panel
  coverage check already on record, and neither computed yet. MODEL_DECISION §3 step 8.
- Full-text/Methods/Supplementary-Table reads for Schulte-Schrepping, Combes, and Lee remain
  title/abstract-level only (Stephenson itself was already upgraded to a Methods/Fig./Table read on
  2026-10-15 — not reopened here). POSITIVE_CONTROL "What this round did NOT do".

## D. Descriptive snapshot from existing frozen small tables (not recomputed here)

From `results/stephenson_monocyte_pilot_design_readiness_donor_visit_subtype_v2.tsv` (the
multi-visit-donor subset only — 15 true donors with >1 visit, 31 data rows; NOT the full 86-donor
cohort table): per-visit monocyte subtype cell counts are highly uneven within a donor (e.g.
`sdrf::C-8821`'s D0 visit has only 10 `CD83_CD14_mono` cells vs. 1,614 at its D12 visit) — "subtype
present at a visit" is not the same as "usable sample size at that visit," already flagged in
PREREG §6's "Second update" note. Any subtype-scope or visit-selection decision (items B.5/B.8
above) should be checked against this table's real per-subtype counts, not assumed adequate from
presence alone.

## E. Literature gene-set signature track — status CORRECTED: 2 engineering gaps still open (not "nothing pending")

A prior version of this file claimed this track was fully closed. That was wrong, and the error
was caught by `review:20261018_115927`, independently re-verified against the
actual code on 2026-10-18 (no file in this track has been touched since `e7b603f`, 2026-10-17 —
the same commit the original review already covered):

- `results/extract_schulte_candidate_markers_v1.py:79-90` `extract()` returns only 9 keys. It is
  still missing `HLA_II_FUNCTIONAL_CANDIDATE_10_RESEARCHER_DERIVED_v1`, which IS present as the
  10th top-level set in the committed `gene_sets/SCHULTE_SCHREPPING_CANDIDATE_MARKERS_v1.json` —
  so the extractor cannot reproduce/verify that 10th set from source, and nothing currently checks
  it was derived correctly.
- `NONCLASSICAL_MONOCYTE_CD16_C1Q_UNION_PAPER_CLUSTER_MARKER_v1` and
  `S100_ALARMINS_UNION_PAPER_CLUSTER_MARKER_v1` are still labeled `"set_type":
  "PAPER_CLUSTER_MARKER"` in the JSON (verified directly, 2026-10-18), even though both are
  cross-cohort/cross-cluster-resolution unions implementation assistant built, not a single paper-reported cluster
  marker list — contradicting the JSON's own `_provenance` definition of that type. They should be
  `RESEARCHER_DERIVED_UNION` like the two HLA-II union sets already correctly are.

Both gaps are the same two items the FIRST review (`review:20261017_152624`)
already named; neither was fixed before this checklist was first written, and this correction
does not reopen a new literature search — it only reclassifies "already closed" back to
"still open, specific, bounded fix". The remaining items — which signature, if any, becomes the
frozen positive control (`pending_owner`, item B.9 above) and the paper's own unresolved biology
(undefined cluster-3 identity, platelet-contaminated cluster 6, S100A8/A9 not significant at this
resolution, Stephenson has no monocyte HLA-II claim at all) — are unchanged.

## What this document does NOT do

Does not decide any `pending_owner` item above. Does not run DE, masking, or GSEA. Does not treat
the recovered scPRINT-2 artifact as analysis-ready. Does not re-litigate the already-confirmed ICU
comparison, donor count, or single-model direction. Does not expand the literature-signature track
into a new search.
