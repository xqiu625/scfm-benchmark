# Monocyte positive-control evidence: three source papers, independently re-checked

Per `review:20261014_204536` §3/§5.A.4, revised per `review:20261014_221155`
R3: the review assembled three candidate papers for a monocyte IFN/ISG (or HLA-II) positive
control, with specific figure/table pointers. This document independently re-verifies each
citation's identity and core claim (title/journal/year via a direct `WebSearch`, 2026-10-14 —
**not** by reading the full paper or its supplementary tables, which remains the owner's
~20-30-minute task per the review's own framing) and separates a **curated pathway collection**
(a fixed, versioned gene set, now actually downloaded and frozen — see "Offline asset" below)
from a **paper-reported signature** (that paper's own specific gene list or cluster-marker set) —
per `DEVELOPMENT.md` rule 1, neither is asserted beyond what was actually checked this round. **R3
correction**: the three papers' differing IFN findings are reported below as a documented
DISAGREEMENT this project must weigh, not asserted as a "direct contradiction under the same
conditions" — different cohorts, sampling windows, and monocyte subsets could each independently
explain part of the difference, and this document does not rule that out.

## Evidence-status upgrade (2026-10-15): Stephenson itself, now read at Methods/Fig.2f/Table-level

Per `review:20261015_104226` and `review:20261015_114223`: Cindy supplied the
Stephenson et al. *Nature Medicine* 2021 PDF and its Supplementary Excel workbook directly, and
this project has now actually read Methods (p.14–15), Fig. 2f, and Supplementary Tables 2/6/7 —
**not** title/abstract-level, unlike the three literature-reference papers below
(Schulte-Schrepping/Combes/Lee), which remain checked at title/abstract level only and are
explicitly NOT upgraded by this note. Source identity, for every claim cited in this section:

- PDF `/Users/xinruqiu/Downloads/s41591-021-01329-2.pdf`, sha256
  `f4b9f482a73627222162e1aa0622b9162ce4fa85ea9f53ebdb24080ef371ee44`.
- Supplementary Excel `/Users/xinruqiu/Downloads/41591_2021_1329_MOESM3_ESM.xlsx`, sha256
  `85cf707ad68d55f93241ebb26735fd838e995c69ed09eea9268903a118de3838`.
- Table exports (first column = original Excel row number, per
  `codex_output/evidence/20261015_stephenson_sources/README.md`):
  `codex_output/evidence/20261015_stephenson_sources/tables_v2/Supp_Table_{2,6,7}.tsv`.

**Three kinds of evidence exist for a monocyte IFN program, and they must never be collapsed into
one number**:

1. **The paper's own published cell-level score difference and adjusted p-value** (Supplementary
   Table 7, PDF Fig. 2f caption, p.4/journal p.907): `CD14_mono` Severe−Mild Δ=−0.5506,
   `p adj`=0.0210 (Excel row 35); `CD83_CD14_mono` Severe−Mild Δ=−0.8532, `p adj`=0.00275 (row 44);
   `CD16_mono` Severe−Mild Δ=−0.7145, `p adj`=0.1552 (row 62, not significant at the authors' own
   threshold); `C1_CD16_mono` (table's own internal name `CD1_CD16_mono`) Severe−Mild Δ=−0.5207,
   `p adj`=0.7099 (row 71, not significant). This is `sc.tl.score_genes` output scored against a
   **99-symbol gene list hardcoded in the authors' own notebook**
   (`myeloid/myeloid_signature_scoring_publish.ipynb`, cell 7, citing the MSigDB
   `GO_RESPONSE_TO_TYPE_I_INTERFERON` card — see
   `results/stephenson_author_code_investigation_v1.md` §1), **not** this project's pseudobulk
   log2FC, not a GSEA NES, and not the `HALLMARK_INTERFERON_ALPHA_RESPONSE` set this project has
   frozen (direct symbol-set comparison: only 24 of the 99 GO-set symbols overlap the 97-symbol
   Hallmark set — see `gene_sets/README.md`).
2. **This project's own eventual pseudobulk DE/GSEA result** on this same Stephenson monocyte
   subset — not yet computed (no real DE backend run has occurred as of this document).
3. **A candidate scFM's masking-score output** (`MODEL_DECISION_v1.md` §3) — a model-behavior
   signal, never itself evidence of transcriptional regulation or causality, and not yet computed.

**Status of the Healthy/Mild/Severe three-group collapse behind Table 7**: Table 6's cell-level
scores retain SIX peripheral labels (Healthy, Asymptomatic, Mild, Moderate, Severe, Critical,
confirmed this round by exact-Sample-ID matching: 66 monocyte-labeled Table 6 samples all match
Table 2 with no label conflicts — Healthy 12, Asymptomatic 1, Mild 14, Moderate 20, Severe 10,
Critical 9), while Fig. 2f/Table 7 show only three groups. A repo-wide, authenticated code search
of the authors' own `scCOVID-19/COVIDPBMC` repository (commit
`9f51fe9a9ffb45bf3fbd4c302fe6e3f8979f5deb`) found the gene-scoring notebook that produces the raw
per-cell `IFN_score`/`TNF_score`/`JAK_score` values, but **no committed script that performs the
6→3 label collapse or the statistical test behind Table 7's own `Δ`/`p adj` numbers** — that step
is not in the public repository (see `results/stephenson_author_code_investigation_v1.md` §2 for
the exact search terms and hit counts). Until that collapse rule is confirmed (from the authors
directly, since public code does not resolve it) or Cindy decides this project will define its own
severe-vs-mild+moderate grouping independently of how Table 7 was built, **"IFN program down in
CD14-related monocytes, Severe vs Mild" is a CANDIDATE expectation, not a frozen positive-control
target** — see "What this means" below, item 1.

**Confirmed, not assumed**: the Methods p.15 Newcastle/London-only site restriction applies to the
myeloid **differential-abundance** figure (`R_plotting_publish.ipynb` cell 10, `#Exclude Cambridge
samples`) and does **not** appear in the IFN/TNF/JAK scoring notebook's own subsetting cell — this
project must not apply a blanket Cambridge exclusion to the monocyte IFN candidate on the strength
of the Methods' differential-abundance sentence alone (full detail:
`results/stephenson_author_code_investigation_v1.md` §4).

**Calibration, not independent validation**: reproducing this same cohort's own already-published
finding is, by construction, a sanity check on this project's own pipeline against a known
result in the SAME data, never an independent biological discovery, and must never be presented as
evidence that any one method — including a candidate scFM — is superior to another. This applies
equally if the eventual pseudobulk result agrees OR disagrees with Table 7's direction.

## What "independently re-checked" means here, precisely

For each paper: confirmed via `WebSearch` (2026-10-14) that the title/journal/year/PMID match
the review's citation, and that the claimed finding is consistent with publicly indexed
abstract/summary text. This is a **title-and-abstract-level identity and direction check**, not
a read of the paper's Methods, Results, or Supplementary Tables — those require either the full
text or the owner's own read, per `review:20261014_204536` §4.

## Candidate 1: Schulte-Schrepping et al., *Cell* 2020 — literature reference ONLY, not a data source

- **Verified** (`WebSearch`, 2026-10-14): "Severe COVID-19 Is Marked by a Dysregulated Myeloid
  Cell Compartment," *Cell* 182(6):1419-1440, DOI `10.1016/j.cell.2020.08.001`. PMID `32810438`
  as cited by the review was not independently confirmed via a PubMed page fetch (PubMed blocked
  the fetch this round), but title/journal/DOI match exactly via other indexed sources.
- **Core finding, confirmed**: HLA-DR^hi CD11c^hi inflammatory monocytes carrying an
  interferon-stimulated-gene (ISG) signature were **elevated in MILD COVID-19**; severe disease
  instead showed neutrophil precursors (emergency myelopoiesis), dysfunctional mature neutrophils,
  and HLA-DR^lo monocytes. This sharpens the review's own framing ("重症出现 HLA-DR 低表达群体") —
  the independent check confirms HLA-DR^lo tracks with severe, and the ISG-bearing population
  specifically tracks with MILD, not merely "inflammatory."
- **Data-reuse status, stated explicitly per the review's own caution**: `docs/decisions.md`
  records this cohort (`EGAS00001004571`, controlled access) as **permanently dropped** from this
  project's own data on 2026-09-19, and `scripts/31_donor_level_revision/code/
  build_monocyte_canonical.py`'s own docstring independently confirms it was never part of the
  monocyte-cell-type extraction either. **This paper is cited here as a literature reference
  point only — this document does not reinstate the cohort as a data source**, and its own
  Table S4 cluster markers (if used) are a paper-specific signature, never to be treated as an
  already-validated "severity signature" without separate checking.

## Candidate 2: Combes et al., *Nature* 2021 — literature reference, consistent direction

- **Verified** (`WebSearch`, 2026-10-14): a whole-blood-preserving single-cell strategy covering
  neutrophils, monocytes, platelets, and lymphocytes; title/journal/year match the review's
  citation (exact PMID not independently re-confirmed this round — not returned by the search
  results used).
- **Core finding, confirmed**: mild COVID-19 patients showed a **coordinated interferon-stimulated
  gene signature across multiple cell types**, largely **missing in severe cases**; the ISG-positive
  monocyte/neutrophil subsets were tied to disease severity in the authors' own figure comparisons
  (mild vs. severe vs. SARS-CoV-2-negative groups). Same direction as Schulte-Schrepping
  (ISG/IFN signature associates with MILD, not severe).
- **Data-reuse status**: this is a COVID-19 cohort (`302005data02` per `docs/cohort-sources.md`,
  16 donors in the platelet-level object) that **IS already part of this project's admitted
  data** — unlike Schulte-Schrepping/Bernardes, it was never dropped. However, **no
  monocyte-cell-type extraction exists for it yet** (per `results/monocyte_pilot_feasibility_v1.tsv`
  — the monocyte-annotated universe this round is exactly {Ren, Stephenson, Qiu-sepsis}). So
  Combes is admissible as DATA in principle but not yet USABLE for a monocyte comparison without
  a new extraction step (out of this round's scope).

## Candidate 3: Lee et al., *Science Immunology* 2020 — the explicit counter-example

- **Verified** (`WebSearch`, 2026-10-14): "Immunophenotyping of COVID-19 and influenza highlights
  the role of type I interferons in development of severe COVID-19," *Sci Immunol* 5(49):eabd1554,
  DOI `10.1126/sciimmunol.abd1554`. Title/journal/year/DOI confirmed; PMID `32651212` as cited by
  the review was not independently re-confirmed via a direct PubMed fetch this round (same
  PubMed-fetch block as Schulte-Schrepping), but multiple independent indexes agree on this PMID
  for this exact title.
- **Core finding, confirmed, and it points the OPPOSITE direction from candidates 1/2**: in
  **severe** COVID-19, the type I interferon response **co-exists with** a TNF/IL-1β-driven
  inflammatory program specifically in **classical monocytes**; this overlap was **not observed in
  milder cases**. The authors argue type I IFN signaling contributes to amplifying inflammation
  in severe disease. This disagrees with the "ISG/IFN program = mild, absent in severe" direction
  from Schulte-Schrepping and Combes above — **a real, independently-confirmed disagreement in
  what each paper reports, worth weighing explicitly, but per `review:20261014_221155`
  R3's correction, NOT asserted here as a "direct contradiction under identical conditions."**
  Lee's cohort (COVID-19 + influenza comparison, GSE149689), sampling design, and the specific
  monocyte subset/timepoint examined all differ from Schulte-Schrepping's and Combes's own, and
  from Ren's/Stephenson's — any of these could independently account for part of the difference;
  this document does not adjudicate which explanation is correct, only that a real disagreement
  exists and must be weighed, not silently resolved in either direction.
- **Data-reuse status**: GSE149689, a mixed COVID-19+influenza cohort, not currently admitted or
  excluded by any decision found this round; not in the monocyte-annotated universe either way.

## What this means for choosing a positive control — a recommendation, not a decision

1. **The three literature-reference papers do not agree on direction — a disagreement to weigh,
   not a resolved contradiction — and Stephenson's own in-cohort finding (above) is a fourth,
   CANDIDATE data point, not yet a frozen target either, pending the 6→3 grouping confirmation.**
   A positive control defined as "ISG/IFN program is HIGH in mild, LOW in severe" is supported by
   Schulte-Schrepping and Combes, by Stephenson's own published CD14_mono/CD83_CD14_mono Table 7
   values (CANDIDATE status only — see the evidence-status upgrade above), but disagreed with by
   Lee's classical-monocyte finding, in a DIFFERENT cohort/sampling/comparison context. Before
   freezing a directional positive control, the owner needs to judge whether Lee's differences
   from Ren/Stephenson (this project's actual monocyte-annotated data) are enough to explain the
   disagreement, or whether this project's own pilot should be framed as "detection-power
   calibration" (does an ISG program show up at all, in either direction) rather than a
   directional pass/fail test — exactly the fallback the review itself names
   (`review:20261014_204536` §3: "若现有文献不足以支持方向性预测，就明确标为'程序检测能力的校准'").
   **Separately (per `review:20261014_221155` R3): if this project's own pilot does NOT
   detect an ISG program at all, that is NOT automatically a pipeline failure** — a real absent
   effect, a sampling-window mismatch, insufficient subset coverage, or insufficient power can
   each independently produce a null result; a null must be reported as `indeterminate`/`not
   recovered`, never silently reinterpreted as "the pipeline is broken" OR upgraded to a finding
   by relaxing the pre-specified direction after seeing the data (e.g. accepting "significant in
   EITHER direction" as success once the real direction turns out to be the Lee-like one — this
   is exactly the kind of post-hoc direction change `MONOCYTE_PILOT_PREREG_v1.md` §9 already
   guards against, restated here for the positive control specifically).
2. **Curated pathway vs. paper signature, kept separate, per the review's own instruction, now
   WITH an actual frozen asset**: `gene_sets/HALLMARK_INTERFERON_ALPHA_RESPONSE_msigdb.json`
   (fetched directly from MSigDB's own JSON export, 2026-10-14; systematic name `M5911`, 97 genes,
   CC-BY-4.0, sha256 recorded in `gene_sets/README.md`) — not a 3-4-gene hand-pick (ISG15/MX1/
   OAS1/RSAD2 are members of this set, never a substitute for it). A PRELIMINARY coverage check
   against `stephenson_monocyte_raw.h5ad`'s own 24,929-gene RAW PANEL found 96 of 97 symbols
   present (99.0%); the one gap (`WARS1`, and its older symbol `WARS` — both checked, both absent)
   is a genuine panel gap, not a naming-alias artifact. **Per `review:20261014_231259` C:
   this raw-panel number is neither of the two coverage numbers that will actually matter** — it
   is NOT coverage against the eventual FINAL DE gene universe (whatever gene set the real DE
   analysis ends up testing), and NOT coverage against whatever genes a candidate scFM can
   actually measure/attribute a masking score for (`MODEL_DECISION_v1.md` §3 step 8) — both of
   those remain to be computed once those universes exist. Per the review's own correction:
   downloading and freezing a public, already-licensed-for-reuse gene set is NOT the "live API
   provider" this project's standing boundary excludes — that boundary is about LIVE
   enrichment/annotation calls at analysis time, not a one-time, versioned, offline asset freeze.
3. **A second candidate, antigen presentation/HLA-II, is weaker evidence this round**: Schulte-
   Schrepping's HLA-DR finding is confirmed, but it is a PROTEIN-level (flow cytometry / CyTOF)
   finding in the original paper per the review's own framing, not necessarily a transcript-level
   `HLA-DRA`/`HLA-DRB1` finding in Ren/Stephenson's own scRNA-seq data — this project's own
   transcript-level evidence for this candidate was not checked this round and remains unverified
   either way.
4. **Neither candidate's direction should be picked by looking at this project's own Ren/
   Stephenson results first** — per `review:20261014_204536`'s circularity caution, the
   positive-control definition and its pass/fail rule must be frozen BEFORE this project's own
   differential-expression numbers are computed, not fit to them afterward.

## What this round did NOT do

- Did not read Schulte-Schrepping's, Combes's, or Lee's full text, Methods, or Supplementary
  Tables — those three remain at title/abstract-level verification only (unchanged from
  2026-10-14). Stephenson's own Methods/Fig.2f/Tables 2/6/7 WERE read this round (2026-10-15,
  see the evidence-status upgrade section above) — this is the one exception to this bullet.
- Did not compute any gene score, enrichment statistic, or correlation against this project's
  own expression/model data (the coverage check above counts gene-panel MEMBERSHIP only, never a
  score; Stephenson's own Table 7 Δ/p-adj values cited above are the AUTHORS' published numbers,
  not recomputed by this project).
- Did not pick a final positive-control signature or its pass/fail direction — that remains the
  owner's decision, informed by the disagreement documented above AND by the still-unresolved
  Fig.2f/Table 7 grouping gap (`results/stephenson_author_code_investigation_v1.md` §2).
- Did not re-use any already-known literature finding (Schulte-Schrepping's/Combes's/Lee's/
  Stephenson's own reported result) as if it were a NEW discovery from this project's own data —
  those findings are source material for calibrating what this project's pilot should expect to
  see, not results this project can claim credit for. Reproducing Stephenson's own finding in the
  same cohort is calibration, never independent validation.
- Did not run the authors' analysis pipeline, fetch the authors' processed data objects, or
  contact the corresponding author about the missing Fig.2f/Table 7 grouping script.
