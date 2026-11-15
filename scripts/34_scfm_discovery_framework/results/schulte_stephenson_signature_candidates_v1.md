# Schulte-Schrepping / Stephenson candidate marker sets — v1 (2026-10-17)

Responds to `codex_output/mono_20261017_142412.md` Sec.3 and `codex_output/mono_20261017_151151.md`
Sec.2-3. **Candidate reference material only** — independent of, and not blocking, the Stephenson
monocyte scPRINT-2 technical pilot (that remains `not_submitted`, tracked separately in
`claude_output/mono_20261017_*.md`). No HPCC run, no real DE/masking/GSEA, no cluster recomputation
was done for this task — everything below is text/JSON extraction and symbol-identity lookups only.

## What this is

- Source: a user-supplied extraction document, hash-verified at
  `codex_output/evidence/mono_20261017_151151_signatures/mono_signatures_source.md`
  (sha256 `e55ed44f61172ea2378f7a5267c5ed0351cc29db0a630fcfb5ddad129cb86dff`) — not re-derived from
  the original supplementary Excel files directly by Claude this round.
- Output: `gene_sets/SCHULTE_SCHREPPING_CANDIDATE_MARKERS_v1.json` (10 named gene sets, generated
  **programmatically** by parsing the source document's own markdown tables, not hand-copied) and
  `gene_sets/signature_symbol_mapping_v1.tsv` (per-symbol evidence for every normalization that
  went into those sets).
- Existing frozen assets (`STEPHENSON_AUTHOR_MYELOID_SIGNATURES_v1.json`,
  `HALLMARK_INTERFERON_ALPHA_RESPONSE_msigdb.json`) were **not** modified, recomputed, or
  duplicated. The pre-existing 99 (IFN) vs. 97 (Hallmark) symbol counts and their 24-symbol overlap
  were independently re-checked this round (`python3`, set intersection on the two committed JSON
  files) and match the values already recorded in `gene_sets/README.md` — this was a verification,
  not a new computation.

## Corrections applied (per the review)

1. **`GAPT` and `OGFRL1` were wrongly grouped with Excel/date contamination.** Checked against the
   live HGNC REST API (`rest.genenames.org`, 2026-10-17): both are Approved, protein-coding gene
   symbols (`HGNC:26588`, `HGNC:21378`). Correcting their identity does **not** add them to any
   candidate set — the source document's own gap list (G4) never recorded which table/cluster they
   came from, so their biological relevance here stays unresolved. See the mapping TSV.
2. **`HLA-DRB6` is a real symbol but a pseudogene** (`HGNC:4954`, locus_type=pseudogene, checked
   2026-10-17). The 11-gene HLA-II union is kept as-is for fidelity to the source tables; a
   separate 10-gene protein-coding-only candidate is offered alongside it. Neither has been chosen
   as a scientific positive control by Cindy.
3. **No blanket "strip trailing digits" rule was used.** Every digit-suffixed symbol that actually
   appears in a candidate set (10 cases: 7 in the HLA-II table, 3 in the severe CD163hi cluster
   table) was individually checked against HGNC before being normalized — see
   `signature_symbol_mapping_v1.tsv`. Symbols outside the candidate sets that the source document
   separately flagged as general Table-S4 contamination examples (pure numbers, an Excel date, a
   handful of other digit-suffixed symbols never used in any signature here) were left exactly as
   the source document already classified them — not re-litigated, per the review's own scope
   limit.
4. **`CD41` confirmed as a CD-antigen alias, not an HGNC gene symbol** — a direct HGNC symbol
   search returns zero hits; an alias_symbol search resolves it to `ITGA2B`. Not included in any
   candidate set (its only mention is in the already-excluded platelet-marker cluster 6, G6).

## The 10 candidate sets (full detail/provenance in the JSON's own `_provenance` block)

| Set | n | Type | Scope |
|---|---|---|---|
| `HLA_II_UNION_ORIGINAL_11_RESEARCHER_DERIVED_v1` | 11 | researcher union | all Schulte monocyte-subcluster HLA-II markers, both cohorts |
| `HLA_II_FUNCTIONAL_CANDIDATE_10_RESEARCHER_DERIVED_v1` | 10 | researcher union | same, minus pseudogene `HLA-DRB6` |
| `NONCLASSICAL_MONOCYTE_CD16_C1Q_UNION_PAPER_CLUSTER_MARKER_v1` | 7 | paper cluster marker | non-classical (CD16+) cluster, both cohorts |
| `SEVERE_HLA_DRLO_S100AHI_COHORT2_CLUSTER0_PAPER_CLUSTER_MARKER_v1` | 15 | paper cluster marker | Cohort 2 cluster 0 |
| `SEVERE_HLA_DRLO_S100AHI_COHORT1_CLUSTER3_PAPER_CLUSTER_MARKER_v1` | 12 | paper cluster marker | Cohort 1 cluster 3 |
| `SEVERE_HLA_DRLO_CD163HI_COHORT2_CLUSTER2_PAPER_CLUSTER_MARKER_v1` | 6 | paper cluster marker | Cohort 2 cluster 2, non-ISG markers |
| `ISG_COHORT2_CLUSTER2_PAPER_CLUSTER_MARKER_v1` | 19 | paper cluster marker | same cluster 2, ISG markers |
| `ISG_COHORT1_CLUSTER2_SINGLE_DONOR_CAVEAT_PAPER_CLUSTER_MARKER_v1` | 41 | paper cluster marker | Cohort 1 cluster 2 — **single-donor caveat, do not use alone** |
| `MILD_HLA_DRHI_CD83HI_COHORT2_CLUSTER1_PAPER_CLUSTER_MARKER_v1` | 9 | paper cluster marker | Cohort 2 cluster 1 |
| `S100_ALARMINS_UNION_PAPER_CLUSTER_MARKER_v1` | 3 | paper cluster marker | mixed clustering resolution — see warning in JSON |

## Unresolved / not decided (left exactly as the source document and review state them)

- **Cluster 3** (Cohort 2) — strong HLA-II marker cluster whose biological identity the original
  paper's main text never defines (source doc G7). Not included as a named candidate set beyond
  its contribution to the HLA-II union.
- **Cluster 6** (Cohort 2) — the paper's main text calls it "dysfunctional monocytes," but its own
  marker table top hits are platelet genes (PPBP, PF4, ITGA2B/CD41, TUBB1) (source doc G6). Excluded
  from every candidate set here.
- **S100A8/A9 at monocyte-subcluster resolution**: not significant in the subcluster marker tables.
  This means "not listed," not "not expressed" — the source document is explicit that it cannot
  determine why (comparison choice, filtering, or genuinely shared expression are all possible).
- **Stephenson monocyte HLA-II claim**: does not exist in the Stephenson supplement or main text
  (source doc G1/G2, independently re-confirmed by the source document's own second-round PMC/
  supplement re-scan). The HLA-II direction above is Schulte-Schrepping's finding only.
- **Which set, if any, becomes this project's positive control, DE weight, or scientific
  contrast**: not decided here. That is Cindy's decision, same as the existing Hallmark/
  Stephenson-author IFN sets' status in `POSITIVE_CONTROL_EVIDENCE_v1.md`.
- Mendeley-vs-cell.com file-identity, full-workbook completeness, and "independently re-verified"
  language in the source document are the *provider's* claims, not independently re-derived by
  Claude this round — only the hash of the already-saved markdown extraction, the gene-set JSON's
  internal consistency, and each normalized symbol's HGNC identity were checked by Claude.
