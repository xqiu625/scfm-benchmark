# Frozen gene set assets

Per `review:20261014_221155` R3: a public, freely-licensed curated gene set, downloaded
and frozen locally as a committed asset — never a live API call at analysis time, and never
re-fetched silently (a new version would get a new filename, never overwrite this one).

## `HALLMARK_INTERFERON_ALPHA_RESPONSE_msigdb.json`

- **Source**: fetched directly, 2026-10-14, via
  `https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_INTERFERON_ALPHA_RESPONSE&fileType=json`
  (the official MSigDB JSON export, not an HTML-scraped or hand-copied list).
- **Identity, read directly from the downloaded file**: `systematicName: "M5911"`,
  `collection: "H"` (Hallmark), `pmid: "26771021"` (Liberzon et al. 2015, the Hallmark collection's
  own reference paper), 97 gene symbols.
- **License**: MSigDB states CC-BY-4.0 for this content (confirmed via the gene set's own page,
  `gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.html`, 2026-10-14).
- **File sha256**: `2a3d93f3170bc9866ced1135bdde9870f18294f7ad3ae6cb3732dc3680483455` — recompute with
  `shasum -a 256 HALLMARK_INTERFERON_ALPHA_RESPONSE_msigdb.json` and compare before trusting an
  older copy of this file.
- **PRELIMINARY coverage check against one cohort's RAW gene panel, done this round**: 96 of 97
  symbols (99.0%) found directly in `stephenson_monocyte_raw.h5ad`'s own `var` index (24,929
  genes, `gene_symbol` namespace, checked via `h5py`, read-only); the one missing symbol is
  `WARS1` (the gene's current HGNC symbol — the OLDER symbol `WARS` was also checked and is ALSO
  absent, so this is a genuine gap in this dataset's own gene panel, not a symbol-aliasing
  artifact). **Per `review:20261014_231259` C: this is NOT the coverage number that will
  matter for the actual analysis** — the final DE gene universe and whatever a candidate scFM can
  actually measure/attribute a masking score for (see `MODEL_DECISION_v1.md` §3) are two
  SEPARATE, not-yet-computed coverage numbers; this raw-panel check is a preliminary asset-sanity
  check only. Ren's/other cohorts' coverage was not separately checked this round.
- **Not yet done**: no score has been computed from this gene set against any real expression
  data; no direction (which group should be "high") has been frozen — see
  `POSITIVE_CONTROL_EVIDENCE_v1.md`. This file is an input asset only.

## This asset is NOT the gene set Stephenson et al. actually used for Fig. 2f/Table 7

Per `review:20261015_114223` §5 item 1: the authors' own public repository
(`myeloid/myeloid_signature_scoring_publish.ipynb`, cell 7, commit
`9f51fe9a9ffb45bf3fbd4c302fe6e3f8979f5deb`) hardcodes a **99-symbol gene list**, sourced from the
MSigDB `GO_RESPONSE_TO_TYPE_I_INTERFERON` card (the GO:0034340-derived set the paper's Methods
p.15 names), used with `sc.tl.score_genes`. This is a FROZEN list in the authors' own notebook,
not a live GO-annotation call at analysis time — but it is a **different asset from
`HALLMARK_INTERFERON_ALPHA_RESPONSE_msigdb.json`** above: direct symbol-set comparison gives only
24 of the 99 GO-set symbols overlapping the 97-symbol Hallmark set. The currently-frozen Hallmark
set remains usable as an independent, pre-defined candidate for "related program calibration," but
must never be described as a reproduction of what Fig. 2f/Table 7 actually scored.

## `STEPHENSON_AUTHOR_MYELOID_SIGNATURES_v1.json` (frozen 2026-10-15)

Per `review:20261015_142903` §5 item 4: the authors' own `IFN_genes`/`TNF_genes`/
`JAK_genes` lists (all three hardcoded in the same notebook cell, not only the IFN one) are now
frozen as their own committed asset — fetched directly from `raw.githubusercontent.com` at the
pinned commit above (notebook file sha256
`6f5b470f044b22b1ff8cd89f2d18707400c0b3f4b0b33dd2931f7ee5c9089211`, matching
`results/stephenson_author_code_investigation_v1.md`'s own citation), the cell's source executed
in an isolated namespace to recover the three literal Python list assignments exactly as written —
never retyped by hand, never re-queried from a live API. Counts independently confirmed:
`IFN_genes`=99, `TNF_genes`=200, `JAK_genes`=87, all internally duplicate-free; the
IFN/Hallmark 24-symbol overlap was independently recomputed from both frozen files and matches.
**Freezing this asset is not itself adopting it as this project's positive-control gene set** —
whether to use it instead of or alongside the Hallmark set remains Cindy's decision, tracked in
`POSITIVE_CONTROL_EVIDENCE_v1.md`, not made by freezing the file. Full detail:
`results/stephenson_author_code_investigation_v1.md` §1.

## `SCHULTE_SCHREPPING_CANDIDATE_MARKERS_v1.json` + `signature_symbol_mapping_v1.tsv` (frozen 2026-10-15/17)

Per `review:20261017_142412` and `review:20261017_151151`: a
user-supplied extraction of Schulte-Schrepping et al. 2020 Cell (Table S4 cluster marker tables,
hash-verified source at
`codex_output/evidence/mono_20261017_151151_signatures/mono_signatures_source.md`,
sha256 `e55ed44f61172ea2378f7a5267c5ed0351cc29db0a630fcfb5ddad129cb86dff`) was parsed programmatically
into 10 candidate gene sets, explicitly split into `PAPER_CLUSTER_MARKER` sets (one specific
cluster's own differential-expression marker list, tied to that cohort/cluster, NOT this project's
own disease contrast) and `RESEARCHER_DERIVED_UNION` sets (a union implementation assistant built across multiple
cluster tables; the original paper never computed this exact combined list). Every digit-suffixed
("Excel-dedup-looking") symbol that made it into a candidate set was individually re-checked
against the live HGNC REST API before being normalized — no blanket "strip trailing digits" rule
was applied — see `signature_symbol_mapping_v1.tsv` for the full per-symbol evidence, including two
symbols (`GAPT`, `OGFRL1`) that a prior pass had wrongly grouped with genuine Excel/date
contamination: both are real, HGNC-Approved, protein-coding gene symbols, just not included in any
candidate set here because the source document never recorded which cluster/table they came from.
`HLA-DRB6` is kept in the paper-faithful 11-gene HLA-II union but is a confirmed HGNC pseudogene
(`HGNC:4954`), so a separate 10-gene protein-coding-only candidate variant is offered alongside it
— neither replaces the other, and neither has been selected as this project's positive control.
**This entire file is CANDIDATE REFERENCE MATERIAL — not reviewed or adopted by Cindy for any
scientific use.** Full detail and the open/unresolved items the source document itself flags (3
clusters whose biological identity the original paper leaves undefined or disputed, a single-donor
ISG cluster, S100A8/A9 not reaching significance at monocyte-subcluster resolution):
`results/schulte_stephenson_signature_candidates_v1.md`.
