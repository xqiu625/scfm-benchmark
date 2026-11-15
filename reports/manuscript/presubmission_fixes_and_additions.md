# Pre-submission fixes and additions

Compiled 2026-09-14 against the current, authoritative `main.tex` (1,904 lines) and
`supplementary.tex` (921 lines) — now synced identically across the local repo, HPCC, and GitHub
`origin/main`. Part A items are independently verified by direct inspection of the manuscript text
(not taken on any external document's word). Part B items point to this revision thread's own
data/code, all git-tracked and current in the local repo.

---

## Part A — Fixes to `main.tex` (confirmed real errors)

| # | Location | Fix |
|---|---|---|
| A1 | Line 461 | Replace `[repository URL]` with the real URL already at line 1419 (`https://github.com/xqiu625/302006_updategit`) |
| A2 | Line 1372 | Either run the canonical-cohort dropout rerun, or delete the "planned for the submission pass" sentence and soften the UCE-dropout-robustness claim at line 1038 to explicitly scope it to the legacy pool |
| A3 | Lines 1887, 1892, 1897 | Add real author lists to the three bibitems (`egsp2026survival`, `empiricalscfm2025cancer`, `tissueformer2025covid`) — these are exactly the three studies Table 1 is built on |
| A4 | Line 1897–1899 | Separately: the TissueFormer bibitem title ("COVID-19 severity classification with a murine-adapted Geneformer") doesn't match the real paper's abstract (population-level phenotype prediction from groups of cells, evaluated on COVID severity *and* mouse cortical-area identity) — re-verify against full text or fix the row |
| A5 | Front matter (~line 44–48) | Add corresponding-author designation, email, ORCIDs |
| A6 | Table 4 (line 989–1001) | scGPT's Δ column (+0.005) doesn't match its own AUC column (0.677 vs. baseline 0.654 implies +0.023). UCE's Δ column (+0.012) has the wrong **sign** — UCE (0.648) is below baseline (0.654), implied Δ=−0.006, but the table shows it as positive. Recompute both cells. |
| A7 | Table 2 (line 539) vs. Discussion (line 1287) | Table 2 lists UCE's corpus platelet content as "absent"; Discussion quotes and relies on "UCE 0.120%" as a known figure central to the corpus-falsification argument (63% of the between-model AUC spread). Pick one status and make both consistent — this affects a Key-Point-level claim. |
| A8 | Lines 577, 590, 598, 1268 | Deltas printed in text (Geneformer +0.064, STATE +0.052, UCE −0.066, STATE−scGPT +0.055) are each 0.001 off the values implied by Table 3's unrounded ensemble-OOF scores (0.065/0.053/−0.065/0.056) — footnote "computed on unrounded scores" or reprint to 4 decimals |
| A9 | Table 4 column header | "Patient-level AUC" doesn't say repeat-mean or ensemble-OOF — decides whether the baseline drop is −0.071 or −0.088 vs. the stated −0.072 |
| A10 | Line 1261–1266 | "exceeds the gap to the only generative model" is undefined since scGPT (0.720–0.739) sits inside the bidirectional-encoder spread (0.661–0.786) it's being compared against — rephrase or cut |
| A11 | Line 737–738 | "largest... assembled to date, and the first evaluated..." — two unsearchable superlatives; hedge to "the largest we are aware of" with a search basis, or drop |

**Checked and found NOT to hold** — don't act on this one: the claim that Limitations' demographics
statement (line 1355) contradicts a supplementary age note. No age-related Note S17 or Figure S5
exists anywhere in the current `supplementary.tex` (verified directly, zero matches).

---

## Part B — Additions, each with its data source

| Addition | Where it goes | Source script | Data directory (local, git-tracked) |
|---|---|---|---|
| **PseudobulkEN, HVG-mean, MarkerScore rungs** | New rows in/alongside Table 3 | `pseudobulk_elasticnet_patient_oof.py`, `hvg_marker_baselines_patient_oof_v2.py` | `scripts/31_donor_level_revision/results_bib_revision_tier2/clinical_utility/baseline_ladder_auc.csv`, `baseline_ladder_bootstrap_ci.csv`, `baseline_ladder_delong_vs_pca200.csv` |
| **Pseudobulk-vs-PCA200 equivalence (TOST)** | Methods/Results, alongside the PseudobulkEN row | `pseudobulk_equivalence_test.py` | `.../clinical_utility/pseudobulk_pca200_equivalence.csv` |
| **Platelet-to-T-cell composition ratio** | New row/rung, with the n=148/18-fatal caveat and the paired-vs-PCA200 comparison | `platelet_tcell_ratio.py`, `pla_t_vs_pca200_paired.py` | `.../results_bib_revision_tier2/pla_t_ratio/pla_t_ratio_per_patient.csv`, `pla_t_ratio_bootstrap_ci.csv`, `pla_t_vs_pca200_paired.csv` |
| **Resolution tool as a named Methods subsection** | New short Methods subsection; cite `hanley1982meaning` in the main text (currently bibliography-only) | `resolution_tool.py` | — (the tool itself; no CSV, it's a function) |
| **Consensus predictor result** | Report the result already promised at line 445 | `clinical_utility_and_consensus.py` | `.../clinical_utility/consensus_predictor.csv` |
| **PR-AUC / decision-curve net benefit** | Alongside existing ROC-AUC/Brier/ECE reporting | `clinical_utility_and_consensus.py` | `.../clinical_utility/clinical_utility.csv` |
| **k≥4 consensus tier (1,918 genes)** | One sentence next to the existing "531 of 3,917" (line 1130) | `consensus_attribution_signature.py` | `scripts/31_donor_level_revision/results_gene_attribution_gsea/cross_model/consensus_gene_set.csv` |
| **CCL5 contamination audit + reframe** | New paragraph near line 1136 (R6), or a Limitations sentence if not adding the full paragraph | `consensus_gene_contamination_audit.py`, `decontx_reaudit.py`, `multicelltype_decontx_reaudit.py` | `05_RESULTS/bib_revision_tier2/contamination_audit/consensus_gene_contamination_audit.csv`, `decontx_flagged_gene_correlations.csv`, `decontx_multicelltype_flagged_gene_correlations.csv` |
| **Ambient-profile / under-identification note** (if reporting the DecontX methodology at all) | Same paragraph or a Methods sentence | `run_decontx.R`, `decontx_ambient_profile.R` | `.../contamination_audit/decontx_io/`, `decontx_multicelltype/` |
| **Software versions** | New Methods sentence/table (Python, scanpy, scikit-learn, gseapy, torch) | — | needs to be pulled from each conda env (`scpeft_env`, `state_env`, `scgpt_finetune`, etc.) — not yet compiled anywhere, would need a short new check |
| **Platelet cell-calling/QC methods** | New Methods subsection (monocytes get one at line 403; platelets don't) | — | depends on the original source datasets' curation — not something this thread computed; would need writing from the existing pipeline documentation, not a data pull |
| **Accessions for Schulte-Schrepping/Bernardes/GSE279452** | Methods, cohort table | — | literature/repository lookup, not local compute |

---

## Suggested order of operations

1. Part A fixes first — small, mechanical, low-risk, and several (A6, A7) are real correctness
   problems independent of anything else in this document.
2. Part B additions next, in the order listed (baseline ladder and consensus predictor are the
   cheapest/highest-value; software versions and platelet cell-calling methods are the only two
   items that need new work rather than a data pull from something already computed).
