# Findings from the donor-level revision thread not in `main.tex` / `supplementary.tex`

Verified directly against `reports/manuscript/main.tex` (1,904 lines) and `supplementary.tex`
(921 lines), the current, authoritative versions — now synced identically across the local repo,
HPCC, and GitHub `origin/main` (2026-09-14; local was ahead by 6 unpushed commits, HPCC was stale
since 2026-08-06 with no divergent content, GitHub was behind — all three now match commit `ef308f6`).

## Genuinely absent from the manuscript

Zero occurrences of any of these terms/results in `main.tex`, confirmed by direct grep:

| Missing | This thread's result | Status |
|---|---|---|
| `pseudobulk`, `elastic net` | PseudobulkEN ensemble-OOF = 0.7421, bounded to ±0.08 of PCA-200 (TOST fails) | computed, not written up |
| `HVG-mean` | 0.4902 — at chance; cleanest evidence the paper's nulls are real nulls | computed, not written up |
| `MarkerScore` | 0.5894 — 95% CI [0.466, 0.707] includes 0.5, not "above chance" | computed, not written up |
| `Pla-T` / platelet-to-T ratio | 0.7406, 95% CI [0.601, 0.860], n=148/18-fatal | computed, not written up |
| `DecontX`, `ambient`, `contamination`, `doublet` | Full contamination audit + DecontX (platelet-only and multi-celltype reruns) on the 13 flagged consensus genes | computed, not written up — §R6 (line 1136) makes the CCL5 claim with no accompanying audit anywhere in the manuscript |
| `TOST`, `equivalence` | The equivalence-testing framework behind the PseudobulkEN claim | computed, not written up |
| `PR-AUC`, `decision curve` | Clinical-utility layer (PR-AUC, net benefit) on the existing OOF scores | computed, not written up |
| `resolvable_delta_auc` | The named resolution tool itself | **the tool's OUTPUT is quoted repeatedly** (abstract, Key Points, Results, Discussion — "~0.10 AUC resolvable at 25 events") but the tool is never named or described in Methods, and `hanley1982meaning` sits in the bibliography uncited (`\cite{hanley1982meaning}` appears zero times in the text) |
| k≥4 consensus tier (1,918 genes) | Only the k=5 tier (531 genes, line 1130) is in the manuscript | partially in |

## Already in the manuscript (do not re-add, integrate against these)

- **Table 3** (`tab:within_covid`, ~line 680–710) — the patient-level 240/25 within-COVID headline
  table this whole thread's baseline ladder is built to sit alongside: PCA-200+LR ensemble-OOF =
  **0.742**, QC-metrics+cell-count (as-fit) = **0.618**. This is the correct table to extend with the
  new baseline-ladder rungs.
- **Consensus predictor**, line 445 — mentioned as computed, result never reported. This thread's
  consensus-predictor code (`clinical_utility_and_consensus.py`) already answers what's missing here.
- **k=5 consensus gene set**, line 1130 — "531 of the 3,917" already in text.
- **CCL5 cross-model attribution**, line 1136 — "CCL5 is negatively attributed in all five [models]"
  already in text, with no contamination caveat attached.

## Recommended integration points

1. **Baseline ladder → Table 3.** Add PseudobulkEN, HVG-mean, MarkerScore, and Pla-T as additional
   rows/rungs in or alongside Table 3, each with its resolution statement (resolvable vs. not, per
   the existing framework at line ~123/211).
2. **Resolution tool → new short Methods subsection.** Name `resolvable_delta_auc()`, give the
   Hanley-McNeil variance formula and the measured inter-model score correlation ρ it uses, and cite
   `hanley1982meaning` from the main text (currently only in the bibliography).
3. **CCL5 (line 1136) → add the contamination-audit caveat/reframe.** The DecontX work shows CCL5's
   raw mortality correlation survives multi-lineage-informed decontamination (p=0.003→0.037) while
   its depth-artifact component is eliminated (ρ 0.388→-0.005) — this *strengthens* the line-1136
   claim rather than undermining it, and is a one-paragraph addition (Results or Discussion) plus,
   if the R6 attribution section stands alone, a one-sentence Limitations note that raw-count
   attributions have no ambient-RNA/doublet correction elsewhere in the object.
4. **Clinical utility → PR-AUC/decision curve.** Add alongside the existing ROC-AUC/Brier/ECE
   reporting; ~20 lines of code on data already in hand.
5. **Consensus predictor (line 445) → report its result**, or cut the clause if not being added this
   round.

## Not independently verified — from external review documents, flagged not confirmed

A chain of external "review" documents supplied earlier in this session raised additional specific
claims (Table 4 Δ-column arithmetic mismatch, Table 2 vs. Discussion UCE-percentage inconsistency,
several other numerical/wording issues) against this same, now-confirmed-correct `main.tex`. Given
this session's track record with documents from that chain (several turned out to cite a stale/wrong
file; at least one turned out to be accurate once checked against the right file), none of those
specific claims have been independently re-verified line-by-line in this pass. Treat them as
leads to check, not confirmed findings, before acting on them.
