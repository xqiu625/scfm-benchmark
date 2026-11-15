# Review-response statistics (stage 30)

Addresses the peer-review consensus that **significance was overstated** because DeLong
was computed on cell-weighted OOF predictions (~42k correlated cells from 322 donors
treated as independent). Every headline comparison is recomputed with the **donor** as the
unit of inference. Run: `code/donor_clustered_reanalysis.py` (donor cluster bootstrap B=5000,
seed 42) on the COVID+Sepsis OOF predictions. No new model training.

Roadmap items addressed: **A** (donor uncertainty), **B** (ranking distinguishability),
**C** (cross-disease CI overlap), **L** (multiple-comparison FDR), **F** (severe-label sensitivity).

## Headline result (good news + required softening)

> **STATE's lead is robust; the precise middle-pack ordering is not.**

- **STATE significantly beats all four other FMs** on within-disease OOF AUC under
  donor-clustered bootstrap **after BH-FDR** (q ≤ 0.017). This survives the proper inference.
- **scPRINT-2 ≈ Geneformer ≈ scGPT are statistically indistinguishable** (donor bootstrap
  p = 0.37–0.58; CIs cross 0). The published ordering "scPRINT-2 > Geneformer > … > scGPT"
  is **noise** and must be presented as a tie.
- **4 of 10** cell-level DeLong "p < 10⁻³" pairs become **non-significant** once clustered by donor.

## A — per-FM donor-clustered AUC CI (LR, COVID+Sepsis OOF)

| FM | AUC | 95% CI (donor) | width |
|---|---|---|---|
| STATE | 0.732 | [0.625, 0.829] | 0.204 |
| scPRINT-2 | 0.679 | [0.572, 0.785] | 0.213 |
| Geneformer | 0.670 | [0.560, 0.776] | 0.217 |
| scGPT | 0.653 | [0.550, 0.748] | 0.199 |
| UCE | 0.604 | [0.468, 0.736] | 0.268 |

CIs are **wide (~0.20–0.27)** — the effective n is 322 donors, not 42,095 cells.

## A+L — pairwise ΔAUC: donor CI vs cell-DeLong, BH-FDR

| Pair | ΔAUC | donor 95% CI | donor p | cell-DeLong p | BH-q | sig? |
|---|---|---|---|---|---|---|
| STATE − scPRINT-2 | 0.056 | [0.013, 0.104] | 0.008 | ~0 | 0.017 | **yes** |
| STATE − Geneformer | 0.065 | [0.027, 0.105] | <0.001 | ~0 | <0.001 | **yes** |
| STATE − scGPT | 0.079 | [0.043, 0.112] | <0.001 | ~0 | <0.001 | **yes** |
| STATE − UCE | 0.126 | [0.072, 0.187] | <0.001 | ~0 | <0.001 | **yes** |
| scPRINT-2 − UCE | 0.069 | [0.018, 0.124] | 0.007 | ~0 | 0.017 | **yes** |
| Geneformer − UCE | 0.061 | [0.015, 0.107] | 0.011 | ~0 | 0.019 | **yes** |
| scPRINT-2 − Geneformer | 0.008 | [−0.023, 0.039] | 0.584 | 0.007 | 0.584 | no |
| scPRINT-2 − scGPT | 0.023 | [−0.032, 0.073] | 0.372 | 2e-12 | 0.466 | no |
| Geneformer − scGPT | 0.015 | [−0.034, 0.061] | 0.546 | 3e-5 | 0.584 | no |
| scGPT − UCE | 0.048 | [−0.014, 0.109] | 0.138 | ~0 | 0.197 | no |

**Defensible partial order:** STATE > {scPRINT-2, Geneformer} > UCE; scGPT is tied with the
middle pack and with UCE. Note especially **scGPT vs UCE is NOT significant** (the manuscript
currently orders scGPT above UCE).

## B — fold-mean t-CI + Nadeau-Bengio (the right test for CV folds)

| FM | fold mean | t-CI (5 folds) |
|---|---|---|
| scPRINT-2 | 0.762 | [0.575, 0.949] |
| Geneformer | 0.738 | [0.581, 0.895] |
| STATE | 0.724 | [0.576, 0.871] |
| scGPT | 0.641 | [0.531, 0.752] |
| UCE | 0.598 | [0.403, 0.793] |

All fold-mean CIs **overlap massively**. Nadeau-Bengio corrected paired test: the **only**
distinguishable pair is **STATE > UCE (p = 0.021)**; every other pair NS (incl. scPRINT-2 vs
STATE p = 0.74). → The per-fold-mean ranking that opens the Results is not statistically real.

## C — cross-disease CI overlap (canonical C002 sources)

- **Forward (COVID→Sepsis):** STATE / UCE / scPRINT-2 / Geneformer all mutually overlap
  (~0.78–0.86). **Only scGPT (0.354, CI [0.217, 0.503]) is non-overlapping below every FM** →
  the scGPT anti-correlation is the **one robust cross-disease claim** (keep it; note the CI
  upper bound 0.503 nears chance).
- **Reverse (Sepsis→COVID):** STATE (0.701) is highest and its CI excludes **scPRINT-2** (0.542)
  — but **overlaps UCE (0.643), Geneformer (0.612), scGPT (0.545)**.

→ "STATE is the **only** bidirectional generalizer" is **not supported as a categorical claim**:
STATE has the highest reverse transfer and significantly exceeds scPRINT-2, but its advantage
over UCE/Geneformer is within the bootstrap CI. Soften to: *"STATE shows the strongest reverse
transfer, the only FM significantly above scPRINT-2; differences vs UCE/Geneformer are within CI."*

## F — severe-negative label sensitivity (eval-side)

| FM | severe **included** | severe **excluded** | Δ |
|---|---|---|---|
| STATE | 0.732 | 0.782 | +0.050 |
| scPRINT-2 | 0.679 | 0.723 | +0.044 |
| Geneformer | 0.670 | 0.720 | +0.050 |
| scGPT | 0.653 | 0.715 | +0.062 |
| UCE | 0.604 | 0.676 | +0.072 |

Excluding the ~18k severe negatives makes the task **easier** (+0.04–0.07 AUC for every FM);
**STATE's lead is preserved either way.** The headline tables use the **harder severe-included**
contrast. (Eval-side only — a full severe-excluded *retrain* is the definitive version; this
bounds the effect and shows the ranking is label-robust.)

## Implied manuscript changes (for the revision)

1. Replace cell-level DeLong p-values with **donor-clustered ΔAUC CIs** (table above); state
   that 4/10 pairwise "significant" results do not survive donor clustering + FDR.
2. **Keep "STATE leads"** — now rigorously supported (FDR-significant over all 4 FMs) — but
   present scPRINT-2/Geneformer/scGPT as a **statistical tie**, and stop ordering scGPT above UCE.
3. Soften the per-fold-mean ranking sentence: top group **indistinguishable** (NB: only STATE>UCE).
4. Soften **"STATE is the only bidirectional generalizer"** per C; keep the scGPT-0.352
   anti-correlation as the one robust cross-disease finding.
5. Add the **severe-label sensitivity** (+0.05, ranking-robust) and state the headline uses the
   harder severe-included contrast.

Evidence CSVs: `results/{donor_bootstrap_auc_ci, donor_bootstrap_delta_auc, fold_mean_ci,
nadeau_bengio_paired, cross_disease_ci_overlap, severe_sensitivity}.csv`.

---

# Additional roadmap analyses (HPCC batch, 2026-06-09)

## Seed-stability (Devil's-Advocate alt #1: "rankings are CV-seed noise")
8 CV seeds × 5 FMs, within-disease LR. `code/seed_stability.py`, `results/seed_stability.csv`.

| FM | mean AUC across seeds | range | rank-1 freq |
|---|---|---|---|
| **STATE** | **0.724 ± 0.033** | 0.654–0.765 | **6/8** |
| scPRINT-2 | 0.685 ± 0.043 | 0.612–**0.762** | 2/8 |
| Geneformer | 0.680 ± 0.028 | 0.641–0.738 | 0/8 |
| scGPT | 0.636 ± 0.015 | 0.609–0.660 | 0/8 |
| UCE | 0.620 ± 0.021 | 0.594–0.658 | 0/8 |

Mean Kendall τ across seed-pairs = **0.571** (rankings reshuffle substantially). STATE is
rank-1 in 6/8 seeds; **scPRINT-2's headline 0.762 is the top of its range at the published
seed 42** (its mean is 0.685). → Replace "scPRINT-2 achieved the highest per-fold mean" with
"STATE has the highest mean AUC across seeds; scPRINT-2's lead is seed-42-specific."

## F-full — severe-excluded RETRAIN (definitive, not eval-side)
`code/severe_excluded_retrain.py`, `results/severe_excluded_mortality.csv`. Negatives =
mild/moderate/recovered (severe excluded from train+test), ~23,899 cells.

| FM | severe-excluded retrain AUC |
|---|---|
| **STATE** | **0.823** |
| scPRINT-2 | 0.782 |
| UCE | 0.770 |
| Geneformer | 0.763 |
| scGPT | 0.741 |

STATE is #1 under the retrain too → **its lead is label-robust**. The middle order shifts
(UCE rises above Geneformer), reinforcing that only STATE's top rank is stable.

## J — zero-shot Geneformer-V1-10M control (**reverses a manuscript claim**)
`code/geneformer_v1_zeroshot.py`, `results/geneformer_v1_zeroshot_mortality.csv`. Same
backbone / cohort / labels / 5-fold patient CV as the LoRA run; frozen base, mean-pooled.

| | AUC |
|---|---|
| Zero-shot Geneformer-**V1-10M** (new control) | **0.592 ± 0.094** |
| LoRA Geneformer-**V1-10M** (committed) | 0.705 ± 0.135 |
| Zero-shot Geneformer-**V2-316M** (manuscript comparator) | 0.699 |

The manuscript's "LoRA adds only +0.006 → domain adaptation is a null" compared LoRA-V1-10M
against a **30× larger** zero-shot backbone (V2-316M). **Within the same backbone, LoRA gives
+0.113** (0.592→0.705). The "+0.006 null" was capacity-confounded (confirms R1 W8 / R2 M4).
**Revised claim:** LoRA-adapting the small V1-10M (+0.113) brings it level with the 30× larger
zero-shot V2-316M (0.705 vs 0.699) — domain adaptation *is* effective; it just doesn't beat
simply using a bigger pretrained model. → **`tab:domainadapt` / abstract / Conclusion 6 need
revising** (this is the one analysis that changes a conclusion, not just softens it).

## H + N — cross-disease subspace alignment (the geometric test the authors proposed)
`code/subspace_alignment.py`, `results/subspace_alignment.csv`.

| FM | cos_1d (COVID·Sepsis dir) | subspace align | fwd transfer AUC |
|---|---|---|---|
| STATE | 0.169 | 0.112 | 0.798 |
| UCE | 0.098 | 0.142 | 0.807 |
| scPRINT-2 | −0.023 | 0.026 | 0.855 |
| Geneformer | 0.088 | 0.053 | 0.773 |
| scGPT | +0.208 | **0.243** (highest) | **0.352** (worst) |

- **H (null):** subspace alignment is uniformly low and **does not predict transfer** — scGPT
  has the highest alignment yet the worst transfer. The proposed subspace-alignment *mechanism*
  for the cross-disease asymmetry is **not supported** → keep the two-axis story *descriptive*,
  not mechanistic (confirms R2 M2 / R3 M2 / DA).
- **N:** scGPT's 0.352 is **not** a sign-flip (cos_1d = +0.208, positive) nor a vocab artifact
  (99.8% gene overlap) → the anti-correlation is a genuine multidimensional mismatch; keep it.

## K — matched-dimensionality ablation (dimensionality-confound test, R2 C1)
`results/ablation_summary_fullpool.csv` (full pool; existing `tab:ablation` cohort).

| FM | full-dim baseline | **PCA-256 (matched)** |
|---|---|---|
| **STATE** | 0.711 | **0.691** |
| Geneformer | 0.702 | 0.678 |
| scPRINT-2 | 0.682 | 0.664 |
| UCE | 0.653 | 0.650 |
| scGPT | 0.649 | 0.649 |

On the full pool, **STATE leads at matched 256 dims just as at full dim** → its advantage is
not a dimensionality artifact (refutes R2 C1: "STATE wins only because it has 2058 dims").
scGPT (~0.649) fills its `tab:ablation` `—`, tied with UCE at the bottom.

## D — auxiliary tables on the COVID+Sepsis matched cohort (cohort-consistency, R2 C2 / R1 W6)
Filtered subset embedding dir + reran multiclassifier / scIB / activation / ablation / fusion.
`results/subset/*.csv`. **scGPT now present in every subset table.**

- **Cohort-robust (subset ≈ full pool):** activation axis — scGPT ρ = **−0.675** (subset) vs −0.669
  (full), STATE +0.457 both → exhaustion/activation signs unchanged. scIB NMI-sev ordering
  unchanged (UCE > STATE > … > scGPT).
- **Cohort-sensitive (the per-fold-mean metric):** subset multiclassifier LR ranks
  scPRINT-2 0.757 > Geneformer 0.733 > TF 0.723 > **STATE 0.716** > scGPT 0.663 > UCE 0.604 —
  i.e. on the subset the *per-fold-mean* metric puts scPRINT-2 on top (full pool put STATE on top).
  Same at matched dims (subset PCA-256: scPRINT-2 0.748 > Geneformer 0.710 > STATE 0.704).
- **Reconciles cleanly with the headline:** the auxiliary tables use *per-fold mean* (seed 42),
  which we showed (seed-stability, fold-CI) is the unstable metric; the headline OOF/donor-
  clustered inference (STATE) is the robust one. Dimensionality is **not** the driver — within
  each cohort the lead is preserved under PCA-256; *which* FM leads tracks cohort × metric × seed,
  exactly as established. → Use the matched-cohort subset tables in the revision (R2 C2), and frame
  per-fold-mean rankings as the cohort/seed-sensitive view alongside the robust OOF ranking.

## Net effect on the manuscript
- **Strengthen + keep:** STATE leads (donor-FDR, 6/8 seeds, label-robust retrain).
- **Soften to ties:** scPRINT-2/Geneformer/scGPT ordering; per-fold "scPRINT-2 wins"; "STATE only
  bidirectional generalizer."
- **Revise (conclusion change):** the LoRA "+0.006 null" → within-backbone +0.113 (J).
- **Keep but mark descriptive:** two-axis biology (H null on mechanism); scGPT 0.352 (N: real).
