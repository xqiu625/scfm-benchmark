#!/usr/bin/env python3
"""Batch-effect audit, Step 2 (confounding half): patient-level batch x
outcome/disease tables from canonical_donor_labels.csv alone. v2.

Needs no capture-level join -- cv_group and well_block are already resolved,
cohort-wide, in canonical_donor_labels.csv (built and verified exact in the
Phase-1 rebuild, 2026-09-15/19; see docs/decisions.md). This script only
cross-tabulates what is already there; it does not compute a new grouping
key and does not select or imply a batch-correction setting (outcome
information is for reporting here, never for tuning a correction -- see
docs/current-status.md "Outcome-data boundary").

v2 corrects three things a 2026-10-09 review (review:20261009_000948)
found in v1, confirmed independently before fixing:

- TERMINOLOGY/INTERPRETATION: v1 called a cv_group containing both fatal and
  survived donors an "unsplittable block" / spoke of it as a demonstrated
  "batch effect confound". A cv_group is a CV DEPENDENCY COMPONENT (donors
  that must move together between folds, per the connected-components
  construction in docs/decisions.md); containing both outcome classes shows
  OUTCOME OVERLAP within a component, which constrains how CV splits can be
  drawn -- it is not, by itself, evidence of a batch-biology ASSOCIATION
  (that needs the Step 4 diagnostics on actual embeddings/expression). v2
  renames the file and columns to say exactly this, and -- per the review --
  ALSO reports the homogeneous-outcome multi-donor components (v1 silently
  dropped them from its "unsplittable" file even though they must stay
  together in CV too; this file now reports all 12 multi-donor components,
  not just the 8 with mixed known outcomes).
- DTYPE FRAGILITY: v1 classified `event` by string-matching literal "1.0" /
  "0.0", which depends on pandas having inferred a float dtype (true for the
  full file, because other rows have NaN) but is not a robust comparison in
  general (a NaN-free subset would read as int64, and "1" != "1.0"). v2
  compares numeric values directly via pd.to_numeric with explicit NaN
  handling, and asserts no unexpected non-{0,1,NaN} value slipped through.
- OVER-INCLUSION: v1's mixed-outcome filter was `nunique() > 1`, which would
  also have flagged a fatal+NA_undefined or survived+NA_undefined component
  as "unsplittable" even though the report's own text claims to be about
  fatal-AND-survived mixing specifically. (In this run it happened not to
  matter -- all 8 v1-flagged components were genuinely fatal+survived -- but
  the code did not guarantee that.) v2 requires the exact set
  {"fatal", "survived"} to both be present.

v2 also adds a Stephenson grouping-provenance check the review flagged as a
genuinely unresolved contradiction between two files this project already
treats as authoritative, rather than resolving it either way:
`stephenson_pool_composition.csv` records BGCV02 as pooling CV0052/CV0059/
CV0068/CV0902, but canonical_donor_labels.csv gives each of those four donors
its own singleton cv_group/well_block. This script reports that contradiction
explicitly; it does not merge the donors into one component and does not
certify them as independent captures.

v3 fixes defects a second 2026-10-09 review (review:20261009_101831)
found in v2, confirmed independently before fixing:

- `classify_event` used `np.isclose` for the 0/1 comparison, which would
  accept an out-of-contract value like 1.000001 as "fatal" instead of
  rejecting it. v3 uses exact equality (`== 1.0` / `== 0.0`) -- safe here
  since these are small integers with exact float64 representations, and
  the contract is exact binary labels, not a continuous score.
- Donor-uniqueness was checked on bare `donor_id` globally, which would
  reject two unrelated cohorts coincidentally reusing a donor string as two
  different biological patients. v3 validates `(cohort, donor_id)`
  uniqueness instead, and separately LOGS (does not reject) donor_id values
  that recur across different cohorts, since that is expected, not an error.
- Component/contradiction tables built from an empty list produced a
  DataFrame with no columns, which is not a loadable/stable empty schema.
  v3 passes an explicit `columns=` for every output, empty or not.

Outputs:
  cohort_x_event_table_v3.csv              -- in-scope cohorts x event (0/1/NA)
  cv_group_dependency_component_summary_v3.csv -- per cohort: singleton vs
                                              multi-donor component counts
  all_multi_donor_dependency_components_v3.csv -- EVERY multi-donor cv_group
                                              (12, not just the 8 mixed ones),
                                              with a mixed_known_outcome flag
  mixed_known_outcome_dependency_components_v3.csv -- the subset with both
                                              fatal and survived donors
  stephenson_pool_vs_cv_group_contradiction_v3.csv -- the unresolved mismatch
                                              described above (see also the
                                              companion resolve_stephenson_
                                              grouping.py, which adds the
                                              primary-source evidence this
                                              file's check does not have)
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("batch_confounding")

IN_SCOPE_COHORTS = [
    "Qiu 2021", "Cillo 2021", "Chongqing sepsis 2025", "Stephenson 2021",
    "Liu 2021", "Wilk 2021 JEM", "Ren 2021", "COMBAT 2022",
]
CANONICAL_LABELS = "reports/manuscript/revision/literature_audit/canonical_donor_labels.csv"
STEPHENSON_POOLS = (
    "reports/manuscript/revision/literature_audit/T3_patient_sample_mapping/"
    "stephenson_pool_composition.csv"
)


def classify_event(event_col: pd.Series) -> pd.Series:
    """Dtype-robust 0/1/NaN classification -- compares numeric VALUES, not
    their string repr, and fails loudly on anything unexpected rather than
    silently folding it into NA_undefined. Uses EXACT equality, not
    np.isclose -- the contract is exact binary labels (0/1/NaN), and a
    tolerance-based comparison would silently accept a value like
    1.000001 as "fatal" instead of rejecting it as unexpected (flagged in
    the 2026-10-09 review; see tests/test_batch_confounding.py for the
    rejection case)."""
    numeric = pd.to_numeric(event_col, errors="raise")
    is_one = numeric == 1.0
    is_zero = numeric == 0.0
    is_na = numeric.isna()
    unexpected = ~(is_one | is_zero | is_na)
    if unexpected.any():
        bad = numeric[unexpected].unique()
        raise ValueError(
            f"classify_event: {unexpected.sum()} value(s) are neither exactly 0, exactly 1, "
            f"nor NaN: {bad}. Refusing to silently classify these as NA_undefined."
        )
    out = pd.Series(np.where(is_one, "fatal", np.where(is_zero, "survived", "NA_undefined")),
                     index=event_col.index)
    return out


def check_stephenson_pool_vs_cv_group(repo_root: Path, scope: pd.DataFrame) -> pd.DataFrame:
    pool_path = repo_root / STEPHENSON_POOLS
    if not pool_path.exists():
        return pd.DataFrame([{"note": f"MISSING: {pool_path}"}])
    pools = pd.read_csv(pool_path)
    steph = scope[scope["cohort"] == "Stephenson 2021"].set_index("donor_id")
    rows = []
    for _, prow in pools.iterrows():
        donors = [d for d in str(prow["donors"]).split(";") if d]
        if len(donors) < 2:
            continue
        cv_groups = {d: steph.loc[d, "cv_group"] if d in steph.index else "DONOR_NOT_IN_SCOPE"
                     for d in donors}
        well_blocks = {d: steph.loc[d, "well_block"] if d in steph.index else "DONOR_NOT_IN_SCOPE"
                       for d in donors}
        n_distinct_cv_groups = len(set(cv_groups.values()))
        rows.append(dict(
            pool=prow["pool"], donors=";".join(donors), n_donors=len(donors),
            cv_groups=";".join(f"{d}={g}" for d, g in cv_groups.items()),
            well_blocks=";".join(f"{d}={b}" for d, b in well_blocks.items()),
            contradiction=(n_distinct_cv_groups > 1),
            note=(
                "pool_composition.csv records these donors as pooled together, but "
                "canonical_donor_labels.csv assigns them distinct singleton cv_groups "
                "-- UNRESOLVED, not merged and not certified independent here."
                if n_distinct_cv_groups > 1 else ""
            ),
        ))
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    canon = pd.read_csv(args.repo_root / CANONICAL_LABELS)
    scope = canon[canon["cohort"].isin(IN_SCOPE_COHORTS)].copy()
    n_excluded = len(canon) - len(scope)
    log.info("canonical_donor_labels.csv: %d total rows, %d in the 8-cohort scope, "
             "%d excluded (Schulte-Schrepping/Bernardes)", len(canon), len(scope), n_excluded)

    dup_composite = scope.duplicated(subset=["cohort", "donor_id"], keep=False)
    if dup_composite.any():
        raise ValueError(
            f"canonical_donor_labels.csv has {int(dup_composite.sum())} rows sharing a "
            f"(cohort, donor_id) key within the in-scope cohorts -- counting rows as "
            f"patients would be wrong until this is resolved:\n"
            f"{scope[dup_composite].to_string()}"
        )
    # donor_id reused ACROSS DIFFERENT cohorts is not itself an error -- two
    # unrelated studies can coincidentally use the same donor string for two
    # different biological patients (2026-10-09 review). Report it, don't raise.
    reused_across_cohorts = (
        scope.groupby("donor_id")["cohort"].nunique().loc[lambda s: s > 1]
    )
    if len(reused_across_cohorts):
        log.info("%d donor_id value(s) are reused across more than one cohort "
                  "(not an error -- different biological patients, different "
                  "cohorts): %s", len(reused_across_cohorts),
                  reused_across_cohorts.index.tolist()[:10])

    scope["event_bin"] = classify_event(scope["event"])

    ct = pd.crosstab(scope["cohort"], scope["event_bin"], margins=True, margins_name="TOTAL")
    ct_out = args.out_dir / "cohort_x_event_table_v3.csv"
    ct.to_csv(ct_out)
    log.info("Wrote %s", ct_out)

    COMPONENT_COLS = ["cohort", "cv_group", "n_donors", "donor_ids", "event_classes",
                       "mixed_known_outcome"]
    block_rows = []
    all_components = []
    mixed_components = []
    for cohort, g in scope.groupby("cohort"):
        sizes = g.groupby("cv_group").size()
        n_singleton = int((sizes == 1).sum())
        multi_groups = sizes[sizes > 1]
        block_rows.append(dict(
            cohort=cohort, n_donors=len(g), n_singleton_cv_groups=n_singleton,
            n_multi_donor_cv_groups=len(multi_groups),
        ))
        for grp, size in multi_groups.items():
            gg = g[g["cv_group"] == grp]
            classes = set(gg["event_bin"].unique())
            is_mixed_known = {"fatal", "survived"}.issubset(classes)
            row = dict(
                cohort=cohort, cv_group=grp, n_donors=int(size),
                donor_ids=";".join(gg["donor_id"].tolist()),
                event_classes=";".join(sorted(classes)),
                mixed_known_outcome=is_mixed_known,
            )
            all_components.append(row)
            if is_mixed_known:
                mixed_components.append(row)

    block_df = pd.DataFrame(block_rows, columns=["cohort", "n_donors", "n_singleton_cv_groups",
                                                   "n_multi_donor_cv_groups"])
    block_out = args.out_dir / "cv_group_dependency_component_summary_v3.csv"
    block_df.to_csv(block_out, index=False)
    log.info("Wrote %s", block_out)

    all_df = pd.DataFrame(all_components, columns=COMPONENT_COLS)
    all_out = args.out_dir / "all_multi_donor_dependency_components_v3.csv"
    all_df.to_csv(all_out, index=False)
    log.info("Wrote %s (%d multi-donor components total, %d mixed-known-outcome)",
             all_out, len(all_df), len(mixed_components))

    mixed_df = pd.DataFrame(mixed_components, columns=COMPONENT_COLS)
    mixed_out = args.out_dir / "mixed_known_outcome_dependency_components_v3.csv"
    mixed_df.to_csv(mixed_out, index=False)
    log.info("Wrote %s", mixed_out)

    STEPH_COLS = ["pool", "donors", "n_donors", "cv_groups", "well_blocks",
                  "contradiction", "note"]
    steph_df = check_stephenson_pool_vs_cv_group(args.repo_root, scope)
    if steph_df.empty and list(steph_df.columns) == []:
        steph_df = pd.DataFrame(columns=STEPH_COLS)
    steph_out = args.out_dir / "stephenson_pool_vs_cv_group_contradiction_v3.csv"
    steph_df.to_csv(steph_out, index=False)
    n_contra = int(steph_df.get("contradiction", pd.Series(dtype=bool)).sum()) if len(steph_df) else 0
    log.info("Wrote %s (%d pools with a cv_group contradiction)", steph_out, n_contra)

    print("\n=== COHORT x EVENT ===")
    print(ct.to_string())
    print("\n=== DEPENDENCY COMPONENT SUMMARY (all cv_groups, by cohort) ===")
    print(block_df.to_string(index=False))
    print("\n=== ALL MULTI-DONOR DEPENDENCY COMPONENTS (%d total, %d mixed-known-outcome) ==="
          % (len(all_df), len(mixed_components)))
    print(all_df.to_string(index=False) if len(all_df) else "(none)")
    print("\n=== STEPHENSON POOL vs cv_group CONTRADICTION CHECK ===")
    print(steph_df.to_string(index=False) if len(steph_df) else "(none)")


if __name__ == "__main__":
    main()
