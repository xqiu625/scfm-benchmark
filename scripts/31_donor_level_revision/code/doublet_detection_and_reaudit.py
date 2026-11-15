#!/usr/bin/env python3
"""
Targeted response to the contamination-audit finding (consensus_gene_
contamination_audit.py): run doublet detection on the canonical platelet
cohort and re-check whether excluding doublets changes the "likely
technical artifact" signature found for CD79A/CD79B/PRF1/SLC25A37/HBA1/
CD3E/GZMA/CCL5 (strong library-depth correlation, weak-to-null direct
mortality correlation).

Scope decision (2026-09-13, user-confirmed): TARGETED re-check only --
doublet-detect + re-audit the flagged genes and the R6 consensus result.
NOT a full pipeline rerun: FM embeddings are not recomputed (doublet
exclusion is applied only to this diagnostic, not threaded back through
the 5 models' frozen embeddings or the primary Table 3 results), since
canonical ambient-RNA decontamination (SoupX/CellBender) was found
infeasible on this data (no raw unfiltered/empty-droplet matrices
available for these already-cell-called public source datasets) and a
full benchmark rerun is an explicitly out-of-scope, multi-week
undertaking the user did not choose.

Method: sc.pp.scrublet with batch_key="Sample_ID" (native per-sample/
per-capture-run doublet simulation -- doublets are a per-capture
technical artifact, so must not be estimated on the pooled multi-study
object). Platelet-specific caveat: these are anucleate, very-low-RNA-
content fragments; Scrublet's doublet-simulation approach is calibrated
on typical nucleated-cell libraries, so results are reported with that
caveat rather than treated as a definitive doublet call.

SEMANTICS FIX 2026-09-22 (P0b). The previous version wrote
`predicted_doublet = False` for every cell in its fallback branch and left it
False for samples with <30 cells and for samples whose Scrublet call raised.
Those cells were NEVER EXAMINED, but False asserts "examined and not a
doublet". The consequences ran through four layers: the per-cell export, the
per-patient rate (unexamined samples entered as a measured 0%), the exclusion
mask (`.fillna(False)`), and the Spearman correlations computed on those rates.

The state is now a nullable boolean `doublet_call` (True / False / <NA>) plus an
explicit `examination_status`. Rules, enforced by the pure functions below and
covered by tests/test_stage31_doublet_semantics.py:

  * <30 cells, runtime failures and never-run samples stay <NA>, never False;
  * only `doublet_call == True` is removed; unknown cells are retained AND stay
    distinguishable in a_clean.obs -- the result is not a "doublet-cleaned" set;
  * per-patient rate uses EXAMINED cells as the denominator, and is NA (not 0%)
    for a patient with no examined cells;
  * the Spearman analyses exclude undefined rates explicitly and report the
    number of contributing patients.

Outputs:
  doublet_scores_per_cell.csv       per-cell doublet_score, doublet_call, examination_status
  contamination_audit_v2_doublet_excluded.csv
                                     same correlations as v1, recomputed
                                     after excluding predicted doublets
  doublet_rate_per_patient.csv      sanity check: does doublet rate itself
                                     correlate with depth/mortality
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

#: Scrublet needs a minimum batch size; below it a sample is NOT EVALUATED.
MIN_CELLS_PER_SAMPLE = 30

#: PCA dimensionality handed to Scrublet. There is NO default.
#:
#: The value is a methodological choice that has not been made
#: (docs/scrublet_dimension_parameter.md is OPEN), so the script must not be
#: runnable without one. A provisional constant recorded in the manifest does
#: not solve this: it still lets anyone produce results under an unapproved
#: parameter, and "provisional" in a comment binds nothing. `--n-prin-comps` is
#: therefore REQUIRED, and Scrublet cannot be reached without it.
#:
#: It was previously not passed at all -- both call sites took scanpy's
#: implicit default -- so the run had a dimensionality nobody had chosen.
N_PRIN_COMPS_UNSET = None


def scrublet_params(n_prin_comps):
    """The parameters bound to EVERY sc.pp.scrublet call. Fails closed.

    A parameter that is not in this dict is not bound, whatever the manifest
    says about it.
    """
    if n_prin_comps is None:
        raise SystemExit(
            "--n-prin-comps is required and has no default. The dimensionality is a "
            "methodological choice; see docs/scrublet_dimension_parameter.md, which is "
            "still OPEN. Pass the value the approved command specifies.")
    if isinstance(n_prin_comps, bool) or not isinstance(n_prin_comps, (int, np.integer)):
        raise SystemExit(f"--n-prin-comps must be an integer, got {n_prin_comps!r}")
    if int(n_prin_comps) < 1:
        raise SystemExit(f"--n-prin-comps must be >= 1, got {n_prin_comps}")
    return {"n_prin_comps": int(n_prin_comps)}

#: Why a cell does or does not carry a call. `examined` is the ONLY status for
#: which doublet_call is a claim about the cell; the rest mean "we do not know".
EXAM_EXAMINED = "examined"
EXAM_TOO_SMALL = "too_small"
EXAM_THRESHOLD_FAILED = "threshold_failed"
EXAM_RUNTIME_FAILED = "runtime_failed"
#: The sample cannot support N_PRIN_COMPS: PCA needs
#: n_prin_comps < min(n_cells, usable_features). This is a FEASIBILITY fact
#: known before Scrublet is called, not a crash, and it is kept separate from
#: runtime_failed so "the detector broke" and "we never had the dimensions to
#: ask" stay distinguishable in the audit.
EXAM_INSUFFICIENT_DIMENSION = "insufficient_dimension"
#: The preprocessing that measures the sample's available dimensions raised, so
#: feasibility could not even be evaluated. Distinct from insufficient_dimension
#: ("we measured, and it is not enough") and from runtime_failed ("the detector
#: itself broke"): here nothing was ever asked of the detector.
EXAM_DIMENSION_MEASUREMENT_FAILED = "dimension_measurement_failed"
EXAM_NOT_EVALUATED = "not_evaluated"
EXAM_UNKNOWN_STATUSES = (EXAM_TOO_SMALL, EXAM_THRESHOLD_FAILED,
                         EXAM_RUNTIME_FAILED, EXAM_INSUFFICIENT_DIMENSION,
                         EXAM_DIMENSION_MEASUREMENT_FAILED, EXAM_NOT_EVALUATED)


class ExaminationStateError(RuntimeError):
    """status and doublet_call contradict each other. Always fatal."""


def strict_boolean_calls(values, n_expected):
    """Return a nullable boolean array, or None with a reason.

    Deliberately NOT a bare numpy `.astype(bool)`, which maps every non-zero
    value to True and would invent calls out of whatever the column held.
    (pandas' own `dtype="boolean"` does reject floats, but it still accepts a
    wrong-length or NA-bearing array without comment.) Only a real boolean
    array, or an integer array strictly in {0, 1}, is a call.
    A single NA anywhere disqualifies the WHOLE sample -- a per-cell NA means the
    detector did not produce a usable result for that sample, and keeping the
    rest would leave `examined` cells with no call.
    """
    try:
        arr = np.asarray(values)
    except Exception as exc:
        return None, f"unreadable ({type(exc).__name__})"
    if arr.ndim != 1:
        return None, f"not 1-D (ndim={arr.ndim})"
    if arr.shape[0] != n_expected:
        return None, f"length {arr.shape[0]} != {n_expected} cells"
    if arr.dtype == object and pd.isna(arr).any():
        return None, "contains NA"
    if arr.dtype == bool:
        return pd.array(arr, dtype="boolean"), None
    if arr.dtype.kind in "iu" and bool(np.all((arr == 0) | (arr == 1))):
        return pd.array(arr.astype(bool), dtype="boolean"), None
    return None, f"not a strict boolean or 0/1 integer (dtype={arr.dtype})"


def scrublet_preprocessing_dimensions(counts_adata):
    """What PCA will actually see, computed by scanpy's OWN scrublet pipeline.

    Not an estimate. scanpy 1.11.3 `pp.scrublet` preprocesses each batch as:

        filter_genes(min_cells=3)  ->  filter_cells(min_genes=3)
        normalize_total  ->  log1p (temp layer)  ->  highly_variable_genes
        subset to highly_variable  ->  ...  ->  PCA(n_components=n_prin_comps)

    Two consequences the raw cell count hides, and both bit this run:

      * `filter_cells(min_genes=3)` changes the sample size, so the n_samples
        PCA sees is NOT the number of cells in the sample;
      * n_features is the HVG count, not the 13,066 genes in the matrix.

    The v1 batched failure reported `min(n_samples, n_features)=17`. That is a
    MINIMUM and does not say which term was 17, or whether both were -- which
    is why both are measured here rather than read off the error message.

    Returns (n_cells_effective, n_usable_features) -- (0, 0) if filtering
    empties the sample.
    """
    ad = counts_adata.copy()
    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.filter_cells(ad, min_genes=3)
    if ad.n_obs == 0 or ad.n_vars == 0:
        return 0, 0
    sc.pp.normalize_total(ad)
    ad.layers["log1p"] = ad.X.copy()
    sc.pp.log1p(ad, layer="log1p")
    sc.pp.highly_variable_genes(ad, layer="log1p")
    return int(ad.n_obs), int(ad.var["highly_variable"].sum())


#: Preflight exit codes. A non-zero code is the only signal SLURM reads: a run
#: that supports no dimensionality decision must not be recorded as COMPLETED
#: while its own manifest says `incomplete`. Fixing the manifest alone left
#: those two signals contradicting each other.
PREFLIGHT_EXIT_OK = 0
PREFLIGHT_EXIT_INSUFFICIENT = 3


def preflight_exit_code(sufficient_all_samples, sufficient_eligible_only):
    """Exit status for a dimension-audit run.

    Zero when at least one distribution is complete enough to support a
    decision. An incomplete all-samples distribution beside a complete
    eligible-only one is still a usable result -- it supports an explicitly
    narrowed decision -- so it exits 0, and the methodological conclusion must
    then record that it used the eligible-only distribution.
    """
    if sufficient_all_samples or sufficient_eligible_only:
        return PREFLIGHT_EXIT_OK
    return PREFLIGHT_EXIT_INSUFFICIENT


def dimension_feasibility(n_cells_effective, n_usable_features, n_prin_comps):
    """Can this sample support n_prin_comps? Returns (ok, reason_or_None).

    sklearn's arpack solver requires n_components STRICTLY less than
    min(n_samples, n_features); equality raises. Checked ahead of the call so
    the outcome is a recorded feasibility decision rather than a caught
    exception, and so it can never be relabelled as a runtime failure.

    There is deliberately no fallback to a smaller n_prin_comps: a per-sample
    dimensionality would give each sample its own detector working point, and
    the per-patient rate aggregates across samples.
    """
    limit = min(n_cells_effective, n_usable_features)
    if n_prin_comps < limit:
        return True, None
    return False, (f"n_prin_comps={n_prin_comps} needs "
                   f"min(n_cells_effective={n_cells_effective}, "
                   f"usable_features={n_usable_features})={limit} to be strictly greater")


def validate_examination_state(obs, *, where="obs"):
    """Enforce the status <-> call invariant. Raises; never repairs silently.

    Closing `fillna(False)` removed the explicit collapse but left an implicit
    one: a cell marked `examined` whose call is NA enters the rate denominator
    while contributing nothing to the numerator, which dilutes unknown into
    non-doublet exactly as before. The two fields must agree.

        examined        <=> doublet_call in {True, False}
        unknown status  <=> doublet_call is NA
    """
    status = obs["examination_status"]
    call = obs["doublet_call"]

    bad_status = sorted(set(status) - {EXAM_EXAMINED, *EXAM_UNKNOWN_STATUSES})
    if bad_status:
        raise ExaminationStateError(f"{where}: unknown examination_status {bad_status}")

    # The column itself must be the nullable boolean dtype. Checking only
    # NA/non-NA is not enough: an object column holding "False", 0.0 or 2 is
    # non-NA and would sail through while meaning nothing.
    if not pd.api.types.is_bool_dtype(call.dtype):
        raise ExaminationStateError(
            f"{where}: doublet_call has dtype {call.dtype}, expected a nullable boolean. "
            "A non-boolean column cannot carry a call -- 'False', 0.0 and 2 are all non-NA "
            "and would pass an NA-only check while meaning nothing.")

    examined = (status == EXAM_EXAMINED)
    n_examined_na = int(call[examined].isna().sum())
    if n_examined_na:
        raise ExaminationStateError(
            f"{where}: {n_examined_na} cells are marked '{EXAM_EXAMINED}' but carry a NA "
            "doublet_call. They would enter the rate denominator without contributing to "
            "the numerator, diluting unknown into non-doublet.")

    n_unknown_called = int(call[~examined].notna().sum())
    if n_unknown_called:
        raise ExaminationStateError(
            f"{where}: {n_unknown_called} cells carry a doublet_call while their status is "
            f"not '{EXAM_EXAMINED}'. A call is only meaningful for an examined cell.")
    return obs


def summarize_examination(obs):
    """Counts by examination status. `n_doublet` counts confirmed doublets only;
    unknown cells are never folded into the non-doublet total.

    Validates the state first: a contradictory input would make these counts
    quietly wrong rather than loudly absent.
    """
    validate_examination_state(obs, where="summarize_examination")
    call = obs["doublet_call"]
    examined = (obs["examination_status"] == EXAM_EXAMINED)
    out = {
        "n_cells": int(len(obs)),
        "n_examined": int(examined.sum()),
        "n_unknown": int((~examined).sum()),
        "n_doublet": int((call == True).sum()),          # noqa: E712 -- nullable
        "n_non_doublet": int((call == False).sum()),     # noqa: E712 -- nullable
        "coverage": float(examined.sum() / len(obs)) if len(obs) else float("nan"),
        "by_status": dict(obs["examination_status"].value_counts()),
    }
    assert out["n_doublet"] + out["n_non_doublet"] == out["n_examined"], (
        f"examined cells must be exactly doublet + non-doublet, got "
        f"{out['n_doublet']} + {out['n_non_doublet']} != {out['n_examined']}")
    return out


def doublet_rate_per_patient(obs):
    """Per-patient doublet rate over EXAMINED cells only.

    A patient with no examined cells gets NaN, never 0.0 -- 0% would read as a
    measured absence of doublets rather than an absent measurement. The examined
    and unknown counts travel with the rate so a reader can see the denominator.
    """
    validate_examination_state(obs, where="doublet_rate_per_patient")
    examined = (obs["examination_status"] == EXAM_EXAMINED)
    rows = []
    for pid, sub in obs.groupby("patient_id", observed=True):
        ex = sub[examined.loc[sub.index]]
        n_ex = int(len(ex))
        rate = float((ex["doublet_call"] == True).sum() / n_ex) if n_ex else float("nan")  # noqa: E712
        rows.append({"patient_id": pid, "doublet_rate": rate,
                     "n_examined": n_ex, "n_unknown": int(len(sub) - n_ex),
                     "n_cells": int(len(sub)),
                     "examined_coverage": float(n_ex / len(sub)) if len(sub) else float("nan")})
    return pd.DataFrame(rows).set_index("patient_id")

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")
# Output directory is REQUIRED, must be new, and must not be the legacy one.
# The old default pointed at the directory holding the four artifacts produced
# under the broken semantics; running this script would have overwritten three
# of them in place, destroying the record of what the defect produced. They are
# registered stale in scripts/31_donor_level_revision/STALE_ARTIFACTS.md.
LEGACY_OUT_DIR = ("/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/"
                  "contamination_audit")

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out-dir", required=True,
                 help="NEW versioned output directory: must not exist or must be empty, "
                      "and must not resolve to the legacy directory")
_ap.add_argument("--integration-only", action="store_true",
                 help="run detection, export and per-patient rates, but SKIP every "
                      "mortality/outcome correlation. Use this to validate the Scanpy, "
                      "AnnData and CSV semantics without executing an outcome analysis. "
                      "Boundary: outcome is excluded at PARSE time from the patient list "
                      "(usecols); the raw h5ad container DOES carry mortality and severity "
                      "columns into memory because read_h5ad cannot select columns, and "
                      "only the columns the loader requires survive -- the rest are dropped "
                      "before use and recorded in the manifest. So the manifest reports "
                      "known_outcome_values_loaded=true with outcome_values_inspected=false. "
                      "This "
                      "is NOT a claim that outcome is never loaded.")
_ap.add_argument("--dimension-audit", action="store_true",
                 help="PREFLIGHT. Measure per-sample n_cells_effective and "
                      "n_usable_features and stop. Does not require --n-prin-comps, does "
                      "not call Scrublet, and produces no call, rate or contamination "
                      "output. This exists because requiring an approved dimensionality "
                      "before running, while saying the decision should rest on this "
                      "run's diagnostics, is circular: the distribution has to be "
                      "obtainable without first picking a value arbitrarily. Every "
                      "sample is measured, including those below MIN_CELLS_PER_SAMPLE.")
_ap.add_argument("--n-prin-comps", type=int, default=None,
                 help="PCA dimensionality for Scrublet. REQUIRED -- there is no default, "
                      "because the value is a methodological choice that has not been "
                      "made (docs/scrublet_dimension_parameter.md is OPEN). It is used by "
                      "the batched path and the per-sample fallback alike and is never "
                      "reduced on failure; a sample that cannot support it is marked "
                      "insufficient_dimension rather than scored at some other dimension.")
_args = _ap.parse_args()
OUT_DIR = _args.out_dir
DIMENSION_AUDIT = _args.dimension_audit

# The two modes are exclusive. The reason is GOVERNANCE, not measurement:
# n_prin_comps does not enter the preprocessing at all, so supplying one would
# not change a single measured number. It would, however, put an unapproved
# value into the command line and the manifest of the run whose purpose is to
# inform that choice, and so read as though the value were already settled or
# applied. An earlier comment here claimed the measurement would be taken
# "under" the value; that was wrong.
if DIMENSION_AUDIT and _args.n_prin_comps is not None:
    raise SystemExit("--dimension-audit takes no --n-prin-comps. The value does not "
                     "affect what is measured; passing one would record an unapproved "
                     "value in the run that exists to inform the choice.")
if DIMENSION_AUDIT and not _args.integration_only:
    raise SystemExit("--dimension-audit must be run with --integration-only: it is a "
                     "technical preflight and has no business loading outcome.")

#: Bound once, here, from the command that was actually run. None in preflight,
#: where Scrublet is never reached.
N_PRIN_COMPS = _args.n_prin_comps
SCRUBLET_PARAMS = None if DIMENSION_AUDIT else scrublet_params(N_PRIN_COMPS)
INTEGRATION_ONLY = _args.integration_only

# Resolve both sides: a symlink, a relative path or a trailing slash must not be
# able to land on the legacy directory by accident.
if os.path.realpath(OUT_DIR) == os.path.realpath(LEGACY_OUT_DIR):
    raise SystemExit(
        f"--out-dir resolves to the legacy directory {LEGACY_OUT_DIR}. Its artifacts are "
        "registered stale and must survive untouched; write to a new versioned directory.")
if os.path.exists(OUT_DIR):
    if not os.path.isdir(OUT_DIR):
        raise SystemExit(f"--out-dir {OUT_DIR} exists and is not a directory")
    existing = os.listdir(OUT_DIR)
    if existing:
        raise SystemExit(
            f"--out-dir {OUT_DIR} is not empty ({existing[:5]}). A rerun writes to a NEW "
            "versioned directory so a previous run's outputs are never mixed with this "
            "one's or silently replaced.")
os.makedirs(OUT_DIR, exist_ok=True)


def write_csv_atomic(df, path, **kwargs):
    """Write to a temp in the same directory, then rename.

    A job killed mid-write would otherwise leave a truncated CSV that looks
    complete -- the failure mode the stale artifacts already demonstrate is
    expensive.
    """
    tmp = path + ".inprogress"
    try:
        df.to_csv(tmp, **kwargs)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


FLAGGED_GENES = {
    "CD79B": "B_cell", "CD79A": "B_cell",
    "PRF1": "NK_cell",
    "S100A8": "Monocyte_Myeloid", "S100A12": "Monocyte_Myeloid",
    "FCER1G": "Monocyte_Myeloid", "TYROBP": "Monocyte_Myeloid",
    "SLC25A37": "Erythroid", "HBA1": "Erythroid",
    "CD3E": "T_cell", "GZMA": "T_cell", "CCL5": "T_cell",
    "GP1BA": "Platelet_reference_positive_control",
}

CELLS_PER_PATIENT = 60
RNG_SEED = 20260910


# ---------------------------------------------------------------------
# THE OUTCOME BOUNDARY under --integration-only.
#
# Two inputs, and they can offer DIFFERENT guarantees. Declaring one guarantee
# for both would be false, so they are declared separately.
#
#   patient list CSV -- STRICT: `y` is excluded by `usecols` at PARSE time, so
#       no outcome value is ever constructed in memory. A previous version read
#       the whole file and then called .drop(columns=["y"]); that proves only
#       that later code does not use y, not that y was never read.
#
#   raw h5ad -- INCIDENTAL LOAD, ACKNOWLEDGED: read_h5ad offers no column
#       selection, so obs arrives whole. This file's obs carries `Death`, and
#       severity labels besides. The container therefore brings clinical
#       outcome into memory unavoidably. Only REQUIRED_OBS_COLUMNS survives the
#       loader; everything else is dropped before any other code touches obs,
#       and the manifest records what was dropped rather than hiding it.
#
# So the claim this script is entitled to make is NOT "outcome is not read". It
# is: outcome is excluded at parse time where the reader allows it, and where
# the reader does not, it is loaded incidentally as part of the input container
# and is never inspected, propagated, exported or analysed. The strict
# guarantee on the h5ad would require a separate outcome-free technical input,
# which does not exist and is not created by this run.
#
# The previous code dropped the wrong column here: it dropped "y", which this
# h5ad does not have, while `Death` -- the actual mortality variable -- sat in
# obs untouched. It neither achieved the boundary it claimed nor did the thing
# it said it was doing.
#
# RETENTION, NOT DELETION (2026-09-22). A deny-list of outcome columns was the
# first fix and it was wrong in the same way: `("y", "Death")` left
# `covid_severity` -- and `mild_vs_severe`, `severity_score`, `three_group`,
# `progression_group` -- propagating into the AnnData, the slice and the
# temporary h5ad round trip, while the code and docs said "outcome" throughout.
# A deny-list is only as good as the classification behind it, and a new
# clinical label added upstream would propagate by default.
#
# So the direction is inverted: the loader needs exactly ONE column from the
# container (`Sample_ID`, to map cells to patients -- everything else it uses,
# it creates). Keeping only that makes the guarantee independent of whether the
# classification below is complete: an unclassified clinical label cannot leak,
# because nothing but Sample_ID survives. The classification is then used only
# to REPORT what was dropped, and columns it cannot adjudicate are recorded as
# unclassified rather than quietly filed as technical.
PATIENT_LIST_BASE_COLS = ["patient_id", "disease",
                          "median_counts", "median_genes", "median_pct_mito"]

#: The only obs columns the loader needs from the container. A missing one is
#: fatal; everything else is dropped whether or not it is recognised below.
REQUIRED_OBS_COLUMNS = ("Sample_ID",)

#: Obs columns known to carry a clinical outcome or severity label, for the
#: manifest's record of what was dropped. Named from the observed obs header of
#: human_platelet_covid_severity_raw_counts.h5ad. This list does NOT gate the
#: drop -- retention does -- so an omission here understates the report, it
#: cannot cause a leak.
KNOWN_OUTCOME_OBS_COLUMNS = (
    "y", "Death",                                    # mortality
    "covid_severity", "mild_vs_severe", "severity_score",
    "three_group", "progression_group",              # severity / progression
)


def load_canonical_cohort():
    """Load the COVID cohort under the outcome boundary declared above."""
    usecols = list(PATIENT_LIST_BASE_COLS)
    if not INTEGRATION_ONLY:
        usecols.append("y")
    # pandas raises if any name in usecols is absent, so a renamed column fails
    # the run rather than silently producing a frame without it.
    patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str},
                               usecols=usecols)
    patient_list = patient_list[patient_list["disease"] == "COVID-19"]
    assert len(patient_list) == 240, f"expected 240 COVID patients, got {len(patient_list)}"

    if INTEGRATION_ONLY:
        assert "y" not in patient_list.columns, \
            "outcome column was parsed despite --integration-only"
        y_by_patient = None
        print("  outcome column excluded at parse time (usecols); never loaded", flush=True)
    else:
        mark_known_outcome_values_loaded("patient_list_csv:y")
        mark_outcome_values_inspected("patient_list_csv:y")
        assert int(patient_list["y"].sum()) == 25, "event count changed; check the patient list"
        y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()

    qc_by_patient = patient_list.set_index("patient_id")[
        ["median_counts", "median_genes", "median_pct_mito"]]
    covid_patient_ids = set(patient_list["patient_id"])

    sample_map = pd.read_csv(SAMPLE_MAP_CSV, dtype={"sample_id": str, "patient_id": str})
    sample_map = sample_map[sample_map["patient_id"].isin(covid_patient_ids)]
    sample_to_patient = sample_map.set_index("sample_id")["patient_id"].to_dict()

    print("Loading raw counts...", flush=True)
    a = sc.read_h5ad(RAW_H5AD)
    # Incidental load, then RETENTION: keep only what the loader needs, drop
    # everything else before any other code can touch obs. This is the step
    # that makes the boundary independent of the classification -- an
    # unrecognised clinical label is dropped because it is not required, not
    # because it was successfully identified as clinical.
    arrived = list(a.obs.columns)
    missing = [c for c in REQUIRED_OBS_COLUMNS if c not in arrived]
    if missing:
        raise SystemExit(f"raw h5ad obs is missing required column(s) {missing}; "
                         f"got {arrived}")

    # ORDER MATTERS, and the two records go on opposite sides of the drop.
    #
    # Loading is already true: read_h5ad has returned and the values are in
    # memory. Record it now -- a job that dies here did load them.
    known_outcome_arrived = [c for c in arrived if c in KNOWN_OUTCOME_OBS_COLUMNS]
    for col in known_outcome_arrived:
        mark_known_outcome_values_loaded(f"raw_h5ad_obs:{col}")

    # The drop has NOT happened yet, so nothing may claim it has. An earlier
    # version wrote known_outcome_dropped before this line: a job dying in
    # between left a manifest asserting columns were gone while they were still
    # on the object -- the same predicted-fact defect as outcome_data_accessed.
    dropped = [c for c in arrived if c not in REQUIRED_OBS_COLUMNS]
    known_outcome_dropped = [c for c in dropped if c in KNOWN_OUTCOME_OBS_COLUMNS]
    unclassified_dropped = [c for c in dropped if c not in KNOWN_OUTCOME_OBS_COLUMNS]
    a.obs = a.obs[list(REQUIRED_OBS_COLUMNS)].copy()
    if set(a.obs.columns) != set(REQUIRED_OBS_COLUMNS):
        raise SystemExit(f"allowlist did not apply: obs still holds "
                         f"{sorted(set(a.obs.columns) - set(REQUIRED_OBS_COLUMNS))}")
    # `kept` is what is OBSERVED on the object now, not the constant that was
    # asked for, so the record cannot be right while the object is wrong.
    record_container_obs_handling(arrived=arrived,
                                  kept=list(a.obs.columns),
                                  known_outcome_dropped=known_outcome_dropped,
                                  unclassified_dropped=unclassified_dropped)
    print(f"  h5ad obs: kept {list(a.obs.columns)}; dropped "
          f"{len(known_outcome_dropped)} known outcome/severity column(s) "
          f"{known_outcome_dropped} and {len(unclassified_dropped)} other(s)", flush=True)
    a.obs["patient_id"] = a.obs["Sample_ID"].astype(str).map(sample_to_patient)
    a = a[a.obs["patient_id"].notna()].copy()
    if y_by_patient is not None:
        # a deliberate, already-recorded access -- not the container's copy
        a.obs["y"] = a.obs["patient_id"].map(y_by_patient).astype(int)
    print(f"  {a.n_obs} cells, {a.obs['patient_id'].nunique()} patients", flush=True)

    rng = np.random.default_rng(RNG_SEED)
    keep_idx = []
    for pid, grp in a.obs.groupby("patient_id"):
        take = min(len(grp), CELLS_PER_PATIENT)
        keep_idx.extend(rng.choice(grp.index.to_numpy(), size=take, replace=False).tolist())
    a = a[keep_idx].copy()
    print(f"  after stratified subsample: {a.n_obs} cells", flush=True)
    return a, y_by_patient, qc_by_patient


# Manifest lifecycle. Written before any computation so a killed job still says
# what it was, then atomically updated at the end. It records what was
# AUTHORIZED up front and what actually RAN only once it has run: an earlier
# version stamped outcome_analysis_executed=true before the analysis, which
# would have been a false record for any job that died first.
OUTCOME_ANALYSIS_RAN = []
MANIFEST_PATH = os.path.join(OUT_DIR, "RUN_MANIFEST.json")
import json as _json
import subprocess as _sp
import datetime as _dt

_commit = "unknown"
try:
    _commit = _sp.run(["git", "rev-parse", "HEAD"],
                      cwd=os.path.dirname(os.path.abspath(__file__)),
                      capture_output=True, text=True, timeout=10).stdout.strip() or "unknown"
except Exception:
    pass

def initial_manifest():
    """The manifest's starting values, as a function so they are testable.

    Built by a function rather than written inline because a test that supplies
    its own starting dict cannot detect a wrong starting value -- which is how a
    predicted `outcome_data_accessed` survived a behavioural test suite.
    """
    return {
        "run_status": "started",
        "started": _dt.datetime.now().isoformat(timespec="seconds"),
        "finished": None,
        "git_commit": _commit,
        "script": os.path.abspath(__file__),
        "out_dir": os.path.realpath(OUT_DIR),
        "integration_only": INTEGRATION_ONLY,
        "mode": "dimension_audit" if DIMENSION_AUDIT else "detection",
        # Two facts, not one. The previous single `scrublet_executed` was set
        # BEFORE the call, so a detector that raised at its entry point still
        # reported having executed -- it recorded an attempt. Attempt is
        # recorded before the call (a crash did attempt it) and completion only
        # after a successful return. The preflight reaches neither.
        "scrublet_attempted": False,
        # "any call", deliberately. In the fallback the first sample that
        # returns sets it, while later samples may still fail, so it means
        # "at least one Scrublet call returned" -- NOT "the detection stage
        # finished". Stage completion is run_status; per-sample outcomes are in
        # per_sample_diagnostics.csv.
        "scrublet_any_call_completed": False,
        # AUTHORIZATION -- known up front, and the only field that may be predicted.
        "outcome_analysis_requested": not INTEGRATION_ONLY,
        # FACTS -- each starts at its "has not happened" value and is flipped only
        # by the code that actually does it, with an immediate atomic rewrite. A
        # previous version wrote outcome_data_accessed = not INTEGRATION_ONLY up
        # front, so a full-mode job that died before the read claimed an access
        # that never occurred.
        #
        # LOADED vs INSPECTED are separate facts, because one field could not
        # express the real state without contradicting itself: an
        # integration-only run recorded `outcome_data_accessed: false` next to
        # `outcome_incidentally_loaded: ["Death"]`, and "loaded but never looked
        # at" is not "not accessed". Under --integration-only the honest record
        # is loaded=true (the h5ad container carried Death in), inspected=false,
        # analyses=[].
        # "known_" is load-bearing: this speaks for KNOWN_OUTCOME_OBS_COLUMNS
        # only. A container carrying a new, unrecognised clinical label leaves
        # it false while container_obs_columns.unclassified_dropped names that
        # label. Retention, not this flag, is what prevents a leak.
        "known_outcome_values_loaded": False,
        "known_outcome_values_loaded_via": [],
        "outcome_values_inspected": False,
        "outcome_values_inspected_via": [],
        "outcome_analyses_executed": [],
        # what the container brought and what survived the loader
        "container_obs_columns": None,
        "inputs": {"raw_h5ad": RAW_H5AD, "patient_list": PATIENT_LIST_CSV,
                   "sample_map": SAMPLE_MAP_CSV},
        "params": {"RNG_SEED": RNG_SEED, "CELLS_PER_PATIENT": CELLS_PER_PATIENT,
                   "MIN_CELLS_PER_SAMPLE": MIN_CELLS_PER_SAMPLE,
                   "N_PRIN_COMPS": N_PRIN_COMPS,
                   "scrublet_params": dict(SCRUBLET_PARAMS) if SCRUBLET_PARAMS else None},
    }


_MANIFEST = initial_manifest()


def write_manifest():
    """Same-directory temp, then rename -- the manifest gets the same
    durability as the CSVs rather than being the one file a kill can truncate.

    Stated as narrowly as the implementation supports:

      After `os.replace()` returns, a subsequent normal or non-zero process
      exit does not remove the replaced manifest, and no reader sees a
      partially written one.

    It does NOT guarantee anything about an interruption BEFORE the replace. An
    interruption there leaves the previous final manifest in place IF ONE
    EXISTS; during the first write there may be no final manifest at all, since
    the output guard requires a new empty directory -- so the failure mode of a
    first write is "no final manifest", not "the previous one". An uncatchable
    interruption (SIGKILL, node failure) may also leave a `.inprogress` file
    behind, because the cleanup below cannot run. And it does not fsync, so it
    offers no durability GUARANTEE across power loss, node failure or storage
    failure -- absence of a guarantee, not a guarantee of loss; without fsync the
    data may still survive through OS writeback or the filesystem's own
    behaviour, which is exactly why it cannot be relied on either way. An
    earlier version of this docstring said no "signal, non-zero status or
    crash" could lose an already-written manifest, which was broader than the
    code: "already-written" was not pinned to "os.replace has returned".

    `write_csv_atomic` gives exactly the same guarantee and the same
    non-guarantees, which is what makes "the same durability as the CSVs" the
    honest comparison.
    """
    tmp = MANIFEST_PATH + ".inprogress"
    try:
        with open(tmp, "w") as fh:
            _json.dump(_MANIFEST, fh, indent=2)
        os.replace(tmp, MANIFEST_PATH)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _mark(flag, via_list, where):
    if where not in _MANIFEST[via_list]:
        _MANIFEST[via_list].append(where)
    _MANIFEST[flag] = True
    write_manifest()


def mark_known_outcome_values_loaded(where):
    """Values of a column on KNOWN_OUTCOME_OBS_COLUMNS are in memory.

    True even when they arrived incidentally inside an input container and
    nothing ever looks at them. "known" is in the name deliberately: this flag
    speaks only for the classified list, so a container carrying nothing but a
    new, unrecognised clinical label would leave it false while
    `container_obs_columns.unclassified_dropped` names that label. The
    anti-leak guarantee does not rest on it -- retention does -- but the field
    must not be read as "no outcome value was ever in memory".
    """
    _mark("known_outcome_values_loaded", "known_outcome_values_loaded_via", where)


def mark_outcome_values_inspected(where):
    """Outcome values were actually read -- summed, mapped, correlated, printed.
    This is the fact an integration-only authorization constrains."""
    _mark("outcome_values_inspected", "outcome_values_inspected_via", where)


def record_container_obs_handling(*, arrived, kept, known_outcome_dropped,
                                  unclassified_dropped):
    """Record what the container's obs brought and what survived the allowlist.

    Called only AFTER the drop has been applied and verified: this is a record
    of an accomplished fact, not of an intention. `kept` is what was observed
    on the object, not the constant that was requested.

    Recorded on both paths including the empty case, so `null` means the load
    was never reached and a populated record means it was. Anything the
    classification cannot adjudicate is listed as unclassified rather than
    filed as technical -- it is dropped either way, because retention, not
    classification, is what enforces the boundary.
    """
    _MANIFEST["container_obs_columns"] = {
        "source": RAW_H5AD,
        "reader": "anndata.read_h5ad (no column selection available)",
        "policy": "retain only the columns the loader requires; drop all others",
        "classification_scope": (
            "known_outcome_dropped covers KNOWN_OUTCOME_OBS_COLUMNS only. "
            "unclassified_dropped is everything else that arrived -- technical "
            "fields and any clinical label nobody has adjudicated alike. Both "
            "were dropped; the split describes what is known about them, not "
            "what was done to them."),
        "arrived": list(arrived),
        "kept": list(kept),
        "known_outcome_dropped": list(known_outcome_dropped),
        "unclassified_dropped": list(unclassified_dropped),
    }
    write_manifest()


def record_outcome_analysis(name):
    """Append one executed outcome analysis and flush AT ONCE.

    Holding the list in memory until the end meant a job that completed the
    mortality correlation and then died in the gene re-audit left a manifest
    showing no outcome analysis at all. run_status stays `started` -- the run
    did not finish -- while the executed list is already true.
    """
    if name not in OUTCOME_ANALYSIS_RAN:
        OUTCOME_ANALYSIS_RAN.append(name)
    _MANIFEST["outcome_analyses_executed"] = list(OUTCOME_ANALYSIS_RAN)
    write_manifest()


write_manifest()
print(f"manifest written: run_status=started, integration_only={INTEGRATION_ONLY}", flush=True)

a, y_by_patient, qc_by_patient = load_canonical_cohort()

# ---------------------------------------------------------------------
# Doublet detection, per-sample (native scanpy batch handling).
# ---------------------------------------------------------------------
print("\nRunning sc.pp.scrublet (per-sample via batch_key)...", flush=True)
a.layers["counts"] = a.X.copy()
sample_sizes = a.obs["Sample_ID"].value_counts()
too_small = sample_sizes[sample_sizes < MIN_CELLS_PER_SAMPLE].index.tolist()
if too_small:
    print(f"  {len(too_small)} samples have <{MIN_CELLS_PER_SAMPLE} cells (Scrublet needs a minimum per batch); "
          f"these cells will get NaN doublet scores: {too_small[:10]}", flush=True)

# Every cell starts UNKNOWN. Nothing is presumed non-doublet.
a.obs["doublet_call"] = pd.array([pd.NA] * a.n_obs, dtype="boolean")
a.obs["examination_status"] = EXAM_NOT_EVALUATED
a.obs["doublet_score"] = np.nan

# `too_small` is decided ONLY by the cell-count rule computed above -- never
# inferred from a NaN score. A NaN score in a large batch means the detector
# failed, which is a different fact and must not be relabelled as "too small".
too_small_set = set(too_small)


def _finite_threshold(uns_for_sample):
    """The automatic threshold for one sample, or None with a reason.

    This is the check the stage-32 failure actually needs. Scanpy WARNS rather
    than raising when automatic thresholding fails, then returns finite scores
    and an all-False prediction -- output that passes every field-level test
    while resting on no threshold at all. An all-False prediction is NOT itself
    evidence of failure (it can be legitimate), so the threshold must be read.
    """
    if not isinstance(uns_for_sample, dict):
        return None, "no scrublet metadata for this sample"
    if "threshold" not in uns_for_sample:
        return None, "scrublet metadata carries no threshold"
    try:
        thr = float(uns_for_sample["threshold"])
    except (TypeError, ValueError):
        return None, f"threshold is not a number ({uns_for_sample['threshold']!r})"
    if not np.isfinite(thr):
        return None, f"automatic threshold is not finite ({thr})"
    return thr, None


def evaluate_sample(sub_obs, sid, n_cells, uns_for_sample):
    """Decide one sample's outcome. PURE -- returns what to write, writes nothing.

    Returns (status, calls_or_None, scores_or_None). Only EXAM_EXAMINED comes
    with a call and a score; every other status returns (status, None, None) so
    a caller cannot accidentally keep stale values alongside a failure.

    Three things must hold before a sample counts as examined: a finite
    automatic threshold from the sample's own metadata, strictly boolean calls
    of the right length, and finite 1-D scores of the right length. Scrublet
    returning without raising is not proof of a usable result.
    """
    thr, why = _finite_threshold(uns_for_sample)
    if thr is None:
        return EXAM_THRESHOLD_FAILED, None, None, why

    calls, why = strict_boolean_calls(sub_obs.get("predicted_doublet"), n_cells)
    if calls is None:
        return EXAM_RUNTIME_FAILED, None, None, f"unusable predicted_doublet -- {why}"

    # Fail closed on a missing column: np.asarray(None, dtype=float) is 0-d, so
    # `.shape[0]` raised IndexError from outside the per-sample try/except and
    # would have killed the whole run instead of marking one sample unknown.
    raw_scores = sub_obs.get("doublet_score")
    if raw_scores is None:
        return EXAM_RUNTIME_FAILED, None, None, "doublet_score column absent"
    try:
        scores = np.asarray(raw_scores, dtype=float)
    except Exception as exc:
        return EXAM_RUNTIME_FAILED, None, None, f"doublet_score unreadable ({type(exc).__name__})"
    # NOT .ravel(): a 2-D return means the interface or schema is wrong, and
    # flattening would guess an element order. Calls are held to ndim == 1, and
    # scores must be too.
    if scores.ndim != 1:
        return EXAM_RUNTIME_FAILED, None, None, f"doublet_score is not 1-D (ndim={scores.ndim})"
    if scores.shape[0] != n_cells:
        return (EXAM_RUNTIME_FAILED, None, None,
                f"doublet_score length {scores.shape[0]} != {n_cells}")
    if not bool(np.all(np.isfinite(scores))):
        return EXAM_RUNTIME_FAILED, None, None, "doublet_score contains non-finite values"

    return EXAM_EXAMINED, calls, scores, None


def reset_detection_state(obs):
    """Return every detection field to its untouched value.

    Called before the per-sample fallback so partial results written by a
    batched run that later raised cannot survive into it. Without this, a
    sample accepted by the batched loop and then failed by the fallback would
    keep its old call alongside a failure status -- a contradiction the
    validator would (correctly) refuse, aborting the whole run.
    """
    obs["doublet_call"] = pd.array([pd.NA] * len(obs), dtype="boolean")
    obs["doublet_score"] = np.nan
    obs["examination_status"] = EXAM_NOT_EVALUATED
    return obs


def _commit(staged):
    """Write staged per-sample outcomes into a.obs in one pass."""
    for idx, status, calls, scores in staged:
        a.obs.loc[idx, "examination_status"] = status
        if calls is not None:
            a.obs.loc[idx, "doublet_call"] = calls
        if scores is not None:
            a.obs.loc[idx, "doublet_score"] = scores


sample_groups = list(a.obs.groupby("Sample_ID").groups.items())

# ---------------------------------------------------------------------
# Dimension feasibility, decided BEFORE Scrublet is called.
#
# Computed once and shared by the batched path and the per-sample fallback, so
# the two cannot disagree about which samples were askable at N_PRIN_COMPS.
# Infeasible samples are withheld from both: they are not a detector failure,
# they are a question the data cannot support at this dimensionality.
# ---------------------------------------------------------------------
print("Measuring per-sample dimensions with scanpy's own scrublet "
      "preprocessing..." if DIMENSION_AUDIT else
      f"Checking dimension feasibility at n_prin_comps={N_PRIN_COMPS} "
      f"(scanpy's own scrublet preprocessing)...", flush=True)
DIAG = {}
insufficient = set()
measurement_failed = set()
#: Per-sample threshold and failure reason as OBSERVED by whichever path ran.
#: Filled by both the batched loop and the fallback so the diagnostics table
#: describes the run that happened, not the one that was planned.
SAMPLE_THRESHOLD = {}
SAMPLE_FAILURE = {}
#: Samples that entered sc.pp.scrublet, whatever came back.
SAMPLE_ATTEMPTED = set()
for sid, idx in sample_groups:
    n_cells = len(idx)
    row = {"Sample_ID": sid, "n_cells": n_cells,
           "n_input_features": int(a.n_vars),
           "n_cells_effective": None, "n_usable_features": None,
           "n_prin_comps_requested": N_PRIN_COMPS,
           # scrublet_attempted records whether the sample ENTERED the
           # detector. A sample attempted at the fixed dimension and then
           # rejected on its threshold did run at that dimension, and the
           # diagnostics must say so: reporting an empty dimension there
           # described a run that did not happen. Only PRE-CALL rejections --
           # too_small and insufficient_dimension -- have no dimension.
           "scrublet_attempted": False, "n_prin_comps_used_for_attempt": None,
           "threshold": None, "examination_status": None,
           "meets_min_cells_rule": None, "preflight_exclusion_reason": None,
           "dimension_measurement_error": None,
           "failure_reason": None, "n_examined": 0, "n_doublet": 0}
    # Measurement is not gated on eligibility. Skipping the <MIN_CELLS samples
    # here made "the per-sample dimension distribution" mean "the distribution
    # over the samples the CURRENT cell-count rule already admits" -- a
    # detection-stage decision leaking into the measurement that is supposed to
    # inform the next one. Every sample is measured; whether it meets the rule
    # is recorded beside the measurement.
    row["meets_min_cells_rule"] = sid not in too_small_set
    if sid in too_small_set:
        row["preflight_exclusion_reason"] = (
            f"n_cells={n_cells} < MIN_CELLS_PER_SAMPLE={MIN_CELLS_PER_SAMPLE}")
        row["failure_reason"] = row["preflight_exclusion_reason"]
    # Bound before the try. Without this, an eligible sample whose measurement
    # raised fell through to dimension_feasibility(n_eff, ...) with the names
    # unbound -- an UnboundLocalError that also destroyed the recorded reason.
    n_eff = n_usable = None
    try:
        n_eff, n_usable = scrublet_preprocessing_dimensions(a[idx])
        row["n_cells_effective"], row["n_usable_features"] = n_eff, n_usable
    except Exception as _e:
        # Measuring a sample can fail in the preprocessing itself. That is a
        # fact about the sample, recorded, not a reason to abort the run.
        n_eff = n_usable = None
        row["dimension_measurement_error"] = f"{type(_e).__name__}: {_e}"
    if sid in too_small_set:
        DIAG[sid] = row
        continue
    if n_eff is None:
        # An ELIGIBLE sample we could not measure. There is nothing to judge
        # feasibility against, so it reaches neither the verdict nor Scrublet.
        measurement_failed.add(sid)
        row["failure_reason"] = row["dimension_measurement_error"]
        DIAG[sid] = row
        continue
    if not DIMENSION_AUDIT:
        # No verdict in the preflight: there is no value to judge against, and
        # inventing one would defeat the point of measuring first.
        ok, why = dimension_feasibility(n_eff, n_usable, N_PRIN_COMPS)
        if not ok:
            insufficient.add(sid)
            row["failure_reason"] = why
    DIAG[sid] = row
if DIMENSION_AUDIT:
    # ---------------------------------------------------------------------
    # PREFLIGHT ENDS HERE. Scrublet is never called, no cell is called, no rate
    # and no contamination table is produced. The output is the distribution
    # the dimensionality decision is supposed to rest on.
    # ---------------------------------------------------------------------
    _audit_cols = ["Sample_ID", "n_cells", "n_input_features",
                   "n_cells_effective", "n_usable_features",
                   "meets_min_cells_rule", "preflight_exclusion_reason",
                   "dimension_measurement_error"]
    _audit = pd.DataFrame([DIAG[sid] for sid, _ in sample_groups])[_audit_cols]
    _measured = _audit[_audit["n_cells_effective"].notna()]
    write_csv_atomic(_audit, os.path.join(OUT_DIR, "dimension_audit.csv"), index=False)
    print(f"SAVED {os.path.join(OUT_DIR, 'dimension_audit.csv')} "
          f"({len(_audit)} samples, {len(_measured)} measured, "
          f"{int(_audit['meets_min_cells_rule'].sum())} meeting the current "
          f"cell-count rule)", flush=True)

    # TWO distributions, reported separately and labelled. They are not
    # interchangeable, and the parameter decision must say which one it used.
    for label, frame in (("ALL measured samples", _measured),
                         ("measured AND meeting the current cell-count rule",
                          _measured[_measured["meets_min_cells_rule"]])):
        if not len(frame):
            print(f"  {label}: none", flush=True)
            continue
        _lim = frame[["n_cells_effective", "n_usable_features"]].min(axis=1)
        print(f"  min(n_cells_effective, n_usable_features) over {len(frame)} "
              f"samples -- {label}:", flush=True)
        print("    " + "  ".join(f"p{q}={np.percentile(_lim, q):.0f}"
                                 for q in (0, 5, 25, 50, 75, 95, 100)), flush=True)
    print("  A fixed n_prin_comps admits a sample only where that minimum is "
          "STRICTLY greater than it.", flush=True)
    # ---------------------------------------------------------------------
    # Acceptance. Catching every per-sample measurement error and then writing
    # `completed` unconditionally meant a preflight in which all 301 samples
    # failed would still exit 0 and look usable, with the errors buried in a
    # CSV column. A parameter decision cannot rest on that.
    #
    # The two distributions have different requirements, so sufficiency is
    # reported per distribution rather than as one boolean.
    # ---------------------------------------------------------------------
    _n_total = len(_audit)
    _n_measured = len(_measured)
    _n_errors = int(_audit["dimension_measurement_error"].notna().sum())
    _elig = _audit[_audit["meets_min_cells_rule"] == True]          # noqa: E712
    _n_elig = len(_elig)
    _n_elig_measured = int(_elig["n_cells_effective"].notna().sum())
    _suff_all = (_n_measured == _n_total) and _n_total > 0
    _suff_elig = (_n_elig_measured == _n_elig) and _n_elig > 0
    _MANIFEST["dimension_audit"] = {
        "n_samples_total": _n_total,
        "n_samples_measured": _n_measured,
        "n_measurement_errors": _n_errors,
        "n_eligible_samples": _n_elig,
        "n_eligible_measured": _n_elig_measured,
        "dimension_audit_sufficient_for_decision": {
            "all_samples_distribution": bool(_suff_all),
            "eligible_only_distribution": bool(_suff_elig),
        },
    }
    # `completed` means every sample was measured. Anything less is incomplete,
    # whatever the exit code would otherwise suggest.
    _MANIFEST["run_status"] = "completed" if _suff_all else "incomplete"
    _MANIFEST["finished"] = _dt.datetime.now().isoformat(timespec="seconds")
    _MANIFEST["scrublet_attempted"] = False
    _MANIFEST["scrublet_any_call_completed"] = False
    write_manifest()
    print(f"\n  measured {_n_measured}/{_n_total} samples "
          f"({_n_errors} measurement error(s)); "
          f"eligible {_n_elig_measured}/{_n_elig}", flush=True)
    print(f"  sufficient for a decision on the ALL-samples distribution: {_suff_all}", flush=True)
    print(f"  sufficient for a decision on the ELIGIBLE-only distribution: {_suff_elig}",
          flush=True)
    # write_manifest() above has already returned, so os.replace() has
    # succeeded and the exit below cannot cost the record of why. The claim is
    # scoped to THIS control flow -- the exit happens after the replace -- and
    # not to signals or crashes generally, which could still arrive before it.
    # write_manifest does not fsync; an earlier comment here claimed it did,
    # carried over from gate2_checkpoint.py without that module's fsync
    # implementation.
    _exit = preflight_exit_code(_suff_all, _suff_elig)
    if _exit != PREFLIGHT_EXIT_OK:
        print("  NEITHER distribution is complete. This preflight does not support a "
              "dimensionality decision.", flush=True)
        print(f"ERROR: dimension audit incomplete -- {_n_errors} of {_n_total} samples "
              f"could not be measured, including {_n_elig - _n_elig_measured} of "
              f"{_n_elig} that meet the cell-count rule. No dimensionality decision may "
              f"be taken from this run. Inspect "
              f"{os.path.join(OUT_DIR, 'dimension_audit.csv')} (column "
              f"dimension_measurement_error) and {MANIFEST_PATH} "
              f"(dimension_audit.dimension_audit_sufficient_for_decision).",
              file=sys.stderr, flush=True)
    elif not _suff_all:
        print("  Only the ELIGIBLE-ONLY distribution is complete. A decision taken from "
              "this run must state that it used that distribution.", flush=True)
    print(f"\nmanifest updated: run_status={_MANIFEST['run_status']}, "
          f"mode=dimension_audit, scrublet_attempted=false, "
          f"scrublet_any_call_completed=false", flush=True)
    print("No doublet call, rate or contamination result was produced. "
          "The dimensionality decision is still open "
          "(docs/scrublet_dimension_parameter.md).", flush=True)
    raise SystemExit(_exit)

if insufficient:
    print(f"  {len(insufficient)} sample(s) cannot support n_prin_comps={N_PRIN_COMPS}; "
          f"marked {EXAM_INSUFFICIENT_DIMENSION}, NOT runtime_failed", flush=True)

#: Cells belonging to samples that are askable at N_PRIN_COMPS. Withholding the
#: rest keeps an infeasible batch from taking the whole batched run down with
#: it -- which is what happened on the v1 run.
askable_idx = [i for sid, idx in sample_groups
               if sid not in too_small_set and sid not in insufficient
               and sid not in measurement_failed
               for i in idx]
if measurement_failed:
    print(f"  {len(measurement_failed)} eligible sample(s) could not be MEASURED; "
          f"marked {EXAM_DIMENSION_MEASUREMENT_FAILED}", flush=True)
print(f"  {len(sample_groups) - len(too_small_set) - len(insufficient) - len(measurement_failed)}"
      f" of {len(sample_groups)} samples are askable ({len(askable_idx)} cells)", flush=True)

try:
    if not askable_idx:
        raise RuntimeError("no sample is askable at the fixed dimensionality")
    scored = a[askable_idx].copy()
    # Every askable sample enters the detector in this one call -- recorded
    # BEFORE it, so a sample is not claimed un-attempted because the call raised.
    SAMPLE_ATTEMPTED.update(sid for sid, _ in sample_groups
                            if sid not in too_small_set and sid not in insufficient
                            and sid not in measurement_failed)
    _MANIFEST["scrublet_attempted"] = True
    write_manifest()
    sc.pp.scrublet(scored, batch_key="Sample_ID", verbose=False, **SCRUBLET_PARAMS)
    _MANIFEST["scrublet_any_call_completed"] = True
    write_manifest()
    # scanpy stores per-batch metadata under uns['scrublet']['batches'][batch]
    # when batch_key is given. Without it there is no sample-level threshold to
    # verify, and a sample with no metadata must stay unknown rather than be
    # accepted on the strength of the pooled run.
    batches = (scored.uns.get("scrublet", {}) or {}).get("batches", {})
    if not batches:
        raise RuntimeError("batched scrublet returned no per-batch metadata; "
                           "sample-level thresholds cannot be verified")
    # Evaluate every sample into a staging list FIRST. Writing as we went meant
    # a later exception dropped us into the fallback with part of a.obs already
    # populated, so a sample accepted here and failed there kept a stale call
    # beside a failure status. Nothing touches a.obs until the whole loop
    # succeeds.
    staged = []
    for sid, idx in sample_groups:
        if sid in too_small_set:
            staged.append((idx, EXAM_TOO_SMALL, None, None))
            continue
        if sid in insufficient:
            staged.append((idx, EXAM_INSUFFICIENT_DIMENSION, None, None))
            continue
        if sid in measurement_failed:
            staged.append((idx, EXAM_DIMENSION_MEASUREMENT_FAILED, None, None))
            continue
        thr, _ = _finite_threshold(batches.get(sid))
        SAMPLE_THRESHOLD[sid] = thr
        status, calls, scores, why = evaluate_sample(
            scored.obs.loc[idx], sid, len(idx), batches.get(sid))
        if why:
            SAMPLE_FAILURE[sid] = why
            print(f"    sample {sid}: {why}; left UNKNOWN", flush=True)
        staged.append((idx, status, calls, scores))
    _commit(staged)
except Exception as e:
    print(f"  batched scrublet failed ({e}); falling back to per-sample loop", flush=True)
    # Discard anything the batched attempt may have written before failing --
    # the per-sample records too, or a threshold from the abandoned batched run
    # would be reported beside a status the fallback assigned.
    reset_detection_state(a.obs)
    SAMPLE_THRESHOLD.clear()
    SAMPLE_FAILURE.clear()
    SAMPLE_ATTEMPTED.clear()
    staged = []
    for sid, idx in sample_groups:
        if sid in too_small_set:
            # NOT evaluated. Stays NA -- writing False here would assert
            # "examined and not a doublet" about cells never shown to Scrublet.
            staged.append((idx, EXAM_TOO_SMALL, None, None))
            continue
        if sid in insufficient:
            # Same feasibility verdict as the batched path -- one fixed
            # dimension, one decision, never re-derived here.
            staged.append((idx, EXAM_INSUFFICIENT_DIMENSION, None, None))
            continue
        if sid in measurement_failed:
            staged.append((idx, EXAM_DIMENSION_MEASUREMENT_FAILED, None, None))
            continue
        sub = a[idx].copy()
        SAMPLE_ATTEMPTED.add(sid)
        if not _MANIFEST["scrublet_attempted"]:
            _MANIFEST["scrublet_attempted"] = True
            write_manifest()
        try:
            sc.pp.scrublet(sub, verbose=False, **SCRUBLET_PARAMS)
            if not _MANIFEST["scrublet_any_call_completed"]:
                _MANIFEST["scrublet_any_call_completed"] = True
                write_manifest()
        except Exception as e2:
            SAMPLE_FAILURE[sid] = f"scrublet raised: {e2}"
            print(f"    sample {sid} ({sub.n_obs} cells) failed: {e2}", flush=True)
            staged.append((idx, EXAM_RUNTIME_FAILED, None, None))
            continue
        _uns = sub.uns.get("scrublet", {}) or {}
        thr, _ = _finite_threshold(_uns)
        SAMPLE_THRESHOLD[sid] = thr
        status, calls, scores, why = evaluate_sample(
            sub.obs, sid, sub.n_obs, _uns)
        if why:
            SAMPLE_FAILURE[sid] = why
            print(f"    sample {sid}: {why}; left UNKNOWN", flush=True)
        staged.append((idx, status, calls, scores))
    _commit(staged)

# Fail closed before anything reads these columns.
validate_examination_state(a.obs, where="after doublet detection")

counts = summarize_examination(a.obs)
print(f"Doublet detection done: {counts['n_doublet']} doublets among "
      f"{counts['n_examined']} EXAMINED cells "
      f"({100 * counts['n_doublet'] / max(counts['n_examined'], 1):.1f}%); "
      f"{counts['n_unknown']} cells UNKNOWN "
      f"(coverage {100 * counts['coverage']:.1f}%)", flush=True)
for st, n in counts["by_status"].items():
    print(f"    {st:<16} {n}", flush=True)

# doublet_call is a nullable boolean: True / False / <NA>. A blank is UNKNOWN
# and must never be read as False.
_cell_cols = ["patient_id", "Sample_ID", "doublet_score", "doublet_call",
              "examination_status"]
if not INTEGRATION_ONLY:
    _cell_cols.insert(2, "y")
write_csv_atomic(a.obs[_cell_cols], os.path.join(OUT_DIR, "doublet_scores_per_cell.csv"))
print(f"SAVED {os.path.join(OUT_DIR, 'doublet_scores_per_cell.csv')}", flush=True)

# ---------------------------------------------------------------------
# Per-sample diagnostics. The dimensions each sample actually had, the
# threshold it actually got, and why it was not examined when it was not --
# so a reviewer can judge the detector's working point per sample instead of
# inferring it from an aggregate rate.
# ---------------------------------------------------------------------
for sid, idx in sample_groups:
    row = DIAG[sid]
    sub = a.obs.loc[idx]
    st = sub["examination_status"].iloc[0]
    row["examination_status"] = st
    row["threshold"] = SAMPLE_THRESHOLD.get(sid)
    # Attempted, not examined: a threshold or runtime failure still ran at the
    # fixed dimension. Only a pre-call rejection leaves the dimension empty.
    row["scrublet_attempted"] = sid in SAMPLE_ATTEMPTED
    if sid in SAMPLE_ATTEMPTED:
        row["n_prin_comps_used_for_attempt"] = N_PRIN_COMPS
    else:
        assert st in (EXAM_TOO_SMALL, EXAM_INSUFFICIENT_DIMENSION,
                      EXAM_DIMENSION_MEASUREMENT_FAILED), (
            f"{sid}: status {st} without entering Scrublet; only pre-call "
            f"rejections may be un-attempted")
    if st == EXAM_EXAMINED:
        row["n_examined"] = int(len(sub))
        row["n_doublet"] = int((sub["doublet_call"] == True).sum())  # noqa: E712
    if row["failure_reason"] is None and st != EXAM_EXAMINED:
        row["failure_reason"] = SAMPLE_FAILURE.get(sid, st)

_diag_cols = ["Sample_ID", "n_cells", "n_input_features", "n_cells_effective",
              "n_usable_features", "n_prin_comps_requested", "scrublet_attempted",
              "n_prin_comps_used_for_attempt", "threshold", "examination_status",
              "failure_reason", "n_examined", "n_doublet"]
_diag = pd.DataFrame([DIAG[sid] for sid, _ in sample_groups])[_diag_cols]

# The diagnostics table and the cell-level state must agree, or one of them is
# describing a run that did not happen.
_by_status = _diag.groupby("examination_status")["n_cells"].sum().to_dict()
for st, n in counts["by_status"].items():
    assert int(_by_status.get(st, 0)) == int(n), (
        f"diagnostics disagree with cell-level counts for {st}: "
        f"{_by_status.get(st, 0)} vs {n}")
assert int(_diag["n_examined"].sum()) == int(counts["n_examined"])
assert int(_diag["n_doublet"].sum()) == int(counts["n_doublet"])

write_csv_atomic(_diag, os.path.join(OUT_DIR, "per_sample_diagnostics.csv"), index=False)
print(f"SAVED {os.path.join(OUT_DIR, 'per_sample_diagnostics.csv')} "
      f"({len(_diag)} samples)", flush=True)

# ---------------------------------------------------------------------
# Doublet rate per patient -- does it itself track depth/mortality?
# ---------------------------------------------------------------------
per_patient = doublet_rate_per_patient(a.obs)
qc = qc_by_patient.reindex(per_patient.index)
# y is built ONLY in full mode, inside the branch that uses it.

# A patient with zero examined cells has an UNDEFINED rate, not a 0% rate.
valid = per_patient["doublet_rate"].notna()
n_total, n_valid = len(per_patient), int(valid.sum())
print(f"\nPatients: {n_total} total, {n_valid} with >=1 examined cell, "
      f"{n_total - n_valid} with NO examined cells (rate undefined, EXCLUDED)", flush=True)

# Depth is a technical covariate and is always computed. The MORTALITY
# correlation is an outcome analysis and is skipped under --integration-only,
# so the semantics of the nullable state can be validated without executing it.
rho_depth, p_depth = spearmanr(per_patient.loc[valid, "doublet_rate"],
                               qc.loc[valid, "median_counts"])
print(f"Doublet rate vs depth     : rho={rho_depth:+.3f} (p={p_depth:.3f}), n={n_valid} patients",
      flush=True)

out = per_patient.copy()
out["median_counts"] = qc["median_counts"].values
out["included_in_spearman"] = valid.values

if INTEGRATION_ONLY:
    print("Doublet rate vs mortality : SKIPPED (--integration-only)", flush=True)
else:
    y = pd.Series({p: y_by_patient[p] for p in per_patient.index})
    rho_y, p_y = spearmanr(per_patient.loc[valid, "doublet_rate"], y[valid])
    print(f"Doublet rate vs mortality : rho={rho_y:+.3f} (p={p_y:.3f}), n={n_valid} patients",
          flush=True)
    out["y"] = y.values
    record_outcome_analysis("doublet_rate_vs_mortality")

write_csv_atomic(out.reset_index(), os.path.join(OUT_DIR, "doublet_rate_per_patient.csv"),
                 index=False)

# ---------------------------------------------------------------------
# Re-audit the flagged genes EXCLUDING predicted doublets.
# ---------------------------------------------------------------------
print("\n=== Re-audit excluding predicted doublets ===", flush=True)
# Remove ONLY confirmed doublets. `.fillna(False)` here previously turned every
# never-examined cell into a non-doublet; the cells were kept either way, but
# a_clean then travelled downstream as a "doublet-cleaned" set in which
# unexamined cells were indistinguishable from examined non-candidates.
# na_value=False states explicitly that an UNKNOWN cell is not dropped.
# It is not a claim that the cell is a non-doublet -- examination_status
# keeps that distinction, and the kept set is not "doublet-cleaned".
drop = (a.obs["doublet_call"] == True).to_numpy(dtype=bool, na_value=False)  # noqa: E712
a_clean = a[~drop].copy()
# Nullable state must survive the slice/copy: a boolean extension dtype can
# quietly degrade to object through AnnData subsetting, which would let a
# downstream consumer see "False"-looking values again.
validate_examination_state(a_clean.obs, where="a_clean after slicing")
assert pd.api.types.is_bool_dtype(a_clean.obs["doublet_call"].dtype), (
    f"doublet_call degraded to {a_clean.obs['doublet_call'].dtype} through the AnnData slice")

# And it must survive a real write/read cycle. Done on a TEMPORARY file that is
# deleted immediately -- this is a check, not a result artifact.
_rt = os.path.join(OUT_DIR, "_nullable_roundtrip_check.h5ad")
try:
    a_clean.write_h5ad(_rt)
    _back = sc.read_h5ad(_rt)
    assert pd.api.types.is_bool_dtype(_back.obs["doublet_call"].dtype), (
        f"doublet_call came back from disk as {_back.obs['doublet_call'].dtype}")
    assert int(_back.obs["doublet_call"].isna().sum()) == \
        int(a_clean.obs["doublet_call"].isna().sum()), "NA count changed across the h5ad round trip"
    assert (_back.obs["examination_status"].astype(str).to_numpy()
            == a_clean.obs["examination_status"].astype(str).to_numpy()).all()
    validate_examination_state(_back.obs, where="a_clean after an h5ad round trip")
    print("  nullable state survives AnnData slice AND h5ad round trip", flush=True)
    del _back
finally:
    if os.path.exists(_rt):
        os.remove(_rt)

n_unknown_kept = int((a_clean.obs["examination_status"] != EXAM_EXAMINED).sum())
print(f"Cells after removing confirmed doublets: {a_clean.n_obs} / {a.n_obs} "
      f"({int(drop.sum())} removed)", flush=True)
print(f"  of which {n_unknown_kept} were NEVER EXAMINED -- retained, and still "
      f"distinguishable via examination_status; this set is NOT 'doublet-cleaned'",
      flush=True)

var_names = a_clean.var_names.astype(str).to_numpy()
gene_to_col = {g: i for i, g in enumerate(var_names)}
found = {g: gene_to_col[g] for g in FLAGGED_GENES if g in gene_to_col}

X = a_clean.X
patient_ids = a_clean.obs["patient_id"].to_numpy()
unique_patients = sorted(set(patient_ids))

rows = []
for gene, col in found.items():
    col_data = X[:, col].toarray().ravel() if hasattr(X[:, col], "toarray") else np.asarray(X[:, col]).ravel()
    per_patient_mean = pd.Series(col_data, index=patient_ids).groupby(level=0).mean()
    per_patient_mean = per_patient_mean.reindex(unique_patients)
    qc_s = qc_by_patient.reindex(unique_patients)

    # Technical covariates always; the mortality correlation only when an
    # outcome analysis was authorized.
    rho_counts, p_counts = spearmanr(per_patient_mean, qc_s["median_counts"])
    rho_mito, p_mito = spearmanr(per_patient_mean, qc_s["median_pct_mito"])
    row = {"gene": gene, "lineage": FLAGGED_GENES[gene],
           "rho_vs_depth_doublet_excluded": rho_counts, "p_vs_depth_doublet_excluded": p_counts,
           "rho_vs_pct_mito_doublet_excluded": rho_mito, "p_vs_pct_mito_doublet_excluded": p_mito}

    if INTEGRATION_ONLY:
        print(f"{gene:10s} ({FLAGGED_GENES[gene]:12s}): "
              f"rho_vs_depth={rho_counts:+.3f} (p={p_counts:.3f})  "
              f"[mortality SKIPPED]", flush=True)
    else:
        y_s = pd.Series({p: y_by_patient[p] for p in unique_patients})
        rho_y, p_y = spearmanr(per_patient_mean, y_s)
        record_outcome_analysis("flagged_gene_vs_mortality")
        row["rho_vs_mortality_doublet_excluded"] = rho_y
        row["p_vs_mortality_doublet_excluded"] = p_y
        print(f"{gene:10s} ({FLAGGED_GENES[gene]:12s}): "
              f"rho_vs_mortality={rho_y:+.3f} (p={p_y:.3f})  "
              f"rho_vs_depth={rho_counts:+.3f} (p={p_counts:.3f})", flush=True)
    rows.append(row)

df = pd.DataFrame(rows)
out_path = os.path.join(OUT_DIR, "contamination_audit_v2_doublet_excluded.csv")
write_csv_atomic(df, out_path, index=False)
print(f"\nSAVED {out_path}", flush=True)


# ---------------------------------------------------------------------
# Everything completed. Only now does the manifest say so, and only the
# analyses that actually ran are listed.
# ---------------------------------------------------------------------
_MANIFEST["run_status"] = "completed"
_MANIFEST["finished"] = _dt.datetime.now().isoformat(timespec="seconds")
_MANIFEST["outcome_analyses_executed"] = list(OUTCOME_ANALYSIS_RAN)
_MANIFEST["examination_summary"] = {k: v for k, v in counts.items() if k != "by_status"}
_MANIFEST["examination_by_status"] = {str(k): int(v) for k, v in counts["by_status"].items()}
write_manifest()
print(f"\nmanifest updated: run_status=completed, "
      f"outcome_analyses_executed={OUTCOME_ANALYSIS_RAN or 'none'}", flush=True)
