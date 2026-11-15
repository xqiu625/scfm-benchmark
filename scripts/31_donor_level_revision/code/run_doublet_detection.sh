#!/bin/bash
# Stage 31 doublet detection and contamination re-audit, under the FIXED
# nullable-state semantics (execution-plan P0b).
#
# The output directory is a REQUIRED argument and is never hardcoded here: the
# script also requires it, so there is one place it can be wrong rather than two
# that can disagree. The script refuses a directory that resolves to the legacy
# contamination_audit path, or that exists and is not empty.
#
#   sbatch run_doublet_detection.sh <NEW_OUT_DIR> <N_PRIN_COMPS> [--integration-only]
#   sbatch run_doublet_detection.sh <NEW_OUT_DIR> --dimension-audit --integration-only
#
# The second form is the PREFLIGHT: it attempts a measurement of
# n_cells_effective and n_usable_features for EVERY sample -- recording
# failures as failures -- writes dimension_audit.csv and stops. It exits 3 if
# neither reported distribution is complete, so the scheduler cannot mark a
# run that supports no dimensionality decision as COMPLETED. It
# takes no dimensionality, calls no detector and produces no call, rate or
# contamination output. It exists because requiring an approved dimensionality
# to run, while saying the choice should rest on this run's diagnostics, is
# circular -- the distribution has to be obtainable without first guessing.
#
# It is not a read-only run: it writes a manifest and one CSV into its own new
# isolated directory. What it leaves untouched is the inputs and every existing
# artifact. It refuses --n-prin-comps for governance reasons, not because the
# value would change what is measured -- it does not enter the preprocessing.
#
# N_PRIN_COMPS is REQUIRED and has no default here or in the script. The value
# is a methodological choice that has not been made
# (docs/scrublet_dimension_parameter.md is OPEN), so neither this wrapper nor
# the script may supply one -- a caller who has not been given a value by that
# decision cannot run the detector at all. Pass exactly what the approved
# command specifies.
#
# --integration-only runs detection, export and per-patient rates but SKIPS every
# mortality correlation, so the Scanpy / AnnData / CSV semantics can be validated
# without executing an outcome analysis. Omitting it means the run DOES execute
# the doublet-rate-versus-mortality analysis, which needs its own authorization.
#
# The outcome boundary that flag gives, stated exactly: outcome is excluded at
# PARSE time from the patient list (pandas usecols), while the raw h5ad
# container DOES carry mortality and severity obs columns into memory --
# read_h5ad has no column selection. Only the columns the loader requires
# survive it; the rest are dropped before any use and recorded in
# RUN_MANIFEST.json under container_obs_columns. The manifest therefore reports
#   known_outcome_values_loaded   true   (the container brought them in)
#   outcome_values_inspected      false  (nothing read their values)
#   outcome_analyses_executed     []
# "known_" is load-bearing: that flag speaks for the classified column list
# only, so an unrecognised clinical label shows up in
# container_obs_columns.unclassified_dropped instead. It is dropped either way
# -- retention, not the classification, is what prevents a leak.
# It is NOT "outcome is never loaded"; the strict form would need a separate
# outcome-free technical input, which does not exist.
#
# Resources: the input h5ad is a dense int64 47,977 x 13,066 matrix (~5 GB in
# memory). Peak comes from holding the full object plus the COVID subset copy,
# the subsample copy and Scrublet's own temporaries; the per-sample scoring
# itself is small (161 eligible samples of at most 60 cells each). 32G is a
# conservative first request -- check MaxRSS afterwards with
#   sacct -j <id> --format=JobID,MaxRSS,Elapsed,State
# and tighten it on the next run.

#SBATCH --job-name=doublet_detect_v2
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
# Log paths are RELATIVE, so SLURM resolves them against the submit directory
# and the logs land beside the checkout that actually ran. They used to name one
# specific checkout absolutely, which sent a run's logs into a tree that was not
# the tree being executed.
#
# Submit from <checkout>/scripts/31_donor_level_revision/ and create logs/ there
# first -- SBATCH directives are evaluated at submit time, so the directory must
# already exist; logs/ is untracked, so a fresh worktree will not have it.
#SBATCH --output=logs/doublet_detect_%j.out
#SBATCH --error=logs/doublet_detect_%j.err

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "usage: sbatch $0 <NEW_OUT_DIR> <N_PRIN_COMPS> [--integration-only]" >&2
    echo "  N_PRIN_COMPS is required and has no default; see" >&2
    echo "  docs/scrublet_dimension_parameter.md (OPEN)." >&2
    exit 2
fi
OUT_DIR="$1"; shift
if [ "$1" = "--dimension-audit" ]; then
    # Preflight: measures the per-sample dimensions and stops. It takes no
    # dimensionality precisely because it exists to inform that choice.
    N_PRIN_COMPS=""
else
    N_PRIN_COMPS="$1"; shift
    case "$N_PRIN_COMPS" in
        ''|*[!0-9]*) echo "N_PRIN_COMPS must be a positive integer, got '$N_PRIN_COMPS'" >&2; exit 2 ;;
    esac
fi

# ---------------------------------------------------------------------
# Which checkout is this run? Resolved, then VERIFIED.
#
# REPO used to be a hardcoded absolute path. On HPCC that path is a checkout
# whose git HEAD is 225 commits behind, with these two files sitting in it as
# UNTRACKED copies delivered by rsync. The bytes matched, but the Python entry
# point derives RUN_MANIFEST.json's git_commit from `git rev-parse HEAD` in its
# own directory -- so every run recorded a commit that did not contain the code
# it ran. Byte equality is not provenance.
#
# SLURM may stage the batch script elsewhere, so $0 is not reliable;
# SLURM_SUBMIT_DIR is where sbatch was invoked and is the contract: submit from
# the checkout you intend to run.
# ---------------------------------------------------------------------
_origin_dir="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO=$(cd "$_origin_dir" && git rev-parse --show-toplevel 2>/dev/null) || REPO=""
if [ -z "$REPO" ]; then
    echo "cannot resolve a git checkout from '$_origin_dir'." >&2
    echo "Submit from inside the checkout you intend to run." >&2
    exit 2
fi

ENTRY_REL=scripts/31_donor_level_revision/code/doublet_detection_and_reaudit.py
WRAPPER_REL=scripts/31_donor_level_revision/code/run_doublet_detection.sh
for _rel in "$WRAPPER_REL" "$ENTRY_REL"; do
    if ! git -C "$REPO" ls-files --error-unmatch -- "$_rel" > /dev/null 2>&1; then
        echo "$_rel is NOT TRACKED in $REPO -- refusing to run code that no commit contains." >&2
        exit 2
    fi
    if [ -n "$(git -C "$REPO" status --porcelain -- "$_rel")" ]; then
        echo "$_rel is modified in $REPO -- the run would not belong to its commit." >&2
        exit 2
    fi
done
REPO_COMMIT=$(git -C "$REPO" rev-parse HEAD)

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

mkdir -p "$REPO/scripts/31_donor_level_revision/logs"

echo "host      : $(hostname)"
echo "submit-dir: $_origin_dir"
echo "repo      : $REPO"
echo "commit    : $REPO_COMMIT"
echo "wrapper   : $(sha256sum "$REPO/$WRAPPER_REL" | cut -d' ' -f1)"
echo "entry     : $(sha256sum "$REPO/$ENTRY_REL" | cut -d' ' -f1)"
echo "out-dir   : $OUT_DIR"
echo "n_prin_comps: ${N_PRIN_COMPS:-(none -- preflight)}"
echo "extra args: $*"
python3 -c "import scanpy, anndata, numpy, pandas; print(
    f'scanpy {scanpy.__version__} anndata {anndata.__version__} '
    f'numpy {numpy.__version__} pandas {pandas.__version__}')"

python3 -B "$REPO/$ENTRY_REL" \
    --out-dir "$OUT_DIR" \
    ${N_PRIN_COMPS:+--n-prin-comps "$N_PRIN_COMPS"} "$@"
