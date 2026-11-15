#!/bin/bash
# Evidence binding a reported test count to a known set of file bytes, and
# those bytes to a commit.
#
# HPCC is NOT a checkout of this branch: its git HEAD is 225 commits behind and
# scripts/31_donor_level_revision/ is untracked there. A remote `git rev-parse`
# identifies the remote repo's baseline, never the code a stage 31 run used.
# Only content hashes do.
#
# Chain, and every link fails the script rather than being reported:
#   HEAD blob == working tree  (exit 74)
#   -> sync -> remote == local (exit 71)
#   -> run -> remote == local AGAIN (exit 71)
# with the remote query's own transport status checked (exit 72) and its row
# count and path set validated (exit 73), because a truncated or failed query
# that still emitted the right text would otherwise read as agreement.
#
#   bash capture_test_run_evidence.sh > 2026-09-22_test_run.txt
#   echo $?     # non-zero => incomplete or void; do not cite
set -euo pipefail
CLUSTER=xinruqiu@cluster.hpcc.ucr.edu
REMOTE_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
ENVDIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
SUITE=tests/test_stage31_doublet_semantics.py
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

FILES=(
  scripts/31_donor_level_revision/code/doublet_detection_and_reaudit.py
  scripts/31_donor_level_revision/code/run_doublet_detection.sh
  tests/test_stage31_doublet_semantics.py
  docs/scrublet_dimension_parameter.md
  scripts/31_donor_level_revision/STALE_ARTIFACTS.md
  DEVELOPMENT.md
)

say () { printf '\n=== %s  [%s] ===\n' "$1" "$(date '+%Y-%m-%d %H:%M:%S %Z')"; }
local_hashes () { for f in "${FILES[@]}"; do printf '%s  %s\n' "$(shasum -a 256 "$f" | cut -d' ' -f1)" "$f"; done; }

# The remote query writes to a FILE, so its exit status is visible. An earlier
# version ran it inside a process substitution -- `diff <(local) <(remote|...)`
# -- where bash does not propagate the substitution's status, so an ssh that
# failed AFTER printing six correct hashes still produced "AGREE".
remote_hashes_to () {
  local raw="$1" rc=0
  ssh -x "$CLUSTER" "cd $REMOTE_ROOT && for f in ${FILES[*]}; do if [ -f \"\$f\" ]; then printf '%s  %s\n' \"\$(sha256sum \"\$f\" | cut -d' ' -f1)\" \"\$f\"; else printf 'MISSING  %s\n' \"\$f\"; fi; done" > "$raw" 2>&1 || rc=$?
  cat "$raw"
  if [ "$rc" -ne 0 ]; then
    printf '\n*** REMOTE HASH QUERY FAILED (ssh exit %s) -- evidence void\n' "$rc"
    exit 72
  fi
  # Structured rows only; banners stay in the printed raw output above.
  grep -E '^([0-9a-f]{64}|MISSING)  ' "$raw" > "$raw.rows" || true
  local n; n=$(wc -l < "$raw.rows" | tr -d ' ')
  if [ "$n" -ne "${#FILES[@]}" ]; then
    printf '\n*** EXPECTED %s hash rows, GOT %s -- evidence void\n' "${#FILES[@]}" "$n"
    exit 73
  fi
  if ! diff <(printf '%s\n' "${FILES[@]}" | sort) <(awk '{print $2}' "$raw.rows" | sort) > /dev/null; then
    printf '\n*** hash rows do not cover exactly the declared file set -- evidence void\n'
    exit 73
  fi
}

agree_or_die () {
  local stage="$1" raw="$2"
  if diff <(local_hashes) "$raw.rows"; then
    echo "AGREE: all ${#FILES[@]} files byte-identical local <-> HPCC ($stage)"
  else
    echo "*** MISMATCH at $stage -- this evidence is void, do not cite it"
    exit 71
  fi
}

cat <<HDR
================================================================================
Test-run evidence: which file bytes the reported pass count belongs to, and
which commit those bytes come from.
Produced by run_evidence/capture_test_run_evidence.sh, committed alongside.
Unfiltered combined stdout+stderr. Exit code non-zero => void, do not cite.
================================================================================
captured from : $(hostname)
captured at   : $(date '+%Y-%m-%d %H:%M:%S %Z %z')
local HEAD    : $(git rev-parse HEAD)
local dirty   : $(git status --porcelain | wc -l | tr -d ' ') file(s) in the repo overall
cluster       : ${CLUSTER#*@}
remote root   : $REMOTE_ROOT
suite         : $SUITE
files under test-identity check: ${#FILES[@]}

NOTE: "local dirty" counts the WHOLE repo. It says nothing about the six files
below, which is why step 0 checks them against HEAD explicitly -- an earlier
version reported dirty=4 and left the reader to assume the six were not among
them.
HDR

say "0. EACH FILE MUST EQUAL ITS HEAD BLOB (exit 74 otherwise)"
for f in "${FILES[@]}"; do
  git ls-files --error-unmatch -- "$f" > /dev/null
  blob=$(git show "HEAD:$f" | shasum -a 256 | cut -d' ' -f1)
  work=$(shasum -a 256 "$f" | cut -d' ' -f1)
  printf '%s  HEAD:%s\n%s  work:%s\n' "$blob" "$f" "$work" "$f"
  if [ "$blob" != "$work" ]; then
    echo "*** $f differs from HEAD -- the run would not belong to this commit"
    exit 74
  fi
done
echo "all ${#FILES[@]} files equal their HEAD blobs at $(git rev-parse HEAD)"

say "1. LOCAL hashes, before sync"
local_hashes
say "2. REMOTE hashes, before sync (raw, banners included)"
# No `|| true` here. It was added so a pre-sync difference would not abort --
# but a missing file yields a MISSING row that passes validation anyway, so the
# only thing it actually suppressed was a genuine transport or truncation
# failure at the first query.
remote_hashes_to "$TMP/pre"
say "3. SYNC"
for f in "${FILES[@]}"; do
  echo "\$ rsync -az --no-perms --no-owner --no-group $f $CLUSTER:$REMOTE_ROOT/$f"
  rsync -az --no-perms --no-owner --no-group "$f" "$CLUSTER:$REMOTE_ROOT/$f"
done
say "4. REMOTE hashes, after sync"
remote_hashes_to "$TMP/post"
say "5. AGREEMENT CHECK (exit 71 if any file differs)"
agree_or_die "after sync, before the run" "$TMP/post"
say "6. ENVIRONMENT"
# Commands go over stdin, not `bash -l -c '...'`. Nesting single quotes inside a
# single-quoted -c argument broke the f-string and the whole run exited 1 --
# which set -e did catch, so the mechanism worked even while the command did not.
#
# Both remote blocks carry `set -euo pipefail` of their own. Without it, a
# failing `cd $REMOTE_ROOT` would leave pytest running in the WRONG DIRECTORY,
# and ssh returns only the last command's status, so the earlier failure would
# vanish. Fail-closed at the hash query was not enough; the same class of gap
# existed here.
echo "\$ ssh $CLUSTER bash -l -s   # with the block below on stdin"
cat <<'ENVBLOCK'
set -euo pipefail
export PATH=$ENVDIR/bin:$PATH
python3 -V
python3 -c "import scanpy, anndata, numpy, pandas, pytest; print('scanpy', scanpy.__version__, 'anndata', anndata.__version__, 'numpy', numpy.__version__, 'pandas', pandas.__version__, 'pytest', pytest.__version__)"
ENVBLOCK
ssh -x "$CLUSTER" 'bash -l -s' <<ENVBLOCK 2>&1
set -euo pipefail
export PATH=$ENVDIR/bin:\$PATH
python3 -V
python3 -c "import scanpy, anndata, numpy, pandas, pytest; print('scanpy', scanpy.__version__, 'anndata', anndata.__version__, 'numpy', numpy.__version__, 'pandas', pandas.__version__, 'pytest', pytest.__version__)"
ENVBLOCK

say "7. TEST RUN"
echo "\$ ssh $CLUSTER bash -l -s   # with the block below on stdin"
cat <<TESTBLOCK
set -euo pipefail
export PATH=$ENVDIR/bin:\$PATH PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
cd $REMOTE_ROOT
python3 -B -m pytest $SUITE -q
TESTBLOCK
ssh -x "$CLUSTER" 'bash -l -s' <<TESTBLOCK 2>&1
set -euo pipefail
export PATH=$ENVDIR/bin:\$PATH PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
cd $REMOTE_ROOT
python3 -B -m pytest $SUITE -q
TESTBLOCK

say "8. REMOTE hashes, AFTER the run -- rules out the files moving under it"
remote_hashes_to "$TMP/final"
say "9. FINAL AGREEMENT CHECK (exit 71 if anything moved under the run)"
agree_or_die "after the run" "$TMP/final"
say "DONE"
