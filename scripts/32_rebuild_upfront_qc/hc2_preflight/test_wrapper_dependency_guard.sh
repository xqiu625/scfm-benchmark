#!/bin/bash
# Regression test for review:20261012_132816 R2: the wrapper's tracked/
# staged/unstaged dependency guard must refuse when gate2_diagnostics.py (actually
# imported by gate2_stress_runner.py, confirmed by grep) is dirty -- staged or
# unstaged -- not just when the preflight script/wrapper itself is dirty.
#
# Revised 2026-10-12 per review:20261012_145004 F1: the previous version
# ran the REAL, unmodified wrapper for the two "must refuse" cases, trusting the
# guard-under-test to refuse before reaching `mkdir`/`python`. If the guard ever
# regressed, that version's own safety -- not just its correctness -- would have
# depended on the very code being tested, and a regressed guard would let the test
# itself create a real production output directory and invoke the real preflight
# script. Every case below now runs ONLY a sentinel copy of the wrapper, with the
# destructive tail (`mkdir -p "$OUT"` / `python "$SCRIPT" ...`) replaced by harmless
# echoes BEFORE any case executes it -- so even a fully disabled guard can only ever
# reach a sentinel echo, never a real directory or a real process. The guard logic
# itself (the `git`-based tracked/staged/unstaged checks) is left byte-identical.
#
# Runs entirely inside an isolated throwaway git repo under $(mktemp -d). Never
# touches the real repo, the real /bigdata deposit, or any real HPCC production
# path, regardless of whether the guard under test is working or broken.
#
#   bash test_wrapper_dependency_guard.sh
#
# Exits non-zero immediately on the first failure (not after accumulating).
set -euo pipefail

# Do not let an inherited SLURM_SUBMIT_DIR redirect the wrapper's `REPO=` resolution
# to some other, real directory -- this test must only ever operate on $TMP below.
unset SLURM_SUBMIT_DIR 2>/dev/null || true

# Resolved relative to this script's own location (one level under
# scripts/32_rebuild_upfront_qc/) rather than via `git rev-parse`, so this also
# runs correctly from a plain rsync'd mirror that is not itself a git repo (the
# synced files' CONTENT is what is under test here, not their git history).
HC2_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_SCRIPT="$HC2_DIR/capture_preflight_GSM5102901_HC2.py"
SRC_WRAPPER="$HC2_DIR/run_capture_preflight_HC2.sh"
SRC_CODE_DIR="$(cd "$HC2_DIR/../code" && pwd)"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
git init -q
git config user.email "test@test"
git config user.name "test"

mkdir -p scripts/32_rebuild_upfront_qc/hc2_preflight scripts/32_rebuild_upfront_qc/code logs
cp "$SRC_SCRIPT" scripts/32_rebuild_upfront_qc/hc2_preflight/
cp "$SRC_WRAPPER" scripts/32_rebuild_upfront_qc/hc2_preflight/
for f in gate2_parameter_ruling.py gate2_stress_design.py gate2_stress_runner.py gate2_diagnostics.py; do
  cp "$SRC_CODE_DIR/$f" scripts/32_rebuild_upfront_qc/code/
done
git add -A
git commit -q -m "clean baseline"

WRAPPER="scripts/32_rebuild_upfront_qc/hc2_preflight/run_capture_preflight_HC2.sh"
DIAG="scripts/32_rebuild_upfront_qc/code/gate2_diagnostics.py"

fail() { echo "FAIL: $1" >&2; [ -f out.log ] && cat out.log >&2; exit 1; }

# Builds a sentinel copy of $1 into $2: the destructive tail is replaced by
# harmless echoes; the guard logic above it is untouched. Verifies EXACTLY the
# 2 targeted lines changed, so this never silently also neuters the guard itself.
build_sentinel() {
  local src="$1" dst="$2"
  sed -e 's|^mkdir -p "\$OUT"$|echo "MKDIR_SENTINEL $OUT"|' \
      -e 's|^python "\$SCRIPT" --out-dir "\$OUT"$|echo "WORKLOAD_SENTINEL_REACHED"|' \
      "$src" > "$dst"
  local n_changed
  n_changed="$(diff "$src" "$dst" | grep -c '^[<>]' || true)"
  [ "$n_changed" = "4" ] || fail "sentinel build changed $n_changed diff-lines from $src, want exactly 4 (2 removed + 2 added) -- refusing to use an under- or over-modified sentinel"
}

SENTINEL="scripts/32_rebuild_upfront_qc/hc2_preflight/run_capture_preflight_HC2.sentinel.sh"
build_sentinel "$WRAPPER" "$SENTINEL"

run_sentinel() { bash "$SENTINEL" > out.log 2>&1; }

pass=0
expect_exit() {
  local want="$1" label="$2" got
  set +e
  run_sentinel; got=$?
  set -e
  [ "$got" = "$want" ] || fail "$label: exit=$got, want=$want"
  pass=$((pass+1))
}

# --- Case 1: unstaged change to gate2_diagnostics.py must be refused ---
echo "# extra line" >> "$DIAG"
expect_exit 2 "unstaged gate2_diagnostics.py change must refuse"
grep -q "REFUSING.*gate2_diagnostics.py.*unstaged" out.log \
  || fail "unstaged refusal message did not name gate2_diagnostics.py"
grep -q "WORKLOAD_SENTINEL_REACHED" out.log \
  && fail "unstaged-dirty case reached the (sentinel) workload -- guard did not refuse"
echo "OK: unstaged refusal message names gate2_diagnostics.py, workload not reached"
git checkout -q -- "$DIAG"

# --- Case 2: staged (but uncommitted) change to gate2_diagnostics.py must be refused ---
echo "# extra line" >> "$DIAG"
git add "$DIAG"
expect_exit 2 "staged gate2_diagnostics.py change must refuse"
grep -q "REFUSING.*gate2_diagnostics.py.*staged" out.log \
  || fail "staged refusal message did not name gate2_diagnostics.py"
grep -q "WORKLOAD_SENTINEL_REACHED" out.log \
  && fail "staged-dirty case reached the (sentinel) workload -- guard did not refuse"
echo "OK: staged refusal message names gate2_diagnostics.py, workload not reached"
git reset -q HEAD -- "$DIAG"
git checkout -q -- "$DIAG"

# --- Case 3: clean tree must pass the guard (reach past it) ---
expect_exit 0 "clean tree must pass the guard and reach the sentinel workload"
grep -q "WORKLOAD_SENTINEL_REACHED" out.log \
  || fail "clean tree did not reach the sentinel workload"
echo "OK: clean tree reaches past the guard to the (sentinel) workload"

# --- Mutation-kill: prove these assertions would actually catch a regressed guard,
# and that even a fully disabled guard still only ever reaches a sentinel echo. ---
MUTANT="scripts/32_rebuild_upfront_qc/hc2_preflight/run_capture_preflight_HC2.mutant.sh"
# Replace the guard's `for` loop body with a no-op comment, i.e. simulate the guard
# never refusing anything -- built FROM the sentinel (so its tail is already safe).
awk '
  /^for f in "\$SCRIPT" "\$WRAPPER" "\$\{IMPORTED\[@\]\}"; do$/ { print; print "  : # MUTATION: guard body disabled"; skip=1; next }
  skip && /^done$/ { print; skip=0; next }
  skip { next }
  { print }
' "$SENTINEL" > "$MUTANT"
grep -q "MUTATION: guard body disabled" "$MUTANT" \
  || fail "mutant construction did not find the guard for-loop to disable -- update the awk pattern to match the current wrapper"

echo "# extra line" >> "$DIAG"
set +e
bash "$MUTANT" > out.log 2>&1
mutant_rc=$?
set -e
[ "$mutant_rc" = "0" ] || fail "mutant (guard disabled) unexpectedly refused (exit=$mutant_rc) -- mutation did not take effect, cannot prove detection"
grep -q "WORKLOAD_SENTINEL_REACHED" out.log \
  || fail "mutant (guard disabled) did not reach even the sentinel workload -- mutation broke the script, not just the guard"
echo "OK: mutant (guard disabled) reaches the SENTINEL workload (not a real one) -- proving"
echo "    (a) the case-1/case-2 assertions above would have failed against this mutant"
echo "        (exit=0 instead of 2, workload reached), i.e. they really detect a"
echo "        regressed guard, and"
echo "    (b) even with the guard fully disabled, no real mkdir/python was ever invoked."
git checkout -q -- "$DIAG"
pass=$((pass+1))

echo "---"
echo "pass=$pass fail=0"
