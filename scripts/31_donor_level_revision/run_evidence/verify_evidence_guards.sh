#!/bin/bash
# Prove that each guard in capture_test_run_evidence.sh actually bites.
#
# Claiming "four injections were run" without the commands and their output is
# the same unbacked assertion this directory exists to stop. This script runs
# them, records the exit codes, and fails if any guard lets its fault through.
#
#   bash verify_evidence_guards.sh > 2026-09-22_guard_verification.txt
#   echo $?    # 0 = every guard bit; non-zero = a guard is decorative
#
# Each injection is applied to a COPY; the committed script is never modified.
#
# TWO REPORTING BUGS THIS HARNESS HAD, both of which turned "the injection did
# not apply" into "the guard is decorative":
#   * the old/new pair was one argument split on "||", which occurs inside the
#     shell code being injected, so split() raised -- and the injector's status
#     was not checked, so an UNMODIFIED copy ran and passed.
#   * the replacement took the first occurrence, which for the pytest command is
#     the ECHOED display block, not the executed one.
# Old and new are separate arguments now, the occurrence index is explicit, and
# a failed injection is fatal. "The guard did not bite" must never be reachable
# from "the injection did not apply".
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
SRC="$HERE/capture_test_run_evidence.sh"
WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
FAILURES=0

echo "=== guard verification  [$(date '+%Y-%m-%d %H:%M:%S %Z %z')] ==="
echo "repo    : $REPO"
echo "HEAD    : $(git -C "$REPO" rev-parse HEAD)"
echo "capture : $SRC"
echo "sha256  : $(shasum -a 256 "$SRC" | cut -d' ' -f1)"

inject () {   # $1 file  $2 old  $3 new  $4 occurrence
  python3 - "$1" "$2" "$3" "$4" <<'PYEOF'
import sys
path, old, new, nth = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
s = open(path).read()
n = s.count(old)
if n < nth:
    sys.exit(f"INJECTION TARGET NOT FOUND ({n} occurrence(s), wanted #{nth}): {old[:70]!r}")
i = -1
for _ in range(nth):
    i = s.index(old, i + 1)
open(path, 'w').write(s[:i] + new + s[i + len(old):])
PYEOF
}

run_probe () {   # $1 name  $2 want  $3 marker
  local name="$1" want="$2" marker="$3" got out
  printf '\n--- %s ---\n' "$name"
  echo "expected exit: $want"
  ( cd "$REPO" && bash "$WORK/cap.sh" ) > "$WORK/out.txt" 2>&1
  got=$?
  echo "actual exit  : $got"
  out=$(grep -m1 -E "$marker" "$WORK/out.txt") || out="(marker NOT FOUND: $marker)"
  echo "marker       : $out"
  if [ "$got" -ne "$want" ]; then
    echo "*** GUARD DID NOT BITE: expected $want, got $got"
    FAILURES=$((FAILURES + 1))
  fi
}

probe () {   # $1 name  $2 want  $3 old  $4 new  $5 nth  $6 marker
  cp "$SRC" "$WORK/cap.sh"
  if ! inject "$WORK/cap.sh" "$3" "$4" "$5"; then
    printf '\n--- %s ---\n*** INJECTION FAILED TO APPLY -- a result here would be meaningless\n' "$1"
    FAILURES=$((FAILURES + 1))
    return
  fi
  run_probe "$1" "$2" "$6"
}

probe "A. transport failure after six correct hashes are printed" 72 \
  'done" > "$raw" 2>&1' \
  'done; exit 55" > "$raw" 2>&1' \
  1 'REMOTE HASH QUERY FAILED'

probe "B. hash query truncated to five rows" 73 \
  'for f in ${FILES[*]}; do if' \
  'for f in ${FILES[*]:0:5}; do if' \
  1 'EXPECTED 6 hash rows'

# Occurrence 2: the first is the ECHOED display block. Injecting there changes
# only the printed text, the run succeeds, and the probe reports a false
# negative -- which is exactly what happened the first time.
probe "C. remote cd fails; pytest must not run in the wrong directory" 1 \
  'cd $REMOTE_ROOT
python3 -B -m pytest' \
  'cd $REMOTE_ROOT/no_such_dir
python3 -B -m pytest' \
  2 'no_such_dir|No such file'

printf '\n--- D. a tracked file differs from its HEAD blob ---\n'
echo "expected exit: 74"
echo '# injected drift, reverted immediately' >> "$REPO/DEVELOPMENT.md"
( cd "$REPO" && bash "$SRC" ) > "$WORK/out.txt" 2>&1; got=$?
git -C "$REPO" checkout -- DEVELOPMENT.md
echo "actual exit  : $got"
grep -m1 "differs from HEAD" "$WORK/out.txt" || echo "(marker NOT FOUND)"
if [ "$got" -ne 74 ]; then echo "*** GUARD DID NOT BITE"; FAILURES=$((FAILURES + 1)); fi

printf '\n--- E. remote content drifts after the run ---\n'
echo "expected exit: 71"
cp "$SRC" "$WORK/cap.sh"
if ! inject "$WORK/cap.sh" 'say "8. REMOTE hashes, AFTER the run' \
     'ssh -x "$CLUSTER" "echo '"'"'# drift'"'"' >> $REMOTE_ROOT/DEVELOPMENT.md"
say "8. REMOTE hashes, AFTER the run' 1; then
  echo "*** INJECTION FAILED TO APPLY"; FAILURES=$((FAILURES + 1))
else
  ( cd "$REPO" && bash "$WORK/cap.sh" ) > "$WORK/out.txt" 2>&1; got=$?
  echo "actual exit  : $got"
  grep -m1 "MISMATCH at after the run" "$WORK/out.txt" || echo "(marker NOT FOUND)"
  if [ "$got" -ne 71 ]; then echo "*** GUARD DID NOT BITE"; FAILURES=$((FAILURES + 1)); fi
fi

printf '\n--- repair: a clean capture run, whose sync restores the drifted remote file ---\n'
( cd "$REPO" && bash "$SRC" ) > "$WORK/clean.txt" 2>&1; got=$?
echo "clean capture exit: $got  (must be 0)"
grep -cE '^AGREE:' "$WORK/clean.txt" | xargs echo "AGREE lines in the clean run:"
if [ "$got" -ne 0 ]; then FAILURES=$((FAILURES + 1)); fi

printf '\n=== %s guard(s) failed to bite ===\n' "$FAILURES"
[ "$FAILURES" -eq 0 ]
