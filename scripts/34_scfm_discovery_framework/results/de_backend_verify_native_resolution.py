"""Per review:20261015_150231 F1: `pip install -c constraints.txt` bounds the
VERSION of any package that gets pulled into a resolution; it is NOT a package-name
whitelist and does not itself refuse an unexpected extra package, nor does installing
`-c` by itself force every listed package to be installed. Comparing a REAL, executed
native-resolver dry-run report's own `install` list against the lock file's pinned set
-- exactly, both directions -- is what actually enforces "nothing outside the approved
closure gets installed, and nothing in it is silently skipped."

Usage: python de_backend_verify_native_resolution.py <pip_dry_run_report.json> <lock_file.txt>
Exit 0 only if the report's install list and the lock file's pinned set are IDENTICAL
(same package names, same versions, both directions). Exit 1 with a detailed mismatch
list otherwise. This does not install anything and does not require pydeseq2.
"""
import json
import sys


def _parse_lock_file(path):
    pins = {}
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if "==" not in line:
                raise SystemExit(f"{path}: line {line!r} is not pinned with '==' -- refusing an "
                                 "unpinned/range-pinned lock file")
            name, version = line.split("==", 1)
            pins[name.strip().lower()] = version.strip()
    return pins


def _parse_dry_run_report(path):
    with open(path) as f:
        report = json.load(f)
    installed = {}
    for item in report["install"]:
        md = item["metadata"]
        installed[md["name"].lower()] = md["version"]
    return installed


def compare(report_path, lock_path):
    locked = _parse_lock_file(lock_path)
    installed = _parse_dry_run_report(report_path)

    missing = sorted(set(locked) - set(installed))
    extra = sorted(set(installed) - set(locked))
    version_mismatches = sorted(
        name for name in (set(locked) & set(installed)) if locked[name] != installed[name])

    problems = []
    if missing:
        problems.append(f"locked but NOT in the native resolver's install list: {missing}")
    if extra:
        problems.append(f"in the native resolver's install list but NOT locked: {extra}")
    if version_mismatches:
        details = {name: (locked[name], installed[name]) for name in version_mismatches}
        problems.append(f"version differs (locked, resolved): {details}")
    return problems


def main():
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} <pip_dry_run_report.json> <lock_file.txt>")
    problems = compare(sys.argv[1], sys.argv[2])
    if problems:
        for p in problems:
            print(f"MISMATCH: {p}", file=sys.stderr)
        raise SystemExit(
            "native resolution does not exactly match the locked closure -- this is a scope "
            "change, not a routine difference; stopping per the execution packet's own rule")
    print("OK: native resolution exactly matches the locked closure.")


if __name__ == "__main__":
    main()
