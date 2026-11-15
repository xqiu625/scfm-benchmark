"""Per review:20261015_150231 F1: "真实测试不能以 skip 算成功" -- a bare
`pytest` exit code of 0, or a log file ending in some success-looking line, does not
distinguish a genuinely PASSED test from one that was SKIPPED (e.g. because pydeseq2
still is not importable). This parses pytest's own machine-readable JUnit XML output
and refuses unless EXACTLY ONE testcase matches the expected name and its outcome is
an actual pass -- no `<skipped>`, `<failure>`, or `<error>` child, and not simply
"found zero matching testcases" (which would otherwise look like success by omission).

Usage: python de_backend_verify_synthetic_test_passed.py <junit.xml> <test_name_substring>
Exit 0 only if exactly one matching, genuinely-passed testcase is found. Exit 1 with
the reason otherwise (not found, more than one match, skipped, failed, errored).
"""
import sys
import xml.etree.ElementTree as ET


def check(xml_path, test_name_substring):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    matches = [
        tc for tc in root.iter("testcase")
        if test_name_substring in tc.get("name", "") or test_name_substring in tc.get("classname", "")
    ]
    if len(matches) != 1:
        return [f"expected exactly 1 testcase matching {test_name_substring!r}, found {len(matches)}"]

    (tc,) = matches
    problems = []
    skipped = tc.find("skipped")
    if skipped is not None:
        problems.append(f"test was SKIPPED, not passed: {skipped.get('message', '')!r}")
    failure = tc.find("failure")
    if failure is not None:
        problems.append(f"test FAILED: {failure.get('message', '')!r}")
    error = tc.find("error")
    if error is not None:
        problems.append(f"test ERRORED: {error.get('message', '')!r}")
    return problems


def main():
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} <junit.xml> <test_name_substring>")
    problems = check(sys.argv[1], sys.argv[2])
    if problems:
        for p in problems:
            print(f"NOT A GENUINE PASS: {p}", file=sys.stderr)
        raise SystemExit("the synthetic backend verification test did not genuinely pass")
    print(f"OK: exactly one testcase matching {sys.argv[2]!r} genuinely passed.")


if __name__ == "__main__":
    main()
