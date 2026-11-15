"""Tests for stable_seed.derive_seed() -- the fix for the hash()-randomization
bug a 2026-10-09 review found in both run_batch_diagnostics.py and
prepare_provisional_diagnostic_input.py's patient-balanced samplers."""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/00_shared/code"))
from stable_seed import derive_seed


def test_deterministic_within_one_process():
    assert derive_seed(0, "patient", "P001") == derive_seed(0, "patient", "P001")


def test_different_parts_give_different_seeds():
    a = derive_seed(0, "patient", "P001")
    b = derive_seed(0, "patient", "P002")
    assert a != b


def test_no_part_boundary_collision():
    """Naive string concatenation could collide ("ab","c") with ("a","bc");
    the NUL-byte separator must prevent that."""
    a = derive_seed("ab", "c")
    b = derive_seed("a", "bc")
    assert a != b


def test_stable_across_fresh_python_processes_with_different_pythonhashseed():
    """The actual bug: hash() is salted per-process. derive_seed() must give
    the SAME value in two freshly spawned interpreters with DIFFERENT
    PYTHONHASHSEED, which hash() would not."""
    script = (
        "import sys; sys.path.insert(0, %r); "
        "from stable_seed import derive_seed; "
        "print(derive_seed(0, 'patient', 'P001'))"
        % str(Path(__file__).parent.parent / "scripts/00_shared/code")
    )
    out1 = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env={"PYTHONHASHSEED": "1", "PATH": __import__("os").environ.get("PATH", "")},
    )
    out2 = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env={"PYTHONHASHSEED": "42", "PATH": __import__("os").environ.get("PATH", "")},
    )
    assert out1.returncode == 0, out1.stderr
    assert out2.returncode == 0, out2.stderr
    assert out1.stdout.strip() == out2.stdout.strip(), (
        f"derive_seed gave different values under different PYTHONHASHSEED: "
        f"{out1.stdout.strip()!r} vs {out2.stdout.strip()!r}"
    )


def test_builtin_hash_is_the_bug_being_fixed():
    """Negative control: confirms the ORIGINAL bug (hash() on a tuple
    containing a str) really does change across PYTHONHASHSEED values in
    this environment -- i.e. this test suite's premise is real, not assumed."""
    script = "print(hash((0, 'patient', 'P001')) % (2**32 - 1))"
    out1 = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env={"PYTHONHASHSEED": "1", "PATH": __import__("os").environ.get("PATH", "")},
    )
    out2 = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env={"PYTHONHASHSEED": "42", "PATH": __import__("os").environ.get("PATH", "")},
    )
    assert out1.stdout.strip() != out2.stdout.strip(), (
        "expected Python's built-in hash() on a str-containing tuple to differ "
        "across PYTHONHASHSEED values -- if it no longer does in this Python "
        "version, the bug this module fixes may not apply the way documented."
    )
