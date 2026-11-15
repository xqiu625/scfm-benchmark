"""Runs the actual runnable example end-to-end, so it stays exercised by the
regular test suite rather than only being a one-off manual demo. Per
`review:20261012_163250`'s "交付缺口": a config template + fixture +
the full validate -> plan -> read_manifest chain, not just a test fixture.
"""

import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(
    REPO_ROOT, "scripts", "34_scfm_discovery_framework", "examples", "run_example.py"
)


def test_run_example_end_to_end(tmp_path):
    out_dir = tmp_path / "example_run"
    result = subprocess.run(
        [sys.executable, EXAMPLE, "--out-dir", str(out_dir)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "OK: read_manifest" in result.stdout
    assert (out_dir / "data" / "dataset_config.json").exists()
    assert (out_dir / "data" / "tiny_monocyte_example.h5ad").exists()
    assert (out_dir / "plan" / "PLAN.json").exists()


def test_run_example_refuses_nonempty_out_dir(tmp_path):
    out_dir = tmp_path / "example_run"
    out_dir.mkdir()
    (out_dir / "junk.txt").write_text("pre-existing")
    result = subprocess.run(
        [sys.executable, EXAMPLE, "--out-dir", str(out_dir)], capture_output=True, text=True
    )
    assert result.returncode != 0
    assert (out_dir / "junk.txt").read_text() == "pre-existing"
