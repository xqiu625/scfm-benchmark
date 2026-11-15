#!/usr/bin/env python3
"""99_audit — verify the project satisfies the audit-ready invariants.

Generates no scientific conclusions; only checks traceability/consistency.
Stdlib only so it can run anywhere. See README.md and architecture §12, §17.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import time
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=str, default=".")
    p.add_argument(
        "--out",
        type=str,
        default="scripts/99_audit/results/audit_summary.json",
    )
    return p.parse_args()


def setup_logging():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


# match backtick-quoted paths like `scripts/06_batch_generalization/`
PATH_RE = re.compile(r"`([A-Za-z0-9_./-]+)`")


def find_referenced_paths(md_path: Path, prefixes=("scripts/", "data/", "results/", "reports/")):
    if not md_path.exists():
        return []
    paths = []
    for m in PATH_RE.finditer(md_path.read_text()):
        cand = m.group(1).rstrip("/")
        if cand.startswith(prefixes):
            paths.append(cand)
    return sorted(set(paths))


def main():
    args = parse_args()
    setup_logging()
    start = time.time()
    root = Path(args.repo_root).resolve()
    logging.info("Auditing repo at %s", root)

    findings = {"critical": [], "major": [], "minor": [], "note": []}

    # 1. Pipeline-index paths exist
    idx = root / "docs/pipeline-index.md"
    missing = [p for p in find_referenced_paths(idx) if not (root / p).exists()]
    for p in missing:
        findings["major"].append(f"pipeline_index references missing path: {p}")

    # 2. Every numbered stage has a README and a manifest
    scripts_dir = root / "scripts"
    stage_dirs = sorted(
        d for d in scripts_dir.glob("[0-9]*_*") if d.is_dir()
    ) if scripts_dir.exists() else []
    for d in stage_dirs:
        if not (d / "README.md").exists():
            findings["minor"].append(f"stage missing README.md: {d.relative_to(root)}")
        manifest = d / "results" / "manifest.json"
        if not manifest.exists():
            findings["note"].append(
                f"stage has no results/manifest.json: {d.relative_to(root)}"
            )

    # 3. Claim/inventory evidence files referenced should exist
    for doc in ("docs/claim-traceability.md", "docs/figure-table-inventory.md"):
        for p in find_referenced_paths(root / doc):
            if "_pending_" in p:
                continue
            if not (root / p).exists():
                findings["note"].append(f"{doc} references not-yet-present path: {p}")

    summary = {
        "stage": "99_audit",
        "status": "success",
        "repo_root": str(root),
        "n_stages": len(stage_dirs),
        "counts": {k: len(v) for k, v in findings.items()},
        "findings": findings,
        "runtime_seconds": round(time.time() - start, 2),
    }

    out = root / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    logging.info(
        "Audit done: %d stages, %s. Wrote %s",
        len(stage_dirs),
        summary["counts"],
        out,
    )


if __name__ == "__main__":
    main()
