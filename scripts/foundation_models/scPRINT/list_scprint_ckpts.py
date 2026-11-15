#!/usr/bin/env python3
"""List downloaded scPRINT checkpoints with sizes.

Run on HPCC after `hf download jkobject/scPRINT`:
    python scripts/foundation_models/scPRINT/list_scprint_ckpts.py
"""
import os
import glob

HF_HOME = os.environ.get("HF_HOME", "/bigdata/godziklab/shared/Xinru/.huggingface")
SNAPSHOT_GLOB = f"{HF_HOME}/hub/models--jkobject--scPRINT/snapshots/*"

PATTERNS = ("*.ckpt", "*.safetensors", "*.bin", "*.pt")


def main() -> None:
    snapshots = sorted(glob.glob(SNAPSHOT_GLOB))
    if not snapshots:
        raise SystemExit(f"No scPRINT snapshots found under {SNAPSHOT_GLOB}")

    for snap in snapshots:
        print(f"snapshot: {snap}")
        files = []
        for pat in PATTERNS:
            files.extend(glob.glob(os.path.join(snap, pat)))
            files.extend(glob.glob(os.path.join(snap, "**", pat), recursive=True))
        files = sorted(set(files))
        if not files:
            print("  (no weight files matched)")
            continue
        for f in files:
            size_gb = os.path.getsize(f) / 1e9
            print(f"  {size_gb:6.2f} GB  {os.path.relpath(f, snap)}")


if __name__ == "__main__":
    main()
