#!/usr/bin/env python3
"""Deterministic seed derivation that does NOT use Python's built-in
`hash()`. `hash()` on str/tuple-of-str is salted per-process (PYTHONHASHSEED)
unless explicitly disabled, so `hash((seed, patient_id))` gives a DIFFERENT
value in every fresh Python process -- a fixed `--seed` CLI argument does not
actually reproduce the same patient-balanced selection across runs. Found by
a 2026-10-09 review of `run_batch_diagnostics.py` and
`prepare_provisional_diagnostic_input.py`, both of which used exactly this
pattern despite the synthetic-fixture generator having already been fixed
for the same defect in its outcome assignment.

Usage:
    from stable_seed import derive_seed
    rng = np.random.RandomState(derive_seed(seed, "patient", patient_id))
"""
from __future__ import annotations

import hashlib


def derive_seed(*parts) -> int:
    """Deterministic, cross-process-stable 32-bit seed from arbitrary parts.
    Canonicalizes each part to str and joins with a separator that cannot
    appear inside a part's own repr by using a NUL byte (never valid in these
    identifiers) -- collisions would otherwise be possible between e.g.
    parts=("ab", "c") and parts=("a", "bc") under naive concatenation."""
    joined = "\x00".join(str(p) for p in parts)
    digest = hashlib.sha256(joined.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], byteorder="big")
