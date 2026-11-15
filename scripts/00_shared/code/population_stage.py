"""Machine-checkable population-stage identifiers, kept separate from the free-text
descriptions already used elsewhere (e.g. `uns["population_stage"]`'s prose).

A 2026-10-09 review (review:20261009_230728) found that the preparer's
free-text `uns["population_stage"]` ("pre-Gate-2 PROVISIONAL, not final (see
docs/current-status.md)") and `PROVISIONAL_JOB_SPEC_v2.md`'s documented
`--population-stage` argument ("pre-Gate-2 decontx_input PROVISIONAL, not final
(2026-10-09)") were different strings for the same stage -- so following the
documented command verbatim triggered the exact mismatch the runner's check was
supposed to catch, and `run_batch_diagnostics.py` only WARNED on any mismatch, never
failed. This module gives both tools one small, fixed vocabulary of stage IDs to
compare instead of free text, and makes a mismatch a hard failure -- mislabeling
provisional data as final-stage (or vice versa) is exactly what this exists to
prevent, so it must not be possible to merely log past it.
"""
from __future__ import annotations

KNOWN_STAGE_IDS = frozenset({
    "synthetic_fixture",
    "pre_gate2_decontx_provisional",
})


def check_stage_id(declared: str, actual: str, *, context: str) -> None:
    """Raises ValueError (never just warns) on an unregistered or mismatched stage
    id. `declared` is what the caller (e.g. a CLI flag) asserts; `actual` is what the
    input data itself records (e.g. `adata.uns["population_stage_id"]`)."""
    if declared not in KNOWN_STAGE_IDS:
        raise ValueError(
            f"{context}: stage id {declared!r} is not a registered KNOWN_STAGE_IDS "
            f"value ({sorted(KNOWN_STAGE_IDS)}) -- refusing an unregistered stage claim."
        )
    if declared != actual:
        raise ValueError(
            f"{context}: declared stage id {declared!r} does not match the input's "
            f"own recorded stage id {actual!r} -- refusing to proceed with a "
            f"mismatched stage claim. This is a hard failure, not a warning: "
            f"mislabeling provisional data as a different stage (or vice versa) is "
            f"exactly what this check exists to prevent."
        )
