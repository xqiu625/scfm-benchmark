"""Atomic ownership of an output path.

A 2026-10-09 review (review:20261009_230728) found both stage-33 CLIs used a
check-then-act output guard: check the destination doesn't exist / is empty, THEN
create it. That is not atomic -- two processes launched close together can both pass
the check before either has written anything, and then race on the same output. This
is exactly what happened in this project's own session history: an orphaned
`nohup`'d process (never actually killed, just invisible to a `ps` check that landed
on the other login node in this cluster's round-robin) raced with a later "clean"
retry, and `save_representation()`'s write-then-reload-verify step caught the other
process mid-write (`[Errno 2] No such file or directory` for a file that existed
moments later). The existing no-clobber check turned that into a diagnosable
`FileExistsError` rather than silent corruption -- but it did not PREVENT the race;
it only made one specific symptom of it visible.

`claim_ownership()` closes the gap with an OS-guaranteed-atomic exclusive create
(`os.O_CREAT | os.O_EXCL`): the loser gets `FileExistsError` directly from the
syscall, never from a race window a human has to reconstruct after the fact. The
ownership record itself carries host, PID, start time and the exact command, so a
later liveness check can be aimed at the correct host -- this cluster's process
table is NOT shared between its two login nodes, so a `ps` check on the wrong host
looks identical to "the process is gone" (this is the same mistake the race above
was originally misdiagnosed as).
"""
from __future__ import annotations

import json
import os
import socket
from datetime import datetime, timezone
from pathlib import Path


class OwnershipError(FileExistsError):
    """Raised when claim_ownership() loses the race; carries the loser's lock path
    and (if readable) the winner's recorded ownership payload for diagnosis."""


def claim_ownership(lock_path: Path, *, command: str) -> dict:
    """Atomically creates `lock_path` and writes a JSON ownership record into it, or
    raises OwnershipError if another process already holds it. Never deletes or
    overwrites an existing lock -- a lock left behind by a genuinely crashed run is
    evidence, not clutter, per this project's standing "never overwrite a result"
    rule; a caller that hits this must pick a fresh output path, same as any other
    already-exists refusal in this project, not force past it."""
    payload = dict(
        host=socket.gethostname(),
        pid=os.getpid(),
        started_at=datetime.now(timezone.utc).isoformat(),
        command=command,
        run_path=str(lock_path.parent),
    )
    try:
        fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        existing = None
        try:
            existing = json.loads(lock_path.read_text())
        except Exception:
            pass
        raise OwnershipError(
            f"{lock_path} is already owned by another run -- refusing to write into "
            f"this output path. Recorded owner: {existing!r}. If that process is "
            f"confirmed not running (check the SAME host named above -- this "
            f"cluster's process table is not shared across login nodes), use a "
            f"fresh output path rather than deleting this lock; the lock itself is "
            f"part of the run record."
        )
    try:
        os.write(fd, json.dumps(payload, indent=2).encode())
    finally:
        os.close(fd)
    return payload
