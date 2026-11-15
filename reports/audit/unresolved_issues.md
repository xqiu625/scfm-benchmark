# Unresolved Issues

Open discrepancies found during review. Resolve or migrate to a stage README.

| ID | Severity | Issue | Owner | Status |
|---|---|---|---|---|
| A001 | note | No stage emits `results/manifest.json` yet (30 stages) | — | open |
| A002 | major | Claim traceability / figure inventory rows seeded but `needs-evidence` | — | open |
| A003 | major | HPCC working copy must be renamed `analyses/→scripts/` to match the mirror | — | open |
| A005 | major | 130 hardcoded `/bigdata/.../302006` paths in HPCC scripts; centralize into `configs/hpcc.yaml` (many are overridable argparse defaults) | — | open |
| A006 | note | `audit.py` stale-output majors (81) are mtime false-positives from the restructure; ignore | — | resolved |
