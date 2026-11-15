# 99_audit

## Purpose
Verify that existing conclusions are traceable, reproducible, and internally
consistent. This stage generates **no scientific conclusions** — it only checks
that the project satisfies the audit-ready invariants (architecture §12, §17).

## Inputs
| Input | Path | Produced by | Required |
|---|---|---|---|
| Pipeline index | `docs/pipeline-index.md` | maintainers | stage path table |
| Stage manifests | `scripts/NN_*/results/manifest.json` | each stage | optional (warns if absent) |
| Claim traceability | `docs/claim-traceability.md` | maintainers | evidence file column |
| Figure/table inventory | `docs/figure-table-inventory.md` | maintainers | source/script columns |

## Outputs
| Output | Path | Description |
|---|---|---|
| Audit summary | `scripts/99_audit/results/audit_summary.json` | machine-readable check results |
| Audit log | `scripts/99_audit/results/run_log.txt` | full run log |

## Environment
Python ≥3.10, stdlib only (no heavy deps). Runs locally or on HPCC.

## Run
```bash
# local
python scripts/99_audit/code/run_audit.py --repo-root .
# HPCC
sbatch scripts/99_audit/submit.sh
```

## Checks
- All stage paths in `docs/pipeline-index.md` exist.
- Every numbered stage has a `README.md`.
- Stage manifests present and report success (warn if missing).
- Claim/inventory evidence files referenced actually exist on disk.
- Flag orphan results, stale outputs (output older than generating code/input).

## Known Issues
- Manifests are not yet emitted by most stages (open item A001).
- Claim/inventory tables are seeded but `needs-evidence`.

## Audit Trail
Supports the historical
`archive/superseded_2026-09-21/docs_legacy/audit_report_2026-06.md` and the
`reports/audit/*` dimension reports.
