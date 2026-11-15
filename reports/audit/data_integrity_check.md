# Data Integrity Check

Validate current data files against `docs/data-contracts.md`.

| Contract | File | Rows/Cells | Cols/Features | Unique IDs | Result |
|---|---|---|---|---|---|
| covid_sepsis_mortality | `data/processed/covid_sepsis_mortality.h5ad` | expect ~42,095 | TODO | donor_id | not run |

Run via `scripts/99_audit/code/check_data_contracts.py`.
