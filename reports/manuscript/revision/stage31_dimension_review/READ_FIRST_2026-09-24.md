# READ THIS FIRST — Stage 31 dimension bundle

**Before opening `README.md`, read `ADDENDUM_2026-09-24.md`.**

The bundle's own `README.md` is frozen at 2026-09-22 and asks you to select a fixed
`n_prin_comps`. **That is no longer the question.** A ruling recorded on 2026-09-23 stopped
dimension selection for the design this bundle describes.

The addendum states what is actually being asked: whether Scrublet should run on the complete,
pre-subsample per-capture/sample raw-count input, with calls and status mapped onto the 60-cell
analytic subset afterward — and six specific points to rule on, each with an explicit
"insufficient evidence, do not proceed" option.

The seven original files are unchanged and remain correct as computation. They are now a
**historical feasibility attachment**, not a pending parameter choice.

| file | role |
|---|---|
| `ADDENDUM_2026-09-24.md` | **the current review task** |
| `README.md` | frozen 2026-09-22 cover; superseded task statement |
| `scrublet_dimension_parameter.md` | the open-decision record |
| `RUN_MANIFEST.json` | run identity for the audit |
| `dimension_audit.csv` | 301 measured samples |
| `dimension_coverage_table.csv` | 59 integer candidates |
| `dimension_excluded_samples_long.csv` | 8,247 exclusion records |
| `build_coverage_table.py` | the script that derives the two tables |

`SHA256SUMS.txt` covers those eight content files. This cover note and the checksum file are not
listed in it.

**External human review of Stage 31 has not happened.** The 2026-09-23 ruling is an internal
decision record, not an external review.
