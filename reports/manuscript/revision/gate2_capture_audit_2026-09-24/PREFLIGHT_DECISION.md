# Technical preflight decision

Date: 2026-09-24

Decision: **approved for one score-free technical preflight only**, after deployment from the committed clean worktree.

The captured test transcript reports 205 passing focused/repository tests and 12/12 killed mutations for the numerical preflight implementation. During final static audit, an automated review found that `JOB_PACKAGE.md` specified 48 GB but the wrapper still requested 24 GB. The resulting change affected only the SLURM resource directive from `--mem=24G` to `--mem=48G`; no Python, parameter, input, or output semantics changed after the recorded test run. The committed wrapper hash therefore intentionally differs from the pre-correction hash in `TEST_TRANSCRIPT.txt`.

Authorization is limited to:

- capture `GSM5102900_HC1`;
- the committed `run_capture_preflight.sh` with `--mem=48G`;
- one job in a clean deployed worktree;
- score-free feasibility/manifold measurement only;
- the new versioned output directory constructed by the wrapper.

It does not authorize doublet scoring, calls, rates, outcome access, production-ruling activation, cell deletion, a retry with altered code/flags, or any additional capture. The neighbour graph remains unchecked by design and must remain an explicit limitation in the later config audit.
