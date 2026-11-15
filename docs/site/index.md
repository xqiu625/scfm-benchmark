# scfm-benchmark

**Donor-aware evaluation of frozen single-cell representations.**

The portable demo aligns cells by identity, keeps donor dependency groups intact,
fits training-fold probes, and produces donor-level predictions and a report.

```bash
git clone https://github.com/xqiu625/scfm-benchmark.git
cd scfm-benchmark
python -m pip install .
scfm-bench demo --out demo-output
```

Open `demo-output/results/report.html`. This is a CPU-only **synthetic** example,
not a reproduction of a biological result or a ranking of actual foundation models.
The development package is `1.2.0.dev0`; no new PyPI release is claimed.

- [CLI and input schema](cli.md)
- [Repository and example output](https://github.com/xqiu625/scfm-benchmark)
- [Current research status](https://github.com/xqiu625/scfm-benchmark/blob/main/docs/current-status.md)

The historical research pages are retained below. Their old headline numbers
predate the QC rebuild and must not be treated as current evidence. The later
scPRINT-2 biological-discovery work remains a separate research workflow.
