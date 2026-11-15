# notebooks/

The runnable [platelet case demo](../examples/platelet/platelet_case_demo.ipynb)
lives with its guide in `examples/platelet/`. It presents seven models' archived
results, five stored UMAPs and an optional synthetic evaluation using the existing
portable functions. The real-data results and synthetic evaluation are separate.

Exploration only. Promote stable logic out of notebooks into `src/scfm_benchmark/`
and production `scripts/NN_*/code/run_*.py`.

Promotion rule (architecture §14):

```
notebook / explore script  ->  reusable function in src/  ->  production run_*.py
```

Do not let a notebook become the only record of an important result. Anything a
manuscript claim depends on must trace to a `scripts/` stage, not a notebook.
