#!/usr/bin/env bash
# Run on HPCC with a new output path; no SLURM, real data or shared-env changes.
set -euo pipefail
repo_root=$(cd "$(dirname "$0")/../.." && pwd)
output=${1:?Usage: verify_portable_package.sh NEW_OUTPUT_DIRECTORY}
base_python=${SCFM_TEST_PYTHON:-/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python}
test ! -e "$output"
mkdir -p "$output"
output=$(cd "$output" && pwd)
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH
"$base_python" -m venv --system-site-packages "$output/venv"
package_python="$output/venv/bin/python"
"$package_python" -m pip install --disable-pip-version-check build hatchling 'ruff>=0.5,<0.10' 'black>=24,<25'
"$package_python" -m black --check "$repo_root/tests/test_portable_benchmark.py"
"$package_python" -m ruff check "$repo_root/src/scfm_benchmark/benchmark.py" "$repo_root/src/scfm_benchmark/demo.py" "$repo_root/tests/test_portable_benchmark.py"
"$package_python" -m build --no-isolation --outdir "$output/dist" "$repo_root"
"$package_python" -m pip install --no-deps --ignore-installed "$output"/dist/*.whl
mkdir "$output/outside-checkout"
cd "$output/outside-checkout"
"$package_python" -c 'import scfm_benchmark; from pathlib import Path; import sys; print(scfm_benchmark.__file__); assert Path(scfm_benchmark.__file__).is_relative_to(Path(sys.prefix))'
"$package_python" -m pytest "$repo_root/tests/test_portable_benchmark.py" --junitxml="$output/pytest.xml"
"$output/venv/bin/scfm-bench" demo --out "$output/demo"
"$output/venv/bin/scfm-bench" benchmark --metadata "$output/demo/inputs/metadata.csv" \
  --embedding "toy_signal=$output/demo/inputs/toy_signal.csv" \
  --embedding "toy_noise=$output/demo/inputs/toy_noise.csv" \
  --embedding "toy_expression_pca=$output/demo/inputs/toy_expression_pca.csv" \
  --pca toy_expression_pca=5 --out "$output/csv-roundtrip"
cmp "$output/demo/results/metrics.csv" "$output/csv-roundtrip/metrics.csv"
cmp "$output/demo/results/donor_predictions.csv" "$output/csv-roundtrip/donor_predictions.csv"
"$package_python" -m pip freeze > "$output/environment.txt"
printf 'Portable wheel, tests, demo and CSV round trip passed: %s\n' "$output"
