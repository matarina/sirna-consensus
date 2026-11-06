#!/usr/bin/env bash
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$package_dir"
export PYTHONPATH="$package_dir/src"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
model_python="${SIRNA_MODEL_PYTHON:-$package_dir/.scorer-venv/bin/python}"
"$model_python" -m sirna_consensus.cli train --project .. \
  --hgnc references/gencode_v50/hgnc_complete_set.txt --out models/development_20261106
"$model_python" benchmarks/scripts/verify_platform_v3.py --model-dir models/development_20261106
"$model_python" benchmarks/scripts/evaluate_platform_v3.py --model-dir models/development_20261106 \
  --out benchmarks/results/platform_v3
bowtie_executable="${SIRNA_BOWTIE:-$package_dir/tools/bowtie-1.3.1/bowtie-1.3.1-linux-x86_64/bowtie}"
"$model_python" benchmarks/scripts/validate_human_screen.py \
  --manifest references/gencode_v50/built/screen_manifest.json \
  --bowtie "$bowtie_executable" --out benchmarks/results/platform_v3 \
  --model-dir models/development_20261106
"$model_python" benchmarks/scripts/build_platform_report_zh.py
