#!/usr/bin/env bash
set -e
export PATH="$HOME/.local/bin:$PATH"

if command -v micromamba >/dev/null 2>&1; then
	RUNNER=(micromamba run -n plantsds-bench uv run)
else
	RUNNER=(uv run)
fi

make clean && make
echo "=== Running SegTrace simulation smoke benchmark ==="
"${RUNNER[@]}" sim_benchmark.py \
        --out-dir benchmark \
        --force \
        --species 3 \
        --chromosomes 2 \
        --chrom-length 2000000 \
        --fragments 10 \
        --threads 8