#!/usr/bin/env bash
# Submit all standard SegTrace jobs except the angiosperm genus analysis.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

bash "$SCRIPT_DIR/submit-validate.sh"
bash "$SCRIPT_DIR/submit-parameter-sweep.sh"
bash "$SCRIPT_DIR/submit-simulation-sweep.sh"
bash "$SCRIPT_DIR/submit-sim-benchmark.sh"