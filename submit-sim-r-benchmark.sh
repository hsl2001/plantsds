#!/usr/bin/env bash
# Submit the Badread-based SegTrace -r benchmark to PBSPro.
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-sim-r-benchmark}"
THREADS="${THREADS:-16}"
MEM="${MEM:-100gb}"
WALLTIME="${WALLTIME:-96:00:00}"
ENV_SETUP="${ENV_SETUP:-}"
OUT_DIR="${OUT_DIR:-sim_r_benchmark}"
SEED="${SEED:-42}"
QUANTITIES="${QUANTITIES:-1x:2x:4x:8x:16x}"
MIN_REPORT_COPIES="${MIN_REPORT_COPIES:-2}"

mkdir -p "$LOGDIR"

qsub -N "$JOB_NAME" \
  -l select=1:ncpus=${THREADS}:mem=${MEM} \
  -l walltime=${WALLTIME} \
  -v WORKDIR="${WORKDIR}",THREADS="${THREADS}",OUT_DIR="${OUT_DIR}",SEED="${SEED}",QUANTITIES="${QUANTITIES}",MIN_REPORT_COPIES="${MIN_REPORT_COPIES}",ENV_SETUP="${ENV_SETUP}" \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
export PATH="$HOME/.local/bin:$PATH"
${ENV_SETUP:+$ENV_SETUP}

if command -v micromamba >/dev/null 2>&1; then
  ENV_NAME="benchmark-segtrace"
  if ! micromamba run -n "$ENV_NAME" python -c 'import edlib' >/dev/null 2>&1; then
    echo "[INFO] Installing edlib in micromamba environment: $ENV_NAME"
    micromamba install -n "$ENV_NAME" -y -c bioconda -c conda-forge edlib
  fi
  micromamba run -n "$ENV_NAME" python -c \
    'import edlib, sys; print(f"[INFO] Python: {sys.executable}"); print(f"[INFO] edlib: {edlib.__file__}")'
  RUNNER=(micromamba run -n benchmark-segtrace uv run)
else
  python3 -c 'import edlib, sys; print(f"[INFO] Python: {sys.executable}"); print(f"[INFO] edlib: {edlib.__file__}")'
  RUNNER=(uv run)
fi

command -v badread >/dev/null 2>&1

make clean && make

echo "=== Running Badread SegTrace -r benchmark ==="
"${RUNNER[@]}" sim_r_benchmark.py \
  --out-dir "$OUT_DIR" \
  --force \
  --threads "$THREADS" \
  --seed "$SEED" \
  --quantities "$QUANTITIES" \
  --min-report-copies "$MIN_REPORT_COPIES"
PBS
