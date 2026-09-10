#!/usr/bin/env bash
# Submit the default SegTrace simulation benchmark to PBSPro via qsub.
# The benchmark itself runs on the compute node and uses sim_benchmark.py's
# defaults for simulation and SegTrace options.
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-sim-benchmark}"
THREADS="${THREADS:-4}"
MEM="${MEM:-480gb}"
WALLTIME="${WALLTIME:-96:00:00}"
# Optional command(s) to prepare the environment on the node.
ENV_SETUP="${ENV_SETUP:-}"
OUT_DIR="${OUT_DIR:-sim_benchmark}"

mkdir -p "$LOGDIR"

qsub -N "$JOB_NAME" \
  -l select=1:ncpus=${THREADS} \
  -l walltime=${WALLTIME} \
  -v WORKDIR="${WORKDIR}",THREADS="${THREADS}",OUT_DIR="${OUT_DIR}",ENV_SETUP="${ENV_SETUP}" \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
${ENV_SETUP:+$ENV_SETUP}

if command -v micromamba >/dev/null 2>&1; then
	RUNNER=(micromamba run -n segtrace-compare uv run)
else
	RUNNER=(uv run)
fi

make clean && make
"${RUNNER[@]}" sim_benchmark.py --out-dir "$OUT_DIR" --force
PBS