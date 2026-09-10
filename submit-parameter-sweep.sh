#!/usr/bin/env bash
# Submit the SegTrace k/w/s parameter sweep to PBSPro via qsub.
# The heavy work (make + genome generation + segtrace grid + evaluation) runs
# on the compute node, not the laptop. The sweep writes one CSV row per run
# (cmd as the last column).
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-sweep}"
THREADS="${THREADS:-128}"
MEM="${MEM:-480gb}"
WALLTIME="${WALLTIME:-96:00:00}"
QUEUE_NODE="${QUEUE_NODE:-node02}"
# Optional command(s) to prepare the environment on the node, e.g.
# ENV_SETUP='micromamba activate segtrace-compare' (needs uv + numpy).
ENV_SETUP="${ENV_SETUP:-}"

# --- Sweep grid (override via env) ---------------------------------------
export K_VALUES="${K_VALUES:-11 13 15 17 19 21}"
export W_VALUES="${W_VALUES:-256 512 1024 2048 4096 8192 9128}"
export S_VALUES="${S_VALUES:-1 2 4 8 16 32 64 100}"

# --- Benchmark genome (override via env) ---------------------------------
export OUT_DIR="${OUT_DIR:-sweep_benchmark}"
export OUT_CSV="${OUT_CSV:-parameter_sweep.csv}"
export SEGTRACE_BIN="${SEGTRACE_BIN:-./segtrace}"
export WORKDIR THREADS ENV_SETUP

mkdir -p "$LOGDIR"

qsub -N "$JOB_NAME" \
    -l select=1:ncpus=${THREADS}:mem=${MEM} \
  -l walltime=${WALLTIME} \
  -v WORKDIR="${WORKDIR}",THREADS="${THREADS}" \
  -V \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
export PATH="$HOME/.local/bin:$PATH"
${ENV_SETUP:+$ENV_SETUP}

if command -v micromamba >/dev/null 2>&1; then
	RUNNER=(micromamba run -n segtrace-compare uv run)
else
	RUNNER=(uv run)
fi

make clean && make

echo "=== Generating sweep benchmark genome in ${OUT_DIR} ==="
"${RUNNER[@]}" sim_benchmark.py \
	--out-dir "${OUT_DIR}" \
	--out-csv "${OUT_CSV}" \
	--force \
    --skip-tools

echo "=== Running k/w/s sweep ==="
"${RUNNER[@]}" --with numpy python - <<'PY'
import csv
import os
import sys
from pathlib import Path

import sim_benchmark as sb

out_dir = Path(os.environ["OUT_DIR"])
csv_path = out_dir / os.environ["OUT_CSV"]
segtrace_bin = os.environ["SEGTRACE_BIN"]
threads = os.environ["THREADS"]
k_values = os.environ["K_VALUES"].split()
w_values = os.environ["W_VALUES"].split()
s_values = os.environ["S_VALUES"].split()

truth_bed = out_dir / "truth.bed"
fastas = sorted(str(path) for path in (out_dir / "fasta").glob("*.fa"))
if not fastas:
    sys.exit(f"[ERROR] no FASTA files under {out_dir / 'fasta'}")

runs_dir = out_dir / "runs"
runs_dir.mkdir(exist_ok=True)

metric_order = [
    "status", "time_perf_seconds", "max_rss_kb",
    "pred_bp", "truth_bp", "intersect_bp", "bp_recall", "bp_precision", "bp_f1",
    "truth_fragments", "pred_fragments", "frag_tp", "frag_fp", "frag_fn",
    "frag_recall", "frag_precision", "frag_f1",
]
header = ["k", "w", "s", *metric_order, "cmd"]

total = len(k_values) * len(w_values) * len(s_values)
done = 0
with csv_path.open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(header)
    for k in k_values:
        for w in w_values:
            for s in s_values:
                label = f"k{k}_w{w}_s{s}"
                prefix = runs_dir / f"segtrace_{label}"
                command = [
                    segtrace_bin, "-k", k, "-s", s, "-w", w,
                    "-t", "0", "-c", "1", "-p", threads,
                    "-o", str(prefix), *fastas,
                ]
                cmd_str = " ".join(command)
                status = "ok"
                wall = 0.0
                rss = 0
                try:
                    profile = sb.profile_command(command, out_dir, label)
                    wall = profile.wall_seconds
                    rss = profile.max_rss_kb
                except Exception as exc:  # noqa: BLE001 - record failures as data
                    status = f"failed:{exc}"

                prediction = prefix.with_suffix(".seg.bed")
                metrics = sb.evaluate_with_numpy(prediction, truth_bed)
                row = {
                    "status": status,
                    "time_perf_seconds": f"{wall:.6f}",
                    "max_rss_kb": rss,
                    **metrics,
                }
                writer.writerow([k, w, s, *(row.get(name, "") for name in metric_order), cmd_str])

                if prediction.exists():
                    prediction.unlink()
                done += 1
                print(f"[{done}/{total}] {label} frag_f1={metrics.get('frag_f1', 0.0):.4f}", flush=True)

print(f"[INFO] sweep CSV: {csv_path}")
PY

echo "=== Sweep complete: ${OUT_DIR}/${OUT_CSV} ==="
PBS
