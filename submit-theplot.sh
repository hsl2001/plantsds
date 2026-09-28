#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$SCRIPT_DIR}"
SELECTED_DIR="${1:-${SELECTED_DIR:-./selected}}"
SEG_BED="${2:-${SEG_BED:-./results/ANGIOSPERM.seg.bed}}"
OUTDIR="${OUTDIR:-./plots/gc_transition}"
BED_NAME="${SEG_BED##*/}"
BED_NAME="${BED_NAME%.seg.bed}"
OUT_PREFIX="${OUT_PREFIX:-$OUTDIR/${BED_NAME}_gc_transition}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-gc-transition}"
NCPUS="${NCPUS:-1}"
MEM="${MEM:-32gb}"
WALLTIME="${WALLTIME:-24:00:00}"

if [[ $# -gt 2 ]]; then
  printf 'Usage: %s [selected_dir] [input.seg.bed]\n' "$0" >&2
  exit 2
fi

cd "$WORKDIR"
if [[ ! -d "$SELECTED_DIR" ]]; then
  printf '[ERROR] selected directory not found: %s\n' "$SELECTED_DIR" >&2
  exit 1
fi
if [[ ! -f "$SELECTED_DIR/angio_wgd_genomes.files" && ! -f "$SELECTED_DIR/genomes.files" ]]; then
  printf '[ERROR] FASTA list not found in %s (expected angio_wgd_genomes.files or genomes.files)\n' \
    "$SELECTED_DIR" >&2
  exit 1
fi
if [[ ! -f "$SEG_BED" ]]; then
  printf '[ERROR] SegTrace BED not found: %s\n' "$SEG_BED" >&2
  exit 1
fi
if [[ ! -f "$WORKDIR/gc_transition_plot.py" ]]; then
  printf '[ERROR] worker script not found: %s/gc_transition_plot.py\n' "$WORKDIR" >&2
  exit 1
fi
if ! command -v qsub >/dev/null 2>&1; then
  printf '[ERROR] qsub is required to submit this PBS job.\n' >&2
  exit 1
fi

mkdir -p "$LOGDIR" "$OUTDIR"
printf '[gc-transition] selected=%s bed=%s output=%s\n' \
  "$SELECTED_DIR" "$SEG_BED" "$OUT_PREFIX"

qsub -N "$JOB_NAME" \
  -l "select=1:ncpus=${NCPUS}:mem=${MEM}" \
  -l "walltime=${WALLTIME}" \
  -v "WORKDIR=${WORKDIR},SELECTED_DIR=${SELECTED_DIR},SEG_BED=${SEG_BED},OUT_PREFIX=${OUT_PREFIX}" \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
uv run "$WORKDIR/gc_transition_plot.py" \
  --selected-dir "$SELECTED_DIR" \
  --bed "$SEG_BED" \
  --output-prefix "$OUT_PREFIX"
PBS
