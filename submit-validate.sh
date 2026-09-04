#!/usr/bin/env bash
# Submit the SegTrace BLAST validation to PBS via `echo ... | qsub`.
# The heavy work (segtrace + makeblastdb + blastn + plotting) runs on the
# compute node, not the laptop. Compute logic lives in validate_segtrace.sh.
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-validate}"
THREADS="${THREADS:-128}"
MEM="${MEM:-480gb}"
WALLTIME="${WALLTIME:-96:00:00}"
# Optional command(s) to prepare the environment on the node, e.g.
# ENV_SETUP='micromamba activate py' (needs blastn, makeblastdb, python+numpy+matplotlib).
ENV_SETUP="${ENV_SETUP:-}"

mkdir -p "$LOGDIR"

echo "#!/usr/bin/env bash
set -euo pipefail
cd '$WORKDIR'
${ENV_SETUP:+$ENV_SETUP}
THREADS='$THREADS' ./validate_segtrace.sh" | qsub -N "$JOB_NAME" \
  -l nodes=node02:ppn=${THREADS} \
  -l mem=${MEM} \
  -l walltime=${WALLTIME} \
  -v WORKDIR="${WORKDIR}",THREADS="${THREADS}" \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log"
