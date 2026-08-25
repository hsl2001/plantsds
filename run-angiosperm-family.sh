#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-angio-family}"
THREADS="${THREADS:-128}"
MEM="${MEM:-480gb}"
WALLTIME="${WALLTIME:-96:00:00}"

mkdir -p "$LOGDIR"
cd "$WORKDIR"
mkdir -p selected results

python3 ./select_angiosperm_family.py \
  --dataset-dir ./eukaryotic_data/ncbi_dataset/data \
  --report ./eukaryotic_data/ncbi_dataset/data/assembly_data_report.jsonl \
  --taxonomy-cache ./selected/taxonomy_rank_lineage.tsv \
  --out ./selected/angiosperm_family.files \
  --summary ./selected/angiosperm_family_min.tsv

qsub -N "$JOB_NAME" \
  -l nodes=node02:ppn=${THREADS} \
  -l mem=${MEM} \
  -l walltime=${WALLTIME} \
  -v WORKDIR="${WORKDIR}",THREADS="${THREADS}" \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
mkdir -p results

mapfile -t FASTAS < ./selected/angiosperm_family.files
printf '[segtrace] input FASTA count: %d\n' "${#FASTAS[@]}"

./time -v ./segtrace \
  -p "$THREADS" \
  -c 1 \
  -o ./results/ANGIOSPERM_FAMILY \
  "${FASTAS[@]}"
PBS
