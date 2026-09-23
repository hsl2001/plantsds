#!/usr/bin/env bash
# Resumable downloader for NCBI genome FASTAs that are absent locally.
# Places each accession under $DATA_DIR/<accession>/ so build_assembly_index()
# in select_angiosperm_genus.py can find it via "*/*_genomic.fna*".
set -euo pipefail

DATA_DIR="${DATA_DIR:-./eukaryotic_data/ncbi_dataset/data}"
API_BASE="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"

# Accessions to fetch: CLI args, else $MISSING_FILE (one per line), else the
# nine that were reported missing by select_angiosperm_genus.py.
DEFAULT_MISSING=(
  GCA_039602395.1 GCA_041429995.1 GCA_965282405.1 GCA_963692815.1
  GCA_965197055.1 GCF_003573695.1 GCA_053477545.1 GCA_977020475.1
  GCA_048183575.1
)

if [[ $# -gt 0 ]]; then
  ACCESSIONS=("$@")
elif [[ -n "${MISSING_FILE:-}" && -f "${MISSING_FILE}" ]]; then
  mapfile -t ACCESSIONS < "${MISSING_FILE}"
else
  ACCESSIONS=("${DEFAULT_MISSING[@]}")
fi

log() { printf '[download] %s\n' "$*"; }

has_fasta() {
  compgen -G "${DATA_DIR}/$1/"'*_genomic.fna*' > /dev/null 2>&1
}

download_via_cli() {
  local acc="$1" zip="$2"
  datasets download genome accession "$acc" \
    --include genome --no-progressbar --filename "$zip"
}

download_via_api() {
  local acc="$1" zip="$2"
  curl -fSL --retry 3 --retry-delay 5 \
    "${API_BASE}/${acc}/download?include_annotation_type=GENOME_FASTA" \
    -o "$zip"
}

fetch_one() {
  local acc="$1"
  if has_fasta "$acc"; then
    log "skip ${acc} (already present)"
    return 0
  fi

  local tmp zip
  tmp="$(mktemp -d)"
  zip="${tmp}/${acc}.zip"
  trap 'rm -rf "$tmp"' RETURN

  log "fetching ${acc} ..."
  if command -v datasets > /dev/null 2>&1; then
    download_via_cli "$acc" "$zip"
  else
    download_via_api "$acc" "$zip"
  fi

  unzip -q -o "$zip" -d "$tmp"

  local src="${tmp}/ncbi_dataset/data/${acc}"
  if [[ ! -d "$src" ]] || ! compgen -G "${src}/"'*_genomic.fna*' > /dev/null 2>&1; then
    log "ERROR: no genomic FASTA in download for ${acc}" >&2
    return 1
  fi

  mkdir -p "${DATA_DIR}/${acc}"
  cp -f "${src}/"*_genomic.fna* "${DATA_DIR}/${acc}/"
  log "done ${acc}"
}

command -v unzip > /dev/null 2>&1 || { echo "unzip is required" >&2; exit 1; }

mkdir -p "$DATA_DIR"
failed=()
for acc in "${ACCESSIONS[@]}"; do
  [[ -z "$acc" ]] && continue
  if ! fetch_one "$acc"; then
    failed+=("$acc")
  fi
done

if [[ ${#failed[@]} -gt 0 ]]; then
  log "FAILED: ${failed[*]}" >&2
  exit 1
fi
log "all accessions present"
