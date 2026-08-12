#!/usr/bin/env bash
# =============================================================================
# download_pangenomes.sh
# Multi-pangenome SD comparison: Download genome assemblies for 8 plant species
# =============================================================================
set -euo pipefail

BASE_DIR="$(pwd)/pangenome_data"
mkdir -p "$BASE_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

DRY_RUN=0

run_wget() {
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY RUN] wget $*"
    else
        wget "$@"
    fi
}

run_tar_if_exists() {
    local file="$1"
    shift
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY RUN] (if $file exists) tar $*"
    elif [ -f "$file" ]; then
        tar "$@" || true
    fi
}


# 2. WATERMELON (Citrullus lanatus) - CuGenDBv2
download_watermelon() {
    log "=== Downloading Watermelon Super-Pangenome ==="
    local DIR="$BASE_DIR/watermelon/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    run_wget -c -r -np -nH --cut-dirs=4 "http://cucurbitgenomics.org/v2/ftp/pan-genome/watermelon/graph_pangenome/assembly/" || true
}

# 3. TOMATO (Solanum lycopersicum) - Zenodo
download_tomato() {
    log "=== Downloading Tomato T2T Super-Pangenome ==="
    local DIR="$BASE_DIR/tomato/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    # Zenodo record 17878268 files are individual instead of a single tar.gz
    python3 -c "
import urllib.request, json
url = 'https://zenodo.org/api/records/17878268'
try:
    d = json.loads(urllib.request.urlopen(url).read())
    for f in d.get('files', []):
        print(f\"{f['links']['self']}\t{f['key']}\")
except Exception as e:
    pass
" | while read -r link key; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done
}

# 4. MARCHANTIA (Marchantia polymorpha) - MarpolBase
download_marchantia() {
    log "=== Downloading Marchantia Pangenome ==="
    local DIR="$BASE_DIR/marchantia/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    run_wget -c "https://marchantia.info/download/pangenome_assemblies.tar.gz" || true
    run_tar_if_exists "pangenome_assemblies.tar.gz" -xzf pangenome_assemblies.tar.gz
}

# 5. GRAPEVINE (Vitis vinifera) - Zenodo
download_grapevine() {
    log "=== Downloading Grapevine Pangenome ==="
    local DIR="$BASE_DIR/grapevine/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    for rec in 10851548 10846425; do
        curl -sL "https://zenodo.org/api/records/$rec" | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    for f in d.get('files', []):
        link = f.get('links', {}).get('self')
        key = f.get('key')
        if link and key and key.endswith('.fa.gz'):
            print(f'{link}\\t{key}')
except Exception:
    pass
" | while read -r link key; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done
    done
}

# 6. CITRUS (Citrus spp.) - HZAU FTP
download_citrus() {
    log "=== Downloading Citrus Pangenome ==="
    local DIR="$BASE_DIR/citrus/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    curl -s "http://citrus.hzau.edu.cn/download.php" | python3 -c "
import sys, re
try:
    html = sys.stdin.read()
    links = set(re.findall(r'href=\"(/data/Genome_info/[^\"]+)\"', html))
    for link in links:
        if link.endswith('.fa') or link.endswith('.fa.gz'):
            print(f'http://citrus.hzau.edu.cn{link}')
except Exception:
    pass
" | while read -r link; do
        if [ -n "$link" ]; then
            run_wget -c "$link" || true
        fi
    done
}

main() {
    for arg in "$@"; do
        if [[ "$arg" == "--dry-run" ]]; then
            DRY_RUN=1
            log "Running in DRY RUN mode"
        fi
    done

    log "Starting data downloads..."
    download_watermelon
    download_tomato
    download_marchantia
    download_grapevine
    download_citrus
    log "All downloads initiated/completed."
}

main "$@"
