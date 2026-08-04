#!/usr/bin/env bash
# =============================================================================
# download_pangenomes.sh
# Multi-pangenome SD comparison: Download genome assemblies for 8 plant species
# =============================================================================
set -euo pipefail

BASE_DIR="$(pwd)/pangenome_data"
mkdir -p "$BASE_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }


# 2. WATERMELON (Citrullus lanatus) - CuGenDBv2
download_watermelon() {
    log "=== Downloading Watermelon Super-Pangenome ==="
    local DIR="$BASE_DIR/watermelon/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    wget -c -r -np -nH --cut-dirs=4 "http://cucurbitgenomics.org/v2/ftp/pan-genome/watermelon/graph_pangenome/assembly/" || true
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
            wget -c "$link" -O "$key" || true
        fi
    done
}

# 4. MARCHANTIA (Marchantia polymorpha) - MarpolBase
download_marchantia() {
    log "=== Downloading Marchantia Pangenome ==="
    local DIR="$BASE_DIR/marchantia/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    wget -c "https://marchantia.info/download/pangenome_assemblies.tar.gz" || true
    [ -f "pangenome_assemblies.tar.gz" ] && tar -xzf pangenome_assemblies.tar.gz
}

# 5. GRAPEVINE (Vitis vinifera) - Zenodo
download_grapevine() {
    log "=== Downloading Grapevine Pangenome ==="
    local DIR="$BASE_DIR/grapevine/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    wget -c "https://zenodo.org/records/10851547/files/Grapepan_v1.0.tar.gz" || true
    wget -c "https://zenodo.org/records/10846425/files/T2T_genomes.tar.gz" || true
    [ -f "Grapepan_v1.0.tar.gz" ] && tar -xzf Grapepan_v1.0.tar.gz
    [ -f "T2T_genomes.tar.gz" ] && tar -xzf T2T_genomes.tar.gz
}

# 6. CITRUS (Citrus spp.) - HZAU FTP
download_citrus() {
    log "=== Downloading Citrus Pangenome ==="
    local DIR="$BASE_DIR/citrus/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    wget -c -r -np -A "*.fa.gz" "http://citrus.hzau.edu.cn/download/assemblies/" || true
}

# 7. ARABIDOPSIS (Arabidopsis thaliana) - GitHub / Edmond
download_arabidopsis() {
    log "=== Downloading Arabidopsis 69 Pangenome ==="
    local DIR="$BASE_DIR/arabidopsis/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    # Download from Edmond Dataverse 10.17617/3.AEOJBL
    python3 -c "
import urllib.request, json
url = 'https://edmond.mpdl.mpg.de/api/datasets/:persistentId/?persistentId=doi:10.17617/3.AEOJBL'
try:
    d = json.loads(urllib.request.urlopen(url).read())
    for f in d['data']['latestVersion']['files']:
        print(f\"https://edmond.mpdl.mpg.de/api/access/datafile/{f['dataFile']['id']}\t{f['dataFile']['filename']}\")
except Exception:
    pass
" | while read -r link key; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            wget -c "$link" -O "$key" || true
        fi
    done
}

# 8. RICE (Oryza sativa) - Figshare / ENA
download_rice() {
    log "=== Downloading Rice 149 Pangenome ==="
    local DIR="$BASE_DIR/rice/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    # Download from Figshare API 25697817
    python3 -c "
import urllib.request, json
url = 'https://api.figshare.com/v2/articles/25697817/files'
try:
    d = json.loads(urllib.request.urlopen(url).read())
    for f in d:
        print(f\"{f['download_url']}\t{f['name']}\")
except Exception:
    pass
" | while read -r link key; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            wget -c "$link" -O "$key" || true
        fi
    done
}

main() {
    log "Starting data downloads..."
    download_watermelon
    download_tomato
    download_marchantia
    download_grapevine
    download_citrus
    download_arabidopsis
    download_rice
    log "All downloads initiated/completed."
}

main "$@"
