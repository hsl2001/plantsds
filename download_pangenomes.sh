#!/usr/bin/env bash
# =============================================================================
# download_pangenomes.sh
# Multi-pangenome SD comparison: Download genome assemblies for 8 plant species
# =============================================================================
set -euo pipefail

BASE_DIR="$(pwd)/pangenome_data"
mkdir -p "$BASE_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# 1. CUCUMBER (Cucumis sativus) - NGDC
download_cucumber() {
    log "=== Downloading Cucumber Pangenome ==="
    local DIR="$BASE_DIR/cucumber/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    # Requires NGDC GSA download links for PRJCA038097, PRJCA043228, PRJCA038675
    # wget -q "https://ngdc.cncb.ac.cn/gsa/browse/PRJCA038097/download"
    log "Cucumber: Check NGDC PRJCA038097 for exact fasta URLs."
}

# 2. WATERMELON (Citrullus lanatus) - CuGenDBv2
download_watermelon() {
    log "=== Downloading Watermelon Super-Pangenome ==="
    local DIR="$BASE_DIR/watermelon/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    wget -c -r -np -nH --cut-dirs=4 "http://cucurbitgenomics.org/v2/ftp/pan-genome/watermelon/graph_pangenome/assemblies/" || true
}

# 3. TOMATO (Solanum lycopersicum) - Zenodo
download_tomato() {
    log "=== Downloading Tomato T2T Super-Pangenome ==="
    local DIR="$BASE_DIR/tomato/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    wget -c "https://zenodo.org/records/17878268/files/Tomato_T2T_assemblies.tar.gz" || true
    [ -f "Tomato_T2T_assemblies.tar.gz" ] && tar -xzf Tomato_T2T_assemblies.tar.gz
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

# 7. ARABIDOPSIS (Arabidopsis thaliana) - GitHub
download_arabidopsis() {
    log "=== Downloading Arabidopsis 69 Pangenome ==="
    local DIR="$BASE_DIR/arabidopsis/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    # Download from GitHub release assets
    # wget -c "https://github.com/qclian/Pan_Ath/releases/download/.../genomes.tar.gz"
    log "Arabidopsis: Check https://github.com/qclian/Pan_Ath/releases for exact asset URLs."
}

# 8. RICE (Oryza sativa) - Figshare / ENA
download_rice() {
    log "=== Downloading Rice 149 Pangenome ==="
    local DIR="$BASE_DIR/rice/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    # Figshare API or direct wget if known
    # wget -c "https://doi.org/10.25452/figshare.plus.25697817"
    log "Rice: Check Figshare 25697817 or ENA PRJEB73710 for exact fasta URLs."
}

main() {
    log "Starting data downloads..."
    download_cucumber
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
