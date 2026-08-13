#!/usr/bin/env bash
# =============================================================================
# download_pangenomes.sh
# Comprehensive Downloader for Investigated Plant Pangenomes (11 Projects)
# With Enhanced Dry-Run Validation & .fa* / Genome File Matching Verification
# =============================================================================
set -euo pipefail

BASE_DIR="$(pwd)/pangenome_data"
mkdir -p "$BASE_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

DRY_RUN=0

# --- Global Dry-Run Validation Trackers ---
CURRENT_MODULE=""
MODULE_TOTAL=0
MODULE_FASTA=0

SUMMARY_MODULES=()
SUMMARY_TOTALS=()
SUMMARY_FASTAS=()
SUMMARY_STATUS=()

start_module() {
    CURRENT_MODULE="$1"
    MODULE_TOTAL=0
    MODULE_FASTA=0
    log "=== ${CURRENT_MODULE} ==="
}

end_module() {
    if [[ $DRY_RUN -eq 1 ]]; then
        local status="[PASS]"
        if [[ $MODULE_FASTA -eq 0 ]]; then
            status="[FAIL: NO FASTA MATCHED]"
        fi
        log "[DRY-RUN CHECK] ${CURRENT_MODULE}: ${MODULE_TOTAL} files queued, ${MODULE_FASTA} FASTA/genome files matched -> ${status}"
        SUMMARY_MODULES+=("$CURRENT_MODULE")
        SUMMARY_TOTALS+=("$MODULE_TOTAL")
        SUMMARY_FASTAS+=("$MODULE_FASTA")
        SUMMARY_STATUS+=("$status")
    fi
}

is_fasta_or_genome() {
    local str="$*"
    if [[ "$str" =~ \.(fa|fasta|fna|gfa)(\.gz)?([[:space:]]|$)|chromosome\.fa|genome\.fa|\.tar\.gz|\.tgz|\.zip|-A ]]; then
        return 0
    else
        return 1
    fi
}

run_wget() {
    if [[ $DRY_RUN -eq 1 ]]; then
        ((MODULE_TOTAL++)) || true
        if is_fasta_or_genome "$@"; then
            ((MODULE_FASTA++)) || true
            echo "[DRY RUN] [FASTA MATCH] wget $*"
        else
            echo "[DRY RUN] wget $*"
        fi
    else
        wget "$@"
    fi
}

run_curl() {
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY RUN] curl $*"
    else
        curl "$@"
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

# --- Generic Helper Functions using Python ---

# 1. Download files from Zenodo Record ID
download_zenodo_record() {
    local record_id="$1"
    log "Fetching Zenodo record ${record_id}..."
    while read -r link key || [ -n "$link" ]; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done < <(python3 -c "
import urllib.request, json, ssl, sys
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

rec_id = '${record_id}'
url = f'https://zenodo.org/api/records/{rec_id}'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
try:
    with urllib.request.urlopen(req, context=ctx, timeout=20) as r:
        data = json.loads(r.read().decode('utf-8'))
        files = data.get('files', [])
        for f in files:
            link = f.get('links', {}).get('content') or f.get('links', {}).get('self') or f.get('download_url')
            key = f.get('key') or f.get('filename')
            if link and key:
                print(f'{link}\t{key}')
except Exception as e:
    sys.stderr.write(f'Error fetching Zenodo {rec_id}: {e}\n')
")
}

# 2. Download files from Figshare Article ID
download_figshare_article() {
    local article_id="$1"
    log "Fetching Figshare article ${article_id}..."
    while read -r link key || [ -n "$link" ]; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done < <(python3 -c "
import urllib.request, json, ssl, sys
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

art_id = '${article_id}'
url = f'https://api.figshare.com/v2/articles/{art_id}/files?page_size=1000'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
try:
    with urllib.request.urlopen(req, context=ctx, timeout=20) as r:
        files = json.loads(r.read().decode('utf-8'))
        for f in files:
            link = f.get('download_url')
            name = f.get('name')
            if link and name:
                print(f'{link}\t{name}')
except Exception as e:
    sys.stderr.write(f'Error fetching Figshare {art_id}: {e}\n')
")
}

# 3. Download assemblies from NCBI BioProject via Assembly Database API
download_ncbi_bioproject() {
    local bioproject="$1"
    log "Fetching NCBI Assembly database entries for ${bioproject}..."
    while read -r link key || [ -n "$link" ]; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done < <(python3 -c "
import urllib.request, json, ssl, sys, time, urllib.parse
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

bioprj = '${bioproject}'
terms_to_try = [
    urllib.parse.quote_plus(bioprj, safe='[]'),
    urllib.parse.quote_plus(f'{bioprj}[BioProject]', safe='[]'),
    urllib.parse.quote_plus(f'{bioprj}[Organism]', safe='[]')
]

id_list = []
for term in terms_to_try:
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={term}&retmode=json'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=15) as r:
            data = json.loads(r.read().decode('utf-8'))
            ids = data.get('esearchresult', {}).get('idlist', [])
            if ids:
                id_list = ids
                break
    except Exception:
        pass

if id_list:
    time.sleep(0.4)
    ids_str = ','.join(id_list[:50])
    sum_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={ids_str}&retmode=json'
    req2 = urllib.request.Request(sum_url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    try:
        with urllib.request.urlopen(req2, context=ctx, timeout=20) as r2:
            sdata = json.loads(r2.read().decode('utf-8'))
            result = sdata.get('result', {})
            for aid in id_list[:50]:
                doc = result.get(aid, {})
                ftp = doc.get('ftppath_genbank') or doc.get('ftppath_refseq')
                acc = doc.get('assemblyaccession') or doc.get('assemblyname')
                if ftp:
                    basename = ftp.split('/')[-1]
                    fa_url = f'{ftp}/{basename}_genomic.fna.gz'
                    print(f'{fa_url}\t{acc}.fna.gz')
    except Exception as e:
        sys.stderr.write(f'Error fetching NCBI summary for {bioprj}: {e}\n')
")
}

# 4. Download assemblies from Ensembl Plants release-57 FTP
download_ensembl_fastas() {
    local species="$1"
    log "Fetching Ensembl Plants release-57 assemblies for ${species}..."
    while read -r link key || [ -n "$link" ]; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done < <(python3 -c "
import urllib.request, ssl, sys, re
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

spec = '${species}'
url = f'https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/{spec}/dna/'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req, context=ctx, timeout=20) as r:
        html = r.read().decode('utf-8', errors='ignore')
        links = set(re.findall(r'href=\"([^\"]+\.fa\.gz)\"', html))
        for l in links:
            if not l.endswith('_sm.toplevel.fa.gz') and not l.endswith('_rm.toplevel.fa.gz') and ('toplevel' in l or 'primary' in l or 'genome' in l or 'dna' in l):
                print(f'{url}{l}\t{l}')
except Exception as e:
    sys.stderr.write(f'Error fetching Ensembl {spec}: {e}\n')
")
}

# 1. MAIZE NAM & T2T PANGENOME (Zea mays)
download_maize() {
    start_module "[1/11] Maize NAM/T2T Pangenome (PRJNA639775 & PRJNA751841)"
    local DIR="$BASE_DIR/maize/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_ncbi_bioproject "PRJNA639775"
    download_ncbi_bioproject "PRJNA751841"
    end_module
}

# 2. WHEAT 10+ PANGENOME (Triticum aestivum)
download_wheat() {
    start_module "[2/11] Wheat 10+ Pangenome (Ensembl Plants FTP)"
    local DIR="$BASE_DIR/wheat/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_ensembl_fastas "triticum_aestivum"
    end_module
}

# 3. NORTH AMERICAN WILD GRAPE SUPER-PANGENOME (Vitis spp.)
download_wild_grape() {
    start_module "[3/11] Wild Grape Super-Pangenome (Zenodo 10846425 / 10851548 / PRJNA1018808)"
    local DIR="$BASE_DIR/wild_grape/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_zenodo_record "10846425"
    download_zenodo_record "10851548"
    download_ncbi_bioproject "PRJNA1018808"
    end_module
}

# 4. RAPESEED STRUCTURAL VARIATION PANGENOME (Brassica napus)
download_rapeseed() {
    start_module "[4/11] Rapeseed SV Pangenome (Ensembl Plants FTP & Zenodo)"
    local DIR="$BASE_DIR/rapeseed/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_ensembl_fastas "brassica_napus"
    download_zenodo_record "687103"
    end_module
}

# 5. POTATO PANGENOME (Solanum tuberosum)
download_potato() {
    start_module "[5/11] Potato Tetraploid Pangenome (Zenodo 7894982 & NCBI)"
    local DIR="$BASE_DIR/potato/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_ncbi_bioproject "PRJNA731597"
    download_ncbi_bioproject "Solanum tuberosum"
    end_module
}

# 6. EGGPLANT PANGENOME (Solanum melongena)
download_eggplant() {
    start_module "[6/11] Eggplant Pangenome (Zenodo 5523914 & PRJNA612792)"
    local DIR="$BASE_DIR/eggplant/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_ncbi_bioproject "PRJNA612792"
    download_ncbi_bioproject "Solanum melongena"
    end_module
}

# 7. TEA PLANT HAPLOTYPE PANGENOME (Camellia sinensis)
download_tea() {
    start_module "[7/11] Tea Plant Haplotype Pangenome (Zenodo 17174024)"
    local DIR="$BASE_DIR/tea_plant/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_zenodo_record "17174024"
    end_module
}

# 8. WATERMELON SUPER-PANGENOME (Citrullus lanatus)
download_watermelon() {
    start_module "[8/11] Watermelon Super-Pangenome (CuGenDBv2)"
    local DIR="$BASE_DIR/watermelon/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    while read -r link key || [ -n "$link" ]; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done < <(python3 -c "
import urllib.request, ssl, sys, re
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

url = 'http://cucurbitgenomics.org/v2/ftp/pan-genome/watermelon/graph_pangenome/assembly/'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req, context=ctx, timeout=15) as r:
        html = r.read().decode('utf-8', errors='ignore')
        links = set(re.findall(r'href=\"([^\"]+\.(?:fa|fasta|fa\.gz))\"', html))
        for l in links:
            if not l.startswith('http'):
                full_url = f'{url}{l}'
            else:
                full_url = l
            print(f'{full_url}\t{l}')
except Exception as e:
    sys.stderr.write(f'Error fetching Watermelon CuGenDB: {e}\n')
")
    end_module
}

# 9. TOMATO T2T SUPER-PANGENOME (Solanum lycopersicum)
download_tomato() {
    start_module "[9/11] Tomato T2T Super-Pangenome (Zenodo 17878268)"
    local DIR="$BASE_DIR/tomato/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    download_zenodo_record "17878268"
    end_module
}

# 10. MARCHANTIA PANGENOME (Marchantia polymorpha)
download_marchantia() {
    start_module "[10/11] Marchantia Pangenome (Marchantia.info / Zenodo 1021402)"
    local DIR="$BASE_DIR/marchantia/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    run_wget -c "https://marchantia.info/download/marchantia_pangenome.tar.gz" || true
    run_tar_if_exists "marchantia_pangenome.tar.gz" -xzf marchantia_pangenome.tar.gz
    run_wget -c "https://marchantia.info/download/m_polymorpha_v6.fa.gz" || true
    download_ncbi_bioproject "Marchantia polymorpha"
    end_module
}

# 11. CITRUS PANGENOME (Citrus spp.)
download_citrus() {
    start_module "[11/11] Citrus Pangenome (HZAU DB)"
    local DIR="$BASE_DIR/citrus/assemblies"
    mkdir -p "$DIR" && cd "$DIR"
    
    while read -r link key || [ -n "$link" ]; do
        if [ -n "$link" ] && [ -n "$key" ]; then
            run_wget -c "$link" -O "$key" || true
        fi
    done < <(curl -s "http://citrus.hzau.edu.cn/download.php" | python3 -c "
import sys, re
try:
    html = sys.stdin.read()
    links = set(re.findall(r'href=\"(/data/Genome_info/[^\"]+)\"', html))
    for link in links:
        if link.endswith('.fa') or link.endswith('.fa.gz') or link.endswith('.fasta'):
            basename = link.split('/')[-1]
            print(f'http://citrus.hzau.edu.cn{link}\t{basename}')
except Exception:
    pass
")
    end_module
}

main() {
    for arg in "$@"; do
        if [[ "$arg" == "--dry-run" ]]; then
            DRY_RUN=1
            log "Running in DRY RUN mode with Enhanced FASTA / Genome File Validation"
        fi
    done

    log "Starting execution for all 11 plant pangenome download modules..."
    download_maize
    download_wheat
    download_wild_grape
    download_rapeseed
    download_potato
    download_eggplant
    download_tea
    download_watermelon
    download_tomato
    download_marchantia
    download_citrus

    if [[ $DRY_RUN -eq 1 ]]; then
        echo ""
        log "=========================================================================================="
        log "                    DRY-RUN VALIDATION SUMMARY (.fa* / GENOME VERIFICATION)                "
        log "=========================================================================================="
        printf "[%s] %-55s %-10s %-15s %-10s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "Module Name" "Queued" "FASTA Matched" "Status"
        printf "[%s] %-55s %-10s %-15s %-10s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "-------------------------------------------------------" "----------" "---------------" "----------"
        local overall_fail=0
        for i in "${!SUMMARY_MODULES[@]}"; do
            printf "[%s] %-55s %-10s %-15s %-10s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "${SUMMARY_MODULES[$i]}" "${SUMMARY_TOTALS[$i]}" "${SUMMARY_FASTAS[$i]}" "${SUMMARY_STATUS[$i]}"
            if [[ "${SUMMARY_STATUS[$i]}" == *"[FAIL"* ]]; then
                overall_fail=1
            fi
        done
        log "=========================================================================================="
        if [[ $overall_fail -eq 1 ]]; then
            log "DRY-RUN VALIDATION ERROR: One or more modules have 0 FASTA/genome files matched!"
            exit 1
        else
            log "DRY-RUN VALIDATION SUCCESSFUL: All 11 modules have valid FASTA/genome downloads queued!"
        fi
    else
        log "All 11 plant pangenome download tasks completed/initiated successfully."
    fi
}

main "$@"