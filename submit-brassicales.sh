#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-brassicales}"
THREADS="${THREADS:-64}"
MEM="${MEM:-120gb}"
WALLTIME="${WALLTIME:-96:00:00}"

mkdir -p "$LOGDIR"
cd "$WORKDIR"
mkdir -p selected/brassicales/clean_fasta selected/brassicales/ncbi results/brassicales

python3 - <<'PY'
import json
import os
import re
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path


taxa = [
    ("Brassicaceae", "Arabidopsis_thaliana", "Arabidopsis thaliana", ()),
    ("Brassicaceae", "Arabidopsis_lyrata", "Arabidopsis lyrata", ()),
    ("Brassicaceae", "Capsella_rubella", "Capsella rubella", ()),
    ("Brassicaceae", "Arabis_alpina", "Arabis alpina", ()),
    ("Brassicaceae", "Brassica_rapa", "Brassica rapa", ()),
    ("Brassicaceae", "Aethionema_arabicum", "Aethionema arabicum", ()),
    ("Cleomaceae", "Gynandropsis_gynandra", "Gynandropsis gynandra", ("Cleome gynandra",)),
    ("Caricaceae", "Carica_papaya", "Carica papaya", ()),
    ("Moringaceae", "Moringa_oleifera", "Moringa oleifera", ()),
    ("Bixaceae", "Bixa_orellana", "Bixa orellana", ()),
]
taxon_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/"
accession_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
root = Path("selected/brassicales")
raw_dir = root / "ncbi"
clean_dir = root / "clean_fasta"
level_score = {"complete genome": 3, "chromosome": 2, "scaffold": 1, "contig": 0}
organelle = re.compile(r"^(?:chr)?(?:c|m|mt|mit|mito|mitochondria|mitochondrion|pt|cp|pltd|plastid)$", re.I)
organelle_words = ("chloroplast", "mitochondrion", "mitochondrial", "plastid", "plastome", "chondriome")
last_request = 0.0


def request(url):
    global last_request
    for attempt in range(5):
        pause = 0.34 - (time.monotonic() - last_request)
        if pause > 0:
            time.sleep(pause)
        last_request = time.monotonic()
        try:
            return urllib.request.urlopen(url, timeout=180)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as error:
            if attempt == 4 or isinstance(error, urllib.error.HTTPError) and error.code not in (429, 500, 502, 503, 504):
                raise
            time.sleep(2 ** attempt)


def species_key(name):
    words = re.findall(r"[A-Za-z]+", name)
    return " ".join(words[:2]).casefold()


def reports_for(species):
    base_url = taxon_url + urllib.parse.quote(species) + "/dataset_report?filters.assembly_version=current&page_size=100"
    reports = {}
    token = ""
    while True:
        url = base_url
        if token:
            url += "&page_token=" + urllib.parse.quote(token, safe="")
        with request(url) as response:
            page = json.load(response)
        for report in page.get("reports", []):
            accession = report.get("accession", "")
            if accession.startswith(("GCA_", "GCF_")):
                reports[accession] = report
        token = page.get("next_page_token")
        if not token:
            break
    return list(reports.values())


def metric(stats, name):
    try:
        return int(stats.get(name) or 0)
    except (TypeError, ValueError):
        return 0


def busco_complete(stats):
    busco = stats.get("busco") or {}
    if not isinstance(busco, dict):
        busco = {"complete": busco}
    try:
        return float(busco.get("complete") or busco.get("buscoScore") or 0)
    except (TypeError, ValueError):
        return 0.0


def quality_key(report):
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    category = str(info.get("refseq_category", "")).lower()
    return (
        busco_complete(stats),
        metric(stats, "contig_n50"),
        metric(stats, "scaffold_n50"),
        level_score.get(str(info.get("assembly_level", "")).lower(), -1),
        int(category in ("reference genome", "representative genome")),
        metric(stats, "total_sequence_length"),
        int(report.get("accession", "").startswith("GCF_")),
        report.get("accession", ""),
    )


def select_assembly(family, label, species, aliases):
    expected = {species_key(species), *(species_key(alias) for alias in aliases)}
    candidates = [
        report for report in reports_for(species)
        if species_key(report.get("organism", {}).get("organism_name", "")) in expected
    ]
    if not candidates:
        raise RuntimeError(f"no current NCBI genome assembly found for {species}")
    report = max(candidates, key=quality_key)
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    return {
        "family": family,
        "label": label,
        "species": species,
        "accession": report["accession"],
        "assembly": info.get("assembly_name", ""),
        "level": info.get("assembly_level", ""),
        "refseq_category": info.get("refseq_category", ""),
        "busco_complete": busco_complete(stats),
        "contig_n50": metric(stats, "contig_n50"),
        "scaffold_n50": metric(stats, "scaffold_n50"),
        "size_bp": metric(stats, "total_sequence_length"),
    }


def genomic_fasta(archive, accession):
    with zipfile.ZipFile(archive) as zipped:
        members = [
            member for member in zipped.namelist()
            if member.startswith(f"ncbi_dataset/data/{accession}/")
            and (member.endswith("_genomic.fna") or member.endswith("_genomic.fna.gz"))
        ]
        if len(members) != 1:
            raise RuntimeError(f"expected one genomic FASTA for {accession}, found {members}")
        with zipped.open(members[0]) as source:
            if members[0].endswith(".gz"):
                import gzip
                with gzip.GzipFile(fileobj=source) as decoded:
                    yield from decoded
            else:
                yield from source


def download(accession, destination):
    url = accession_url + accession + "/download?include_annotation_type=GENOME_FASTA"
    with tempfile.TemporaryDirectory(dir=raw_dir) as temporary:
        archive = Path(temporary) / "genome.zip"
        extracted = Path(temporary) / "genome.fna"
        for attempt in range(4):
            try:
                with request(url) as response, archive.open("wb") as output:
                    shutil.copyfileobj(response, output)
                with extracted.open("wb") as output:
                    for chunk in genomic_fasta(archive, accession):
                        output.write(chunk)
                if extracted.stat().st_size == 0:
                    raise RuntimeError(f"empty genomic FASTA for {accession}")
                os.replace(extracted, destination)
                return
            except (OSError, zipfile.BadZipFile, RuntimeError):
                extracted.unlink(missing_ok=True)
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)


def keep_nuclear(header):
    sequence_id = header[1:].split(None, 1)[0]
    description = header.lower()
    return not organelle.fullmatch(sequence_id) and not any(word in description for word in organelle_words)


def clean_fasta(source, destination, accession):
    with tempfile.NamedTemporaryFile(mode="w", dir=clean_dir, suffix=".fna", delete=False) as temporary:
        temporary_path = Path(temporary.name)
        kept = 0
        kept_bases = 0
        write = False
        try:
            with source.open(encoding="utf-8") as original:
                for line in original:
                    if line.startswith(">"):
                        write = keep_nuclear(line)
                        if write:
                            kept += 1
                            temporary.write(f">{accession}_{line[1:]}")
                    elif write:
                        kept_bases += len(line.strip())
                        temporary.write(line)
            if kept == 0 or kept_bases == 0:
                raise RuntimeError(f"no nuclear FASTA records for {accession}: {source}")
            os.replace(temporary_path, destination)
        finally:
            temporary_path.unlink(missing_ok=True)
    return kept, kept_bases


records = []
for family, label, species, aliases in taxa:
    record = select_assembly(family, label, species, aliases)
    raw = raw_dir / f"{record['accession']}_genomic.fna"
    clean = clean_dir / f"{record['accession']}_{label}.nuclear.fna"
    if not raw.is_file() or raw.stat().st_size == 0:
        download(record["accession"], raw)
    if not clean.is_file() or clean.stat().st_size == 0 or clean.stat().st_mtime < raw.stat().st_mtime:
        clean_fasta(raw, clean, record["accession"])
    record["raw_fasta"] = str(raw.resolve())
    record["used_fasta"] = str(clean.resolve())
    records.append(record)
    print(
        f"{family}\t{species}\t{record['accession']}\t{record['assembly']}\t{record['level']}\t"
        f"contig_N50={record['contig_n50']}\tscaffold_N50={record['scaffold_n50']} -> {clean}",
        file=sys.stderr,
    )

with (root / "genomes.tsv").open("w") as handle:
    columns = (
        "family", "species", "accession", "assembly", "level", "refseq_category",
        "busco_complete", "contig_n50", "scaffold_n50", "size_bp", "raw_fasta", "used_fasta",
    )
    handle.write("\t".join(columns) + "\n")
    for record in records:
        handle.write("\t".join(str(record[column]) for column in columns) + "\n")
with (root / "genomes.files").open("w") as handle:
    for record in records:
        handle.write(record["used_fasta"] + "\n")
print(f"selected_genomes={len(records)}")
PY

qsub -N "$JOB_NAME" \
  -l "select=1:ncpus=${THREADS}:mem=${MEM}" \
  -l "walltime=${WALLTIME}" \
  -v "WORKDIR=${WORKDIR},THREADS=${THREADS}" \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
mapfile -t FASTAS < selected/brassicales/genomes.files
if [[ ${#FASTAS[@]} -ne 10 ]]; then
    echo "Expected ten Brassicales FASTAs" >&2
    exit 1
fi

printf -v genome_count '%02d' "${#FASTAS[@]}"
output="results/brassicales/BRASSICALES_${genome_count}_Bixaceae_Bixa_orellana"
printf '[segtrace] one run with %d genomes: %s\n' "${#FASTAS[@]}" "$output"
./segtrace -p "$THREADS" -c 1 -o "$output" "${FASTAS[@]}"
PBS