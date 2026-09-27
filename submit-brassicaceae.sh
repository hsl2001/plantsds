#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-brassicaceae}"
THREADS="${THREADS:-64}"
MEM="${MEM:-120gb}"
WALLTIME="${WALLTIME:-96:00:00}"

mkdir -p "$LOGDIR"
cd "$WORKDIR"
mkdir -p selected/brassicaceae/clean_fasta selected/brassicaceae/ncbi results/brassicaceae

python3 - <<'PY'
import json
import os
import re
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


assemblies = [
    ("Col-0", "Arabidopsis thaliana", "GCA_978657495.1", "TAIR12"),
    # ("Ler-0", "Arabidopsis thaliana", "GCA_946406525.1", "Ler-0.7213.PacbioHiFiAssembly"),
    # ("Arabidopsis_lyrata", "Arabidopsis lyrata", "GCA_975146905.1", "PU6_v3"),
    # ("Capsella_rubella", "Capsella rubella", "GCA_059238525.1", "Crubella_RAM_v1.0"),
    # ("Arabis_alpina", "Arabis alpina", "GCA_900128785.1", "MPIPZ.v5"),
    # ("Brassica_rapa", "Brassica rapa", "GCA_060414235.1", "ASM6041423v1"),
    ("Aethionema_arabicum", "Aethionema arabicum", "GCA_000411095.1", "VEGI_AA_v_1.0"),
]
base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
root = Path("selected/brassicaceae")
raw_dir = root / "ncbi"
clean_dir = root / "clean_fasta"
organelle = re.compile(r"^(?:chr)?(?:c|m|mt|mit|mito|mitochondria|mitochondrion|pt|cp|pltd|plastid)$", re.I)
organelle_words = ("chloroplast", "mitochondrion", "mitochondrial", "plastid", "plastome", "chondriome")


def request(url):
    for attempt in range(4):
        try:
            return urllib.request.urlopen(url, timeout=180)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as error:
            if attempt == 3 or isinstance(error, urllib.error.HTTPError) and error.code not in (429, 500, 502, 503, 504):
                raise
            time.sleep(2 ** attempt)


def report_for(accession, species, assembly_name):
    with request(base_url + accession + "/dataset_report") as response:
        reports = json.load(response).get("reports", [])
    match = next((entry for entry in reports if entry.get("accession") == accession), None)
    if not match:
        raise RuntimeError(f"NCBI assembly not found: {accession}")
    name = match.get("organism", {}).get("organism_name", "")
    info = match.get("assembly_info", {})
    if not name.startswith(species) or info.get("assembly_name") != assembly_name:
        raise RuntimeError(f"NCBI assembly identity changed: {accession}: {name}, {info.get('assembly_name')}")
    return match


def genomic_fasta(archive, accession):
    with zipfile.ZipFile(archive) as zipped:
        members = [member for member in zipped.namelist()
                   if member.startswith(f"ncbi_dataset/data/{accession}/")
                   and (member.endswith("_genomic.fna") or member.endswith("_genomic.fna.gz"))]
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
    url = base_url + accession + "/download?include_annotation_type=GENOME_FASTA"
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
for label, species, accession, assembly_name in assemblies:
    report = report_for(accession, species, assembly_name)
    raw = raw_dir / f"{accession}_genomic.fna"
    clean = clean_dir / f"{accession}_{label}.nuclear.fna"
    if not raw.is_file() or raw.stat().st_size == 0:
        download(accession, raw)
    if not clean.is_file() or clean.stat().st_size == 0 or clean.stat().st_mtime < raw.stat().st_mtime:
        clean_fasta(raw, clean, accession)
    records.append((label, species, accession, assembly_name,
                    report.get("assembly_info", {}).get("assembly_level", ""),
                    str(raw.resolve()), str(clean.resolve())))
    print(f"{label}: {accession} ({assembly_name}) -> {clean}", file=sys.stderr)

with (root / "genomes.tsv").open("w") as handle:
    handle.write("label\tspecies\taccession\tassembly\tlevel\traw_fasta\tused_fasta\n")
    for record in records:
        handle.write("\t".join(record) + "\n")
with (root / "genomes.files").open("w") as handle:
    for record in records:
        handle.write(record[-1] + "\n")
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
mapfile -t FASTAS < selected/brassicaceae/genomes.files

output="results/brassicaceae/BRASSICACEAE_07_Aethionema_arabicum"
printf '[segtrace] one run with %d genomes: %s\n' "${#FASTAS[@]}" "$output"
./segtrace -p "$THREADS" -o "$output" "${FASTAS[@]}"
PBS