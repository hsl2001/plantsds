#!/usr/bin/env python3
import argparse
import glob
import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ANGIOSPERM_TAXID = 3398
VALID_CHR_LEVELS = {"complete genome", "chromosome"}

FORCED_MAJOR_PLANTS = {
    "Am": ("GCF_000471905.2", "Amborella trichopoda"),
    "At": ("GCF_000001735.4", "Arabidopsis thaliana"),
    "Bd": ("GCF_000005505.3", "Brachypodium distachyon"),
    "Br": ("GCF_000309985.2", "Brassica rapa"),
    "Os": ("GCF_034140825.1", "Oryza sativa"),
    "Pt": ("GCF_000002775.5", "Populus trichocarpa"),
    "Sl": ("GCF_036512215.1", "Solanum lycopersicum"),
    "Vv": ("GCF_030704535.1", "Vitis vinifera"),
}

ORGANELLE_KEYWORDS = [
    "chloroplast",
    "mitochondrion",
    "mitochondrial",
    "chloroplastic",
    "plastid",
    "plastome",
    "chondriome",
    "[location=chloroplast]",
    "[location=mitochondrion]",
    "[location=plastid]",
]

ORGANELLE_ID_PATTERN = re.compile(
    r"^(chr)?(m|mt|mit|mito|mitochondria|mitochondrion|pt|cp|pltd|plastid)$",
    re.IGNORECASE,
)


def is_organelle_header(header_line: str) -> bool:
    line_clean = header_line.lstrip(">").strip()
    if not line_clean:
        return False
    seq_id = line_clean.split()[0]
    if ORGANELLE_ID_PATTERN.match(seq_id):
        return True
    header_lower = line_clean.lower()
    return any(kw in header_lower for kw in ORGANELLE_KEYWORDS)


def filter_organelle_from_fasta(in_path: Path, out_path: Path) -> tuple[int, int, list[str]]:
    """Filters out organelle sequences from a FASTA file.
    Returns (kept_count, discarded_count, discarded_names).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    open_in = gzip.open if str(in_path).endswith(".gz") else open

    kept = 0
    discarded = 0
    discarded_names = []

    is_writing = False
    with open_in(in_path, "rt", encoding="utf-8", errors="replace") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                if is_organelle_header(line):
                    is_writing = False
                    discarded += 1
                    seq_id = line.strip().split()[0].lstrip(">")
                    discarded_names.append(seq_id)
                else:
                    is_writing = True
                    kept += 1
                    fout.write(line)
            elif is_writing:
                fout.write(line)

    return kept, discarded, discarded_names


def completeness_key(record):
    level = record.get("level", "").lower()
    level_score = 2 if level == "complete genome" else (1 if level == "chromosome" else 0)
    is_ref = 1 if record.get("category", "").lower() in ("reference genome", "representative genome") else 0
    busco = record.get("busco_complete", 0.0)
    scaffold_n50 = record.get("scaffold_n50", 0)
    contig_n50 = record.get("contig_n50", 0)
    size = record.get("size", 0)
    return (level_score, is_ref, busco, scaffold_n50, contig_n50, size)


def read_records(report_path):
    records = []
    taxids = set()
    with open(report_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            report = json.loads(line)
            tax_id = report.get("organism", {}).get("taxId")
            stats = report.get("assemblyStats", {})
            size = int(stats.get("totalSequenceLength") or stats.get("totalUngappedLength") or 0)
            accession = report.get("accession") or report.get("currentAccession")
            if not accession or not tax_id or size <= 0:
                continue
            level = report.get("assemblyInfo", {}).get("assemblyLevel", "")
            category = report.get("assemblyInfo", {}).get("assemblyCategory", "")
            scaffold_n50 = int(stats.get("scaffoldN50") or 0)
            contig_n50 = int(stats.get("contigN50") or 0)
            busco_obj = stats.get("busco", {})
            busco_complete = float(busco_obj.get("complete") or busco_obj.get("buscoScore") or 0.0)

            records.append({
                "accession": accession,
                "tax_id": int(tax_id),
                "name": report.get("organism", {}).get("organismName", ""),
                "size": size,
                "level": level,
                "category": category,
                "scaffold_n50": scaffold_n50,
                "contig_n50": contig_n50,
                "busco_complete": busco_complete,
            })
            taxids.add(int(tax_id))
    return records, sorted(taxids)


def fetch_taxonomy(taxids, cache_path):
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as out:
        out.write("tax_id\tscientific_name\ttaxid_rank_name\n")
        for start in range(0, len(taxids), 250):
            batch = taxids[start:start + 250]
            params = urllib.parse.urlencode({
                "db": "taxonomy",
                "id": ",".join(map(str, batch)),
                "retmode": "xml",
            })
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
            root = None
            for attempt in range(5):
                try:
                    with urllib.request.urlopen(url, timeout=60) as response:
                        root = ET.fromstring(response.read())
                    break
                except Exception as exc:
                    if attempt == 4:
                        raise RuntimeError(f"failed to fetch taxonomy batch {start}: {exc}")
                    time.sleep(1.2 * (attempt + 1))
            for taxon in root.findall("./Taxon"):
                tax_id = taxon.findtext("TaxId") or ""
                sci = (taxon.findtext("ScientificName") or "").replace("\t", " ")
                nodes = []
                lineage = taxon.find("LineageEx")
                if lineage is not None:
                    for node in lineage.findall("./Taxon"):
                        node_id = node.findtext("TaxId") or ""
                        rank = (node.findtext("Rank") or "").replace("\t", " ")
                        name = (node.findtext("ScientificName") or "").replace("\t", " ")
                        nodes.append(f"{node_id}|{rank}|{name}")
                rank = (taxon.findtext("Rank") or "").replace("\t", " ")
                nodes.append(f"{tax_id}|{rank}|{sci}")
                out.write(f"{tax_id}\t{sci}\t{';'.join(nodes)}\n")
            print(f"[select] taxonomy {min(start + 250, len(taxids))}/{len(taxids)}", file=sys.stderr)
            time.sleep(0.34)


def read_taxonomy(cache_path):
    ranks_by_tax = {}
    lineage_by_tax = {}
    with open(cache_path) as handle:
        next(handle)
        for line in handle:
            tax_id_s, _, nodes_s = line.rstrip("\n").split("\t")
            tax_id = int(tax_id_s)
            ranks = {}
            lineage = set()
            for node in nodes_s.split(";"):
                parts = node.split("|", 2)
                if len(parts) != 3:
                    continue
                node_id_s, rank, name = parts
                if node_id_s:
                    node_id = int(node_id_s)
                    lineage.add(node_id)
                    ranks.setdefault(rank, (node_id, name))
            ranks_by_tax[tax_id] = ranks
            lineage_by_tax[tax_id] = lineage
    return ranks_by_tax, lineage_by_tax


def assembly_path(dataset_dir, accession):
    matches = sorted(glob.glob(str(dataset_dir / accession / "*_genomic.fna*")))
    return matches[0] if matches else ""


def available_forced_record(preferred_accession, organism_name, records,
                            records_by_accession, dataset_dir):
    preferred = records_by_accession.get(preferred_accession)
    if preferred and assembly_path(dataset_dir, preferred_accession):
        return preferred

    candidates = [
        record for record in records
        if record["name"].lower().startswith(organism_name.lower()) and
        assembly_path(dataset_dir, record["accession"])
    ]
    if not candidates:
        genus = organism_name.split()[0].lower()
        candidates = [
            record for record in records
            if record["name"].lower().startswith(genus + " ") and
            assembly_path(dataset_dir, record["accession"])
        ]
    if candidates:
        chr_candidates = [c for c in candidates if c.get("level", "").lower() in VALID_CHR_LEVELS]
        pool = chr_candidates if chr_candidates else candidates
        fallback = max(pool, key=completeness_key)
        print(f"[select] {preferred_accession} unavailable; using "
              f"{fallback['accession']} for {organism_name}", file=sys.stderr)
        return fallback
    return preferred


def main():
    parser = argparse.ArgumentParser(
        description="Select the most complete chromosome-level assembly per angiosperm family and discard organelle genomes."
    )
    parser.add_argument("--dataset-dir", default="ncbi_dataset/data")
    parser.add_argument("--report", default=None)
    parser.add_argument("--taxonomy-cache", default="selected/taxonomy_rank_lineage.tsv")
    parser.add_argument("--out", default="selected/angiosperm_family.files")
    parser.add_argument("--summary", default="selected/angiosperm_family_complete.tsv")
    parser.add_argument("--clean-dir", default="selected/clean_fasta",
                        help="Directory to save organelle-filtered FASTA files.")
    parser.add_argument("--no-discard-organelle", action="store_true",
                        help="Do not filter out organelle genomes.")
    parser.add_argument("--force-accession", action="append", default=[],
                        help="Extra assembly accession to force include. May be repeated.")
    parser.add_argument("--allow-missing", action="store_true",
                        help="Do not fail when selected NCBI FASTA files are absent.")
    parser.add_argument("--allow-lower-levels", action="store_true",
                        help="Allow scaffold/contig level assemblies if no chromosome level is available.")
    args = parser.parse_args()

    dataset_dir = Path(args.dataset_dir)
    report_path = Path(args.report) if args.report else dataset_dir / "assembly_data_report.jsonl"
    taxonomy_cache = Path(args.taxonomy_cache)
    out_path = Path(args.out)
    summary_path = Path(args.summary)
    clean_dir = Path(args.clean_dir)
    discard_organelles = not args.no_discard_organelle

    records, taxids = read_records(report_path)
    records_by_accession = {record["accession"]: record for record in records}
    if not taxonomy_cache.exists():
        fetch_taxonomy(taxids, taxonomy_cache)
    ranks_by_tax, lineage_by_tax = read_taxonomy(taxonomy_cache)

    best = {}
    for record in records:
        lineage = lineage_by_tax.get(record["tax_id"], set())
        if ANGIOSPERM_TAXID not in lineage:
            continue
        family = ranks_by_tax.get(record["tax_id"], {}).get("family")
        if not family:
            continue
        
        # Chr level 이상 (Chromosome, Complete Genome) 만 사용
        if not args.allow_lower_levels and record.get("level", "").lower() not in VALID_CHR_LEVELS:
            continue

        if family not in best or completeness_key(record) > completeness_key(best[family]):
            best[family] = record

    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    if discard_organelles:
        clean_dir.mkdir(parents=True, exist_ok=True)

    selected = {record["accession"]: ("family_complete", family, record)
                for family, record in best.items()}
    for label, (accession, organism_name) in FORCED_MAJOR_PLANTS.items():
        record = available_forced_record(accession, organism_name, records,
                                         records_by_accession, dataset_dir)
        if record is None:
            raise SystemExit(f"forced accession for {label} is absent from report: {accession}")
        selected[record["accession"]] = (f"forced_{label}", None, record)
    for accession in args.force_accession:
        record = records_by_accession.get(accession)
        if record is None:
            raise SystemExit(f"forced accession is absent from report: {accession}")
        selected[accession] = ("forced_user", None, record)

    missing = []
    fasta_paths = []
    total_organelles_discarded = 0

    with open(summary_path, "w") as summary:
        summary.write(
            "reason\tfamily_tax_id\tfamily\taccession\tsize_bp\tassembly_level\t"
            "scaffold_n50\tcontig_n50\tbusco_complete\torganism\traw_fasta\tused_fasta\tdiscarded_organelle_count\n"
        )
        for accession, (reason, family, record) in sorted(selected.items(), key=lambda item: item[1][2]["name"]):
            raw_fasta = assembly_path(dataset_dir, record["accession"])
            used_fasta = ""
            discarded_count = 0
            if raw_fasta:
                if discard_organelles:
                    stem = Path(raw_fasta).name
                    if stem.endswith(".gz"):
                        stem = stem[:-3]
                    if stem.endswith(".fna") or stem.endswith(".fasta") or stem.endswith(".fa"):
                        stem = Path(stem).stem
                    clean_fasta = clean_dir / f"{record['accession']}_{stem}.nuclear.fna"
                    kept, discarded_count, discarded_names = filter_organelle_from_fasta(Path(raw_fasta), clean_fasta)
                    used_fasta = str(clean_fasta)
                    total_organelles_discarded += discarded_count
                    if discarded_count > 0:
                        print(f"[organelle] {record['accession']}: discarded {discarded_count} organelle sequences: {','.join(discarded_names)}", file=sys.stderr)
                else:
                    used_fasta = raw_fasta
                fasta_paths.append(used_fasta)
            else:
                missing.append(record["accession"])

            family_tax_id = family[0] if family else ""
            family_name = family[1] if family else ""
            summary.write(
                f"{reason}\t{family_tax_id}\t{family_name}\t{record['accession']}\t{record['size']}\t"
                f"{record['level']}\t{record.get('scaffold_n50', 0)}\t{record.get('contig_n50', 0)}\t"
                f"{record.get('busco_complete', 0.0)}\t{record['name']}\t{raw_fasta}\t{used_fasta}\t{discarded_count}\n"
            )

    with open(out_path, "w") as out:
        for path in fasta_paths:
            out.write(path + "\n")

    total_bp = sum(record["size"] for record in best.values())
    selected_bp = sum(record["size"] for _, _, record in selected.values())
    print(f"selected_families={len(best)}")
    print(f"selected_angiosperm_bp={total_bp}")
    print(f"selected_angiosperm_Gbp={total_bp / 1e9:.3f}")
    print(f"forced_major_plants={len(FORCED_MAJOR_PLANTS)}")
    print(f"selected_total_bp={selected_bp}")
    print(f"selected_total_Gbp={selected_bp / 1e9:.3f}")
    print(f"written_fastas={len(fasta_paths)}")
    print(f"discarded_organelles_total={total_organelles_discarded}")
    print(f"missing_ncbi_fastas={len(missing)}")
    if missing and not args.allow_missing:
        print("missing accessions: " + ",".join(missing[:20]), file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
