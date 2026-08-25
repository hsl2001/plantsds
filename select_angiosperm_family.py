#!/usr/bin/env python3
import argparse
import glob
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ANGIOSPERM_TAXID = 3398

FORCED_MAJOR_PLANT_ACCESSIONS = {
    "Am": "GCF_000471905.2",  # Amborella trichopoda
    "At": "GCF_000001735.4",  # Arabidopsis thaliana
    "Bd": "GCF_000005505.3",  # Brachypodium distachyon
    "Br": "GCF_000309985.2",  # Brassica rapa
    # "Cr": "GCF_000002595.2",  # Chlamydomonas reinhardtii
    # "Cv": "GCA_023343905.1",  # Chlorella vulgaris
    # "Mp": "GCA_037833965.1",  # Marchantia polymorpha
    "Os": "GCF_034140825.1",  # Oryza sativa Japonica Group
    # "Pp": "GCF_000002425.5",  # Physcomitrium patens
    "Pt": "GCF_000002775.5",  # Populus trichocarpa
    "Sl": "GCF_036512215.1",  # Solanum lycopersicum
    # "Sm": "GCF_000143415.4",  # Selaginella moellendorffii
    # "Vc": "GCF_000143455.1",  # Volvox carteri
    "Vv": "GCF_030704535.1",  # Vitis vinifera
}


def read_records(report_path):
    records = []
    taxids = set()
    with open(report_path) as handle:
        for line in handle:
            report = json.loads(line)
            tax_id = report.get("organism", {}).get("taxId")
            stats = report.get("assemblyStats", {})
            size = int(stats.get("totalSequenceLength") or stats.get("totalUngappedLength") or 0)
            accession = report.get("accession") or report.get("currentAccession")
            if not accession or not tax_id or size <= 0:
                continue
            records.append({
                "accession": accession,
                "tax_id": int(tax_id),
                "name": report.get("organism", {}).get("organismName", ""),
                "size": size,
                "level": report.get("assemblyInfo", {}).get("assemblyLevel", ""),
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


def main():
    parser = argparse.ArgumentParser(description="Select the smallest assembly per angiosperm family.")
    parser.add_argument("--dataset-dir", default="ncbi_dataset/data")
    parser.add_argument("--report", default=None)
    parser.add_argument("--taxonomy-cache", default="selected/taxonomy_rank_lineage.tsv")
    parser.add_argument("--out", default="selected/angiosperm_family_min.files")
    parser.add_argument("--summary", default="selected/angiosperm_family_min.tsv")
    parser.add_argument("--force-accession", action="append", default=[],
                        help="Extra assembly accession to force include. May be repeated.")
    parser.add_argument("--allow-missing", action="store_true",
                        help="Do not fail when selected NCBI FASTA files are absent.")
    args = parser.parse_args()

    dataset_dir = Path(args.dataset_dir)
    report_path = Path(args.report) if args.report else dataset_dir / "assembly_data_report.jsonl"
    taxonomy_cache = Path(args.taxonomy_cache)
    out_path = Path(args.out)
    summary_path = Path(args.summary)

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
        if family not in best or record["size"] < best[family]["size"]:
            best[family] = record

    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    selected = {record["accession"]: ("family_min", family, record)
                for family, record in best.items()}
    for label, accession in FORCED_MAJOR_PLANT_ACCESSIONS.items():
        record = records_by_accession.get(accession)
        if record is None:
            raise SystemExit(f"forced accession for {label} is absent from report: {accession}")
        selected[accession] = (f"forced_{label}", None, record)
    for accession in args.force_accession:
        record = records_by_accession.get(accession)
        if record is None:
            raise SystemExit(f"forced accession is absent from report: {accession}")
        selected[accession] = ("forced_user", None, record)

    missing = []
    fasta_paths = []
    with open(summary_path, "w") as summary:
        summary.write("reason\tfamily_tax_id\tfamily\taccession\tsize_bp\tassembly_level\torganism\tfasta\n")
        for accession, (reason, family, record) in sorted(selected.items(), key=lambda item: item[1][2]["name"]):
            fasta = assembly_path(dataset_dir, record["accession"])
            if fasta:
                fasta_paths.append(fasta)
            else:
                missing.append(record["accession"])
            family_tax_id = family[0] if family else ""
            family_name = family[1] if family else ""
            summary.write(
                f"{reason}\t{family_tax_id}\t{family_name}\t{record['accession']}\t{record['size']}\t"
                f"{record['level']}\t{record['name']}\t{fasta}\n")

    with open(out_path, "w") as out:
        for path in fasta_paths:
            out.write(path + "\n")

    total_bp = sum(record["size"] for record in best.values())
    selected_bp = sum(record["size"] for _, _, record in selected.values())
    print(f"selected_families={len(best)}")
    print(f"selected_angiosperm_bp={total_bp}")
    print(f"selected_angiosperm_Gbp={total_bp / 1e9:.3f}")
    print(f"forced_major_plants={len(FORCED_MAJOR_PLANT_ACCESSIONS)}")
    print(f"selected_total_bp={selected_bp}")
    print(f"selected_total_Gbp={selected_bp / 1e9:.3f}")
    print(f"written_fastas={len(fasta_paths)}")
    print(f"missing_ncbi_fastas={len(missing)}")
    if missing and not args.allow_missing:
        print("missing accessions: " + ",".join(missing[:20]), file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
