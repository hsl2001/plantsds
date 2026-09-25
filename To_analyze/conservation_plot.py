#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.colors import ListedColormap


ACCESSION_RE = re.compile(r"GC[AF]_\d+\.\d+")
GENOME_RE = re.compile(
    r"^(.+?)-(?:chr|chromosome|contig|scaffold|unplaced|seq\d+|tig\d+|ptg\d+|[A-Z]{2,4}_\d+\.\d+)",
    re.IGNORECASE,
)
NON_CHROMOSOMAL_RE = re.compile(r"(?:contig|scaffold|unplaced)", re.IGNORECASE)
CHROMOSOME_NUMBER_RE = re.compile(r"(?:chr|chromosome)[_-]?(\d+|[xy])", re.IGNORECASE)
GENOME_ALIASES = {"tair12": "Col-0"}
TIP_ORDER = ("Aethionema arabicum", "Arabis alpina", "Brassica rapa", "Capsella rubella", "Arabidopsis lyrata", "Ler-0", "Col-0")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot shared reference segments as a genome-by-sequence heatmap."
    )
    parser.add_argument("bed", type=Path, help="SegTrace BED with chrom, start, end, cluster_id columns")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome label, accession, or species name")
    parser.add_argument("--metadata", type=Path, help="TSV containing species and accession columns for species-name lookup")
    parser.add_argument("-o", "--output", type=Path, help="Output SVG path")
    parser.add_argument("--all-reference-sequences", action="store_true", help="Include contigs and scaffolds as panels")
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def resolve_reference(reference, metadata):
    if metadata is None:
        return reference
    with metadata.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"species", "accession"}.issubset(reader.fieldnames):
            raise ValueError("metadata must contain species and accession columns")
        matches = {row["accession"] for row in reader if row["species"].casefold() == reference.casefold()}
    if not matches:
        return reference
    if len(matches) > 1:
        raise ValueError(f"metadata maps {reference!r} to multiple accessions: {', '.join(sorted(matches))}")
    return matches.pop()


def load_genome_metadata(metadata):
    if metadata is None:
        return {}
    with metadata.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"species", "accession"}.issubset(reader.fieldnames):
            raise ValueError("metadata must contain species and accession columns")
        names = {}
        for row in reader:
            strain = next((row.get(column, "").strip() for column in ("strain", "ecotype", "cultivar", "isolate") if row.get(column, "").strip()), "")
            label = row["species"].strip()
            if strain:
                label = f"{label} ({strain})"
            names[row["accession"].strip()] = label
    return names


def display_genome_name(sequence_name, metadata_names):
    identifier = sequence_genome_id(sequence_name)
    if identifier in metadata_names:
        return metadata_names[identifier]

    aliases = {
        "t2t_nip": "Oryza sativa Nipponbare",
        "nip_hifi": "Oryza sativa Nipponbare (HiFi)",
    }
    if identifier.casefold() in aliases:
        return aliases[identifier.casefold()]

    species = re.search(r"([A-Z][a-z]+_[a-z]+)(?=[._-])", sequence_name)
    if species:
        return species.group(1).replace("_", " ")

    strain = re.search(r"GC[AF]_\d+\.\d+_([A-Za-z0-9]+(?:-[A-Za-z0-9]+)*)", sequence_name)
    if strain and strain.group(1).casefold() not in {"gca", "gcf"}:
        return strain.group(1)

    label = re.sub(r"\.transfer\.merge\.chr$", "", identifier, flags=re.IGNORECASE)
    label = re.sub(r"\.transfer\.merge$", "", label, flags=re.IGNORECASE)
    if ACCESSION_RE.fullmatch(label):
        return "Unmapped genome"
    return label.replace("_", " ")


def sequence_genome_id(sequence_name):
    accession = ACCESSION_RE.search(sequence_name)
    if accession:
        return accession.group(0)
    match = GENOME_RE.match(sequence_name)
    return match.group(1) if match else sequence_name


def genome_id(sequence_name):
    return GENOME_ALIASES.get(sequence_genome_id(sequence_name), sequence_genome_id(sequence_name))


def is_reference_sequence(sequence_name, reference):
    return sequence_genome_id(sequence_name) == reference or sequence_name == reference


def excluded_genomes(reference):
    return {"Col-0"} if reference.casefold() == "tair12" else set()


def chromosome_key(sequence_name):
    match = CHROMOSOME_NUMBER_RE.search(sequence_name)
    if not match:
        return (1, sequence_name)
    value = match.group(1).casefold()
    return (0, int(value) if value.isdigit() else {"x": 23, "y": 24}[value])


def external_sort(source, destination, unique=False):
    command = ["sort", "-t", "\t", "-k1,1", "-k2,2"]
    if unique:
        command.append("-u")
    command.extend([str(source), "-o", str(destination)])
    subprocess.run(command, check=True)


def collect_records(bed, reference, workdir, metadata_names):
    pairs = workdir / "cluster_genome.tsv"
    reference_segments = workdir / "reference_segments.tsv"
    excluded = excluded_genomes(reference)
    genome_names = {}
    with bed.open() as source, pairs.open("w") as pair_output, reference_segments.open("w") as reference_output:
        header = source.readline().rstrip("\n").split("\t")
        expected = ["#chrom", "start", "end", "cluster_id"]
        if header != expected:
            raise ValueError(f"expected BED header {' '.join(expected)!r}, found {' '.join(header)!r}")
        for line_number, line in enumerate(source, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"{bed}:{line_number} has {len(fields)} columns; expected 4")
            sequence_name, start, end, cluster_id = fields
            if sequence_genome_id(sequence_name) in excluded:
                continue
            genome = genome_id(sequence_name)
            pair_output.write(f"{cluster_id}\t{genome}\n")
            genome_names.setdefault(genome, display_genome_name(sequence_name, metadata_names))
            if is_reference_sequence(sequence_name, reference):
                reference_output.write(f"{cluster_id}\t{sequence_name}\t{start}\t{end}\n")
    label_counts = Counter(genome_names.values())
    label_seen = defaultdict(int)
    for genome, label in list(genome_names.items()):
        if label_counts[label] > 1:
            label_seen[label] += 1
            genome_names[genome] = f"{label} (sample {label_seen[label]})"
    return pairs, reference_segments, genome_names


def collect_cluster_genomes(pairs, reference_segments, workdir):
    unique_pairs = workdir / "cluster_genome.unique.tsv"
    cluster_genomes = workdir / "cluster_genomes.tsv"
    reference_clusters = set()
    with reference_segments.open() as source:
        for line in source:
            reference_clusters.add(line.split("\t", 1)[0])

    external_sort(pairs, unique_pairs, unique=True)
    with unique_pairs.open() as source, cluster_genomes.open("w") as output:
        previous_cluster = None
        genomes = []

        def write_cluster():
            if previous_cluster in reference_clusters:
                output.write(previous_cluster + "\t" + "\t".join(genomes) + "\n")

        for line in source:
            cluster_id, genome = line.rstrip("\n").split("\t", 1)
            if previous_cluster is not None and cluster_id != previous_cluster:
                write_cluster()
                genomes = []
            previous_cluster = cluster_id
            genomes.append(genome)
        if previous_cluster is not None:
            write_cluster()
    return cluster_genomes


def load_reference_segments(reference_segments, cluster_genomes, genome_names, include_nonchromosomal):
    members_by_cluster = {}
    with cluster_genomes.open() as source:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            members_by_cluster[fields[0]] = tuple(genome_names[genome] for genome in fields[1:])

    by_sequence = defaultdict(list)
    with reference_segments.open() as source:
        for line in source:
            cluster_id, sequence_name, start, end = line.rstrip("\n").split("\t")
            if include_nonchromosomal or not NON_CHROMOSOMAL_RE.search(sequence_name):
                by_sequence[sequence_name].append((int(start), int(end), members_by_cluster[cluster_id]))
    return dict(sorted(by_sequence.items(), key=lambda item: chromosome_key(item[0])))


def plot_segments(by_sequence, genome_names, reference, output, dpi):
    if not by_sequence:
        raise ValueError(f"no reference segments found for {reference!r}")
    labels = sorted(set(genome_names.values()), key=lambda name: (TIP_ORDER.index(name) if name in TIP_ORDER else len(TIP_ORDER), name.casefold()))
    figure, axis = plt.subplots(figsize=(16, max(7.5, 0.14 * len(labels) + 2)), layout="constrained")
    row_by_label = {label: index for index, label in enumerate(labels)}
    sequence_lengths = {
        name: max(end for _, end, _ in segments)
        for name, segments in by_sequence.items()
    }
    gap = max(100_000, max(sequence_lengths.values()) // 100)
    total_length = sum(sequence_lengths.values()) + gap * (len(sequence_lengths) - 1)
    bin_size = max(1, (total_length + 3999) // 4000)
    matrix = np.zeros((len(labels), (total_length + bin_size - 1) // bin_size), dtype=np.uint8)
    sequence_midpoints = []
    offset = 0

    for sequence_name, segments in by_sequence.items():
        sequence_length = sequence_lengths[sequence_name]
        sequence_midpoints.append((offset + sequence_length / 2, sequence_name))
        for start, end, members in segments:
            for label in members:
                left = (offset + start) // bin_size
                right = max(left + 1, (offset + end + bin_size - 1) // bin_size)
                matrix[row_by_label[label], left:right] = 1
        offset += sequence_length + gap

    image = axis.imshow(
        matrix,
        cmap=ListedColormap(["#f0f3f2", "#167d78"]),
        interpolation="nearest",
        aspect="auto",
        extent=(0, total_length, len(labels) - 0.5, -0.5),
        vmin=0,
        vmax=1,
    )
    axis.set_yticks(range(len(labels)), labels=labels, fontsize=7)
    axis.set_xticks([position for position, _ in sequence_midpoints], [name for _, name in sequence_midpoints], rotation=45, ha="right", fontsize=6)
    for midpoint, sequence_name in sequence_midpoints[:-1]:
        axis.axvline(midpoint + sequence_lengths[sequence_name] / 2, color="#ffffff", linewidth=1.2)
    axis.set_xlabel("Concatenated reference sequence position (bp; color = segment present)")
    figure.colorbar(image, ax=axis, ticks=[0.25, 0.75], label="Shared segment")
    figure.axes[-1].set_yticklabels(["Absent", "Present"])
    figure.savefig(output, dpi=dpi)
    plt.close(figure)


def main():
    args = parse_args()
    if not args.bed.is_file():
        raise SystemExit(f"BED file not found: {args.bed}")
    if args.metadata is not None and not args.metadata.is_file():
        raise SystemExit(f"metadata file not found: {args.metadata}")
    if shutil.which("sort") is None:
        raise SystemExit("required command not found: sort")

    reference = resolve_reference(args.reference, args.metadata)
    metadata_names = load_genome_metadata(args.metadata)
    output = args.output or args.bed.with_name(f"{args.bed.stem}.{reference.replace(' ', '_')}.conservation.svg")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="conservation_plot_") as temporary:
        workdir = Path(temporary)
        pairs, reference_segments, genome_names = collect_records(args.bed, reference, workdir, metadata_names)
        cluster_genomes = collect_cluster_genomes(pairs, reference_segments, workdir)
        by_sequence = load_reference_segments(reference_segments, cluster_genomes, genome_names, args.all_reference_sequences)
        plot_segments(by_sequence, genome_names, reference, output, args.dpi)
    print(output)


if __name__ == "__main__":
    main()