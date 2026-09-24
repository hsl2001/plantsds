#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.collections import LineCollection
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator


ACCESSION_RE = re.compile(r"GC[AF]_\d+\.\d+")
GENOME_RE = re.compile(
    r"^(.+?)-(?:chr|chromosome|contig|scaffold|unplaced|seq\d+|tig\d+|ptg\d+)",
    re.IGNORECASE,
)
NON_CHROMOSOMAL_RE = re.compile(r"(?:contig|scaffold|unplaced)", re.IGNORECASE)
CHROMOSOME_NUMBER_RE = re.compile(r"(?:chr|chromosome)[_-]?(\d+|[xy])", re.IGNORECASE)
GENOME_ALIASES = {"tair12": "Col-0"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot the number of unique genomes represented by each reference segment cluster."
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


def collect_records(bed, reference, workdir):
    pairs = workdir / "cluster_genome.tsv"
    reference_segments = workdir / "reference_segments.tsv"
    excluded = excluded_genomes(reference)
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
            if genome_id(sequence_name) not in excluded:
                pair_output.write(f"{cluster_id}\t{genome_id(sequence_name)}\n")
            if is_reference_sequence(sequence_name, reference):
                reference_output.write(f"{cluster_id}\t{sequence_name}\t{start}\t{end}\n")
    return pairs, reference_segments


def count_genomes(pairs, workdir):
    unique_pairs = workdir / "cluster_genome.unique.tsv"
    cluster_counts = workdir / "cluster_counts.tsv"
    external_sort(pairs, unique_pairs, unique=True)
    with unique_pairs.open() as source, cluster_counts.open("w") as output:
        previous_cluster = None
        count = 0
        for line in source:
            cluster_id = line.split("\t", 1)[0]
            if previous_cluster is not None and cluster_id != previous_cluster:
                output.write(f"{previous_cluster}\t{count}\n")
                count = 0
            previous_cluster = cluster_id
            count += 1
        if previous_cluster is not None:
            output.write(f"{previous_cluster}\t{count}\n")
    return cluster_counts


def load_reference_segments(reference_segments, cluster_counts, include_nonchromosomal):
    counts = {}
    with cluster_counts.open() as source:
        for line in source:
            cluster_id, count = line.rstrip("\n").split("\t")
            counts[cluster_id] = int(count)

    by_sequence = defaultdict(list)
    with reference_segments.open() as source:
        for line in source:
            cluster_id, sequence_name, start, end = line.rstrip("\n").split("\t")
            if include_nonchromosomal or not NON_CHROMOSOMAL_RE.search(sequence_name):
                by_sequence[sequence_name].append((int(start), int(end), counts.get(cluster_id, 0)))
    return dict(sorted(by_sequence.items(), key=lambda item: chromosome_key(item[0])))


def plot_segments(by_sequence, reference, output, dpi):
    if not by_sequence:
        raise ValueError(f"no reference segments found for {reference!r}")
    n_panels = len(by_sequence)
    figure, axes = plt.subplots(
        n_panels,
        1,
        figsize=(12, max(2.6, 2.1 * n_panels)),
        squeeze=False,
        layout="constrained",
    )
    maximum_count = max(count for segments in by_sequence.values() for _, _, count in segments)
    norm = Normalize(vmin=0, vmax=max(1, maximum_count))
    cmap = plt.get_cmap("viridis")

    for axis, (sequence_name, segments) in zip(axes[:, 0], by_sequence.items()):
        lines = [[(start, count), (end, count)] for start, end, count in segments]
        fills = [[(start, 0), (start, count), (end, count), (end, 0)] for start, end, count in segments]
        colors = cmap(norm([count for _, _, count in segments]))
        axis.add_collection(PolyCollection(fills, facecolors=colors, edgecolors="none", alpha=0.6))
        axis.add_collection(LineCollection(lines, colors=colors, linewidths=0.35, alpha=0.9))
        axis.set_xlim(0, max(end for _, end, _ in segments))
        axis.set_ylim(0, maximum_count + 0.5)
        axis.yaxis.set_major_locator(MaxNLocator(integer=True))
        axis.set_ylabel("Genomes")
        axis.text(0.01, 0.95, sequence_name, transform=axis.transAxes, ha="left", va="top", fontsize=8)
    axes[-1, 0].set_xlabel("Position (bp)")
    figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes[:, 0], label="Genomes", pad=0.01)
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
    output = args.output or args.bed.with_name(f"{args.bed.stem}.{reference.replace(' ', '_')}.conservation.svg")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="conservation_plot_") as temporary:
        workdir = Path(temporary)
        pairs, reference_segments = collect_records(args.bed, reference, workdir)
        cluster_counts = count_genomes(pairs, workdir)
        by_sequence = load_reference_segments(reference_segments, cluster_counts, args.all_reference_sequences)
        plot_segments(by_sequence, reference, output, args.dpi)
    print(output)


if __name__ == "__main__":
    main()