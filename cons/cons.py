# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8"]
# ///

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.ticker import MaxNLocator

def natural_key(value):
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", value)
    )

def load_reference_segments(bed_path, reference):
    cluster_genomes = defaultdict(set)
    reference_segments = []

    with bed_path.open(encoding="utf-8") as bed_file:
        for line_number, line in enumerate(bed_file, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{bed_path}:{line_number}: expected at least 4 BED columns"
                )

            chrom = fields[0]
            genome, separator, sequence = chrom.rpartition("-")
            if not separator:
                raise ValueError(
                    f"{bed_path}:{line_number}: expected genome-sequence chrom name"
                )

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{bed_path}:{line_number}: BED start/end must be integers"
                ) from error
            if start < 0 or end <= start:
                raise ValueError(f"{bed_path}:{line_number}: invalid BED interval")

            cluster_id = fields[3]
            cluster_genomes[cluster_id].add(genome)
            if reference in genome:
                if sequence.startswith(reference + "_"):
                    sequence = sequence[len(reference) + 1 :]
                reference_segments.append((sequence, start, end, cluster_id))

    if not reference_segments:
        raise ValueError(
            f"No BED records for reference {reference!r} in {bed_path}"
        )

    return reference_segments, cluster_genomes

def make_events(reference_segments, cluster_genomes):
    events_by_chrom = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    chrom_lengths = defaultdict(int)

    for chrom, start, end, cluster_id in reference_segments:
        genomes = cluster_genomes[cluster_id]
        chrom_lengths[chrom] = max(chrom_lengths[chrom], end)
        for genome in genomes:
            events_by_chrom[chrom][start][genome] += 1
            events_by_chrom[chrom][end][genome] -= 1

    return events_by_chrom, chrom_lengths

def plot_conservation(reference, bed_path, output_path):
    reference_segments, cluster_genomes = load_reference_segments(
        bed_path, reference
    )
    events_by_chrom, chrom_lengths = make_events(
        reference_segments, cluster_genomes
    )

    chromosomes = sorted(events_by_chrom, key=natural_key)
    figure, axes = plt.subplots(
        len(chromosomes),
        1,
        figsize=(14, max(3, 2.6 * len(chromosomes))),
        sharey=True,
        squeeze=False,
        layout="constrained",
    )

    for chrom, axis in zip(chromosomes, axes.flat):
        active_genomes = defaultdict(int)
        coordinates = []
        counts = []

        for position in sorted(events_by_chrom[chrom]):
            for genome, delta in events_by_chrom[chrom][position].items():
                active_genomes[genome] += delta
                if active_genomes[genome] == 0:
                    del active_genomes[genome]
            coordinates.append(position / 1_000_000)
            counts.append(len(active_genomes))

        bar_widths = [
            right - left for left, right in zip(coordinates, coordinates[1:])
        ]
        axis.bar(
            coordinates[:-1],
            counts[:-1],
            width=bar_widths,
            align="edge",
            color="#187c72",
            linewidth=0,
        )
        chrom_length = chrom_lengths[chrom]
        axis.set_xlim(0, chrom_length / 1_000_000)
        axis.set_title(chrom, loc="left", fontsize=10)
        axis.set_xlabel("Coordinate (Mb)")
        axis.yaxis.set_major_locator(MaxNLocator(integer=True))
        axis.set_ylim(bottom=0)
        axis.grid(axis="y", color="#d9e2e1", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)

    figure.supylabel("Genomes with a detected segment")
    figure.suptitle(f"Segment conservation on {reference}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)
    print(
        f"Wrote {output_path} from {bed_path} "
        f"({len(reference_segments)} reference segments, "
        f"{len(chrom_lengths)} chromosomes, "
        f"{max(map(len, cluster_genomes.values()))} genomes max)",
        file=sys.stderr,
    )

def main():
    parser = argparse.ArgumentParser(
        description="Plot the number of genomes with segments across a reference."
    )
    parser.add_argument("-r", "--reference", required=True, help="reference accession")
    parser.add_argument(
        "-b", "--bed", type=Path, required=True, help="input SegTrace BED"
    )
    parser.add_argument("-o", "--output", type=Path, help="output plot path")
    arguments = parser.parse_args()

    try:
        output_path = arguments.output or Path(
            "conservation_plot", f"{arguments.reference}_segment_count.png"
        )
        plot_conservation(arguments.reference, arguments.bed, output_path)
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))

if __name__ == "__main__":
    main()