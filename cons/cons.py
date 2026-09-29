# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8"]
# ///

import argparse
import csv
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

def top_conserved_segments(reference_segments, cluster_genomes, limit=5):
    return sorted(
        reference_segments,
        key=lambda segment: (
            -len(cluster_genomes[segment[3]]),
            -(segment[2] - segment[1]),
            natural_key(segment[0]),
            segment[1],
            segment[2],
            segment[3],
        ),
    )[:limit]

def write_overlapping_features(gff_path, segments, cluster_genomes, output_path):
    ranked_segments = [
        (rank, segment, len(cluster_genomes[segment[3]]))
        for rank, segment in enumerate(segments, start=1)
    ]
    feature_count = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gff_path.open(encoding="utf-8") as gff_file, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as output_file:
        writer = csv.writer(output_file, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "segment_rank",
                "segment_seqid",
                "segment_start",
                "segment_end",
                "conservation_genomes",
                "cluster_id",
                "feature_seqid",
                "feature_source",
                "feature_type",
                "feature_start",
                "feature_end",
                "feature_score",
                "feature_strand",
                "feature_phase",
                "feature_attributes",
                "overlap_bp",
            ]
        )

        for line_number, line in enumerate(gff_file, start=1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    f"{gff_path}:{line_number}: expected 9 GFF3 columns"
                )
            try:
                feature_start = int(fields[3]) - 1
                feature_end = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"{gff_path}:{line_number}: GFF start/end must be integers"
                ) from error
            if feature_start < 0 or feature_end <= feature_start:
                raise ValueError(f"{gff_path}:{line_number}: invalid GFF interval")

            for rank, segment, conservation in ranked_segments:
                segment_seqid, segment_start, segment_end, cluster_id = segment
                if segment_seqid != fields[0]:
                    continue
                overlap_start = max(segment_start, feature_start)
                overlap_end = min(segment_end, feature_end)
                if overlap_start >= overlap_end:
                    continue
                writer.writerow(
                    [
                        rank,
                        segment_seqid,
                        segment_start + 1,
                        segment_end,
                        conservation,
                        cluster_id,
                        *fields[:3],
                        feature_start + 1,
                        feature_end,
                        *fields[5:9],
                        overlap_end - overlap_start,
                    ]
                )
                feature_count += 1

    return feature_count

def plot_conservation(reference, bed_path, output_path):
    reference_segments, cluster_genomes = load_reference_segments(
        bed_path, reference
    )
    events_by_chrom, chrom_lengths = make_events(
        reference_segments, cluster_genomes
    )
    max_reference_conservation = max(
        (len(cluster_genomes[cluster_id]) for _, _, _, cluster_id in reference_segments),
        default=0,
    )

    chromosomes = sorted(events_by_chrom, key=natural_key)
    figure, axes = plt.subplots(
        1,
        len(chromosomes),
        figsize=(max(14, 3.2 * len(chromosomes)), 4.5),
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
        chrom_length = chrom_lengths[chrom]
        axis.set_xlim(0, chrom_length / 1_000_000)
        figure.canvas.draw()
        axis.bar(
            coordinates[:-1],
            counts[:-1],
            width=bar_widths,
            align="edge",
            color="#187c72",
            linewidth=0,
            rasterized=True,
        )
        narrow_centers = []
        narrow_heights = []
        for left, width, count in zip(coordinates, bar_widths, counts):
            left_px = axis.transData.transform((left, 0))[0]
            right_px = axis.transData.transform((left + width, 0))[0]
            if abs(right_px - left_px) < 1:
                narrow_centers.append(left + width / 2)
                narrow_heights.append(count)
        if narrow_centers:
            axis.vlines(
                narrow_centers,
                0,
                narrow_heights,
                color="#187c72",
                linewidth=0.8,
                rasterized=True,
            )
        axis.set_title(chrom, loc="left", fontsize=10)
        axis.set_xlabel("Coordinate (Mb)")
        axis.yaxis.set_major_locator(MaxNLocator(integer=True))
        axis.set_ylim(bottom=0)
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
        f"{max_reference_conservation} genomes max on reference)",
        file=sys.stderr,
    )
    return reference_segments, cluster_genomes

def main():
    parser = argparse.ArgumentParser(
        description="Plot the number of genomes with segments across a reference."
    )
    parser.add_argument("-r", "--reference", required=True, help="reference accession")
    parser.add_argument(
        "-b", "--bed", type=Path, required=True, help="input SegTrace BED"
    )
    parser.add_argument(
        "-g", "--gff", type=Path, help="GFF3 annotation for overlapping features"
    )
    parser.add_argument("-o", "--output", type=Path, help="output plot path")
    arguments = parser.parse_args()

    try:
        output_path = arguments.output or Path(
            "conservation_plot", f"{arguments.reference}_segment_count.svg"
        )
        reference_segments, cluster_genomes = plot_conservation(
            arguments.reference, arguments.bed, output_path
        )
        if arguments.gff:
            segments = top_conserved_segments(reference_segments, cluster_genomes)
            feature_output = output_path.with_name(
                f"{output_path.stem}_top5_features.tsv"
            )
            feature_count = write_overlapping_features(
                arguments.gff, segments, cluster_genomes, feature_output
            )
            print(
                f"Wrote {feature_output} from top {len(segments)} conserved segments "
                f"({feature_count} overlapping features)",
                file=sys.stderr,
            )
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))

if __name__ == "__main__":
    main()