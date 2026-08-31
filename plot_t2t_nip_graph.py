#!/usr/bin/env python3
"""Render text-free visual summaries of the T2T Nipponbare clusters."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET


PALETTE = ["#f28e2b", "#159bd0", "#2ca25f", "#8e5aa9"]
NUCLEAR_PREFIX = "NC_089"


def fasta_lengths(path: Path) -> list[tuple[str, int]]:
    lengths: list[tuple[str, int]] = []
    name = None
    length = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths.append((name, length))
                name = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        lengths.append((name, length))
    return lengths


def read_bed(path: Path) -> list[tuple[str, int, int, int]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, cluster = line.rstrip().split()[:4]
            if chrom.rsplit("-", 1)[-1].startswith(NUCLEAR_PREFIX):
                rows.append((chrom, int(start), int(end), int(cluster)))
    return rows


def largest_clusters(rows: list[tuple[str, int, int, int]]) -> list[tuple[int, int]]:
    members = Counter(cluster for _, _, _, cluster in rows)
    return sorted(members.items(), key=lambda item: (-item[1], item[0]))[:4]


def svg_root(width: int, height: int) -> ET.Element:
    return ET.Element(
        "svg",
        {
            "xmlns": "http://www.w3.org/2000/svg",
            "width": str(width),
            "height": str(height),
            "viewBox": f"0 0 {width} {height}",
            "role": "img",
        },
    )


def write_svg(root: ET.Element, path: Path) -> None:
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def chromosome_key(chrom: str) -> str:
    return chrom.rsplit("-", 1)[-1]


def write_cluster_fasta(
    fasta_path: Path,
    rows: list[tuple[str, int, int, int]],
    cluster: int,
    output: Path,
) -> int:
    targets: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for chrom, start, end, row_cluster in rows:
        if row_cluster == cluster:
            targets[chromosome_key(chrom)].append((start, end, row_cluster))
    for loci in targets.values():
        loci.sort()

    output.parent.mkdir(parents=True, exist_ok=True)
    record_count = 0
    current_chrom = ""
    sequence: list[str] = []
    with fasta_path.open() as source, output.open("w") as destination:
        def emit() -> None:
            nonlocal record_count
            for start, end, _ in targets.get(current_chrom, []):
                subsequence = "".join(sequence)[start:end]
                record_count += 1
                destination.write(
                    f">cluster{cluster}_locus{record_count:03d}|"
                    f"{current_chrom}:{start}-{end}\n"
                )
                for offset in range(0, len(subsequence), 80):
                    destination.write(subsequence[offset:offset + 80] + "\n")

        for line in source:
            if line.startswith(">"):
                if current_chrom:
                    emit()
                current_chrom = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
        if current_chrom:
            emit()
    return record_count


def make_track_svg(
    rows: list[tuple[str, int, int, int]], lengths: list[tuple[str, int]], path: Path
) -> None:
    length_map = {name: length for name, length in lengths}
    selected = largest_clusters(rows)
    selected_colors = {cluster: PALETTE[rank] for rank, (cluster, _) in enumerate(selected)}
    by_chrom: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for chrom, start, end, cluster in rows:
        name = chromosome_key(chrom)
        if cluster in selected_colors and name in length_map:
            by_chrom[name].append((start, end, cluster))

    ordered = [(name, length) for name, length in lengths]
    width = 1800
    left = 40
    right = 40
    track_width = width - left - right
    top = 32
    row_height = 34
    track_height = 16
    height = top + row_height * len(ordered) + 24
    root = svg_root(width, height)
    ET.SubElement(root, "rect", {"width": str(width), "height": "100%", "fill": "#f7f8f5"})

    for row_index, (name, length) in enumerate(ordered):
        y = top + row_index * row_height
        ET.SubElement(
            root,
            "rect",
            {
                "x": str(left),
                "y": str(y),
                "width": str(track_width),
                "height": str(track_height),
                "rx": "8",
                "fill": "#dfe3dc",
            },
        )
        for start, end, cluster in sorted(by_chrom[name]):
            x = left + track_width * start / length
            segment_width = max(1.0, track_width * (end - start) / length)
            ET.SubElement(
                root,
                "rect",
                {
                    "x": f"{x:.2f}",
                    "y": str(y),
                    "width": f"{segment_width:.2f}",
                    "height": str(track_height),
                    "fill": selected_colors[cluster],
                    "fill-opacity": "0.82",
                },
            )

    write_svg(root, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, default=Path("data/t2t_nip.fasta"))
    parser.add_argument("--bed", type=Path, default=Path("t2t_nip.seg.bed"))
    parser.add_argument("--tracks", type=Path, default=Path("t2t_nip_chromosomes.svg"))
    parser.add_argument("--cluster-fasta", type=Path, default=None)
    args = parser.parse_args()

    lengths = [item for item in fasta_lengths(args.fasta)
               if item[0].startswith(NUCLEAR_PREFIX)]
    rows = read_bed(args.bed)
    make_track_svg(rows, lengths, args.tracks)
    selected = largest_clusters(rows)
    top_cluster, top_count = selected[0]
    cluster_fasta = args.cluster_fasta or Path("data") / f"cluster{top_cluster}.fasta"
    extracted = write_cluster_fasta(args.fasta, rows, top_cluster, cluster_fasta)
    print(f"clusters={len(selected)}")
    print(f"loci={sum(count for cluster, count in selected)}")
    print(f"selected={','.join(str(cluster) for cluster, _ in selected)}")
    print(f"tracks={args.tracks}")
    print(f"top_cluster={top_cluster}")
    print(f"top_cluster_loci={top_count}")
    print(f"cluster_fasta={cluster_fasta}")
    print(f"cluster_fasta_records={extracted}")


if __name__ == "__main__":
    main()
