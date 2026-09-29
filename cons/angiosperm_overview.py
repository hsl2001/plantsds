# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8", "networkx>=3.2"]
# ///

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib
import networkx as nx

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_overview(bed_path):
    cluster_genomes = defaultdict(set)
    genome_segments = Counter()
    row_count = 0

    with bed_path.open(encoding="utf-8") as bed_file:
        for line_number, line in enumerate(bed_file, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{bed_path}:{line_number}: expected at least 4 BED columns"
                )

            genome, separator, _ = fields[0].rpartition("-")
            if not separator or not genome:
                raise ValueError(
                    f"{bed_path}:{line_number}: expected genome-sequence chrom name"
                )
            cluster_id = fields[3]
            cluster_genomes[cluster_id].add(genome)
            genome_segments[genome] += 1
            row_count += 1

    if not row_count:
        raise ValueError(f"No BED records found in {bed_path}")

    return row_count, cluster_genomes, genome_segments


def load_genome_orders(genomes, metadata_path, taxonomy_path):
    accession_to_species = {
        row["accession"]: row["species"]
        for row in csv.DictReader(metadata_path.open(encoding="utf-8"), delimiter="\t")
    }
    species_orders = {}
    genus_orders = defaultdict(set)
    with taxonomy_path.open(encoding="utf-8") as taxonomy_file:
        for row in csv.DictReader(taxonomy_file, delimiter="\t"):
            lineage = [
                part.split("|", 2)
                for part in row["taxid_rank_name"].split(";")
                if len(part.split("|", 2)) == 3
            ]
            order = next((name for _, rank, name in lineage if rank == "order"), None)
            genus = next((name for _, rank, name in lineage if rank == "genus"), None)
            if order:
                species_orders[row["scientific_name"]] = order
                if genus:
                    genus_orders[genus].add(order)

    # Older genome records sometimes use names absent from the lineage table.
    species_order_overrides = {
        "Thellungiella parvula": "Brassicales",
        "Streptochaeta angustifolia": "Poales",
        "Mesua ferrea": "Malpighiales",
        "Pennisetum glaucum": "Poales",
    }
    genome_orders = {}
    for genome in genomes:
        match = re.match(r"((?:GCA|GCF)_\d+\.\d+)", genome)
        accession = match.group(1) if match else ""
        species = accession_to_species.get(accession, "")
        order = species_orders.get(species) or species_order_overrides.get(species)
        if not order and species:
            genus_candidates = genus_orders.get(species.split()[0], set())
            if len(genus_candidates) == 1:
                order = next(iter(genus_candidates))
        genome_orders[genome] = order
    return genome_orders


def plot_overview(bed_path, output_prefix, metadata_path, taxonomy_path):
    row_count, cluster_genomes, genome_segments = load_overview(bed_path)
    genome_count = len(genome_segments)
    genome_orders = load_genome_orders(
        genome_segments, metadata_path, taxonomy_path
    )
    graph = nx.Graph()
    for genome, segment_count in genome_segments.items():
        graph.add_node(genome, segments=segment_count)
    for genomes in cluster_genomes.values():
        for left, right in combinations(sorted(genomes), 2):
            if graph.has_edge(left, right):
                graph[left][right]["weight"] += 1
            else:
                graph.add_edge(left, right, weight=1)

    for _, _, data in graph.edges(data=True):
        data["layout_weight"] = data["weight"] ** 0.45

    positions = nx.spring_layout(
        graph,
        k=0.62,
        iterations=300,
        weight="layout_weight",
        seed=42,
    )
    order_colors = {
        "Brassicales": "#00a8e8",
        "Solanales": "#f58220",
        "Poales": "#1fa667",
        "Fagales": "#9b59b6",
    }
    other_color = "#888888"

    background = "#ffffff"
    figure = plt.figure(figsize=(12, 12), facecolor=background)
    network_axis = figure.add_axes([0, 0, 1, 1], facecolor=background)

    edges = list(graph.edges(data=True))

    def draw_edges(edge_records):
        max_weight = max(data["weight"] for _, _, data in edges)
        widths = [
            0.8 + 2.2 * (data["weight"] / max_weight) ** 0.35
            for _, _, data in edge_records
        ]
        nx.draw_networkx_edges(
            graph,
            positions,
            edgelist=[(left, right) for left, right, _ in edge_records],
            width=widths,
            edge_color="#000000",
            alpha=1.0,
            ax=network_axis,
        )

    draw_edges(edges)

    nx.draw_networkx_nodes(
        graph,
        positions,
        node_size=150,
        node_color=[order_colors.get(genome_orders[genome], other_color) for genome in graph],
        edgecolors="black",
        linewidths=1.2,
        ax=network_axis,
    )

    network_axis.set_aspect("equal")
    network_axis.margins(0.04)
    network_axis.axis("off")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = output_prefix.with_suffix(".svg")
    png_path = output_prefix.with_suffix(".png")
    figure.savefig(svg_path, facecolor=figure.get_facecolor(), pad_inches=0)
    figure.savefig(png_path, dpi=220, facecolor=figure.get_facecolor(), pad_inches=0)
    plt.close(figure)
    print(f"Wrote {svg_path} and {png_path}", file=sys.stderr)
    print(
        f"{row_count:,} segments, {genome_count:,} genomes, "
        f"{len(cluster_genomes):,} clusters, {graph.number_of_edges():,} links, "
        f"{sum(order in order_colors for order in genome_orders.values())} genomes "
        "in the four requested orders",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Draw a genome network from shared clusters in a SegTrace BED."
    )
    parser.add_argument(
        "-b",
        "--bed",
        type=Path,
        default=Path(__file__).with_name("ANGIOSPERM.seg.bed"),
        help="input SegTrace BED",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=Path,
        default=Path(__file__).parent.parent
        / "conservation_plot"
        / "angiosperm_genome_overview",
        help="output path prefix for SVG and PNG",
    )
    parser.add_argument(
        "--genomes",
        type=Path,
        default=Path(__file__).with_name("angio_wgd_genomes.tsv"),
        help="genome accession-to-species metadata TSV",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        default=Path(__file__).with_name("taxonomy_rank_lineage.tsv"),
        help="taxonomy lineage TSV",
    )
    arguments = parser.parse_args()

    try:
        plot_overview(
            arguments.bed,
            arguments.output_prefix,
            arguments.genomes,
            arguments.taxonomy,
        )
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()