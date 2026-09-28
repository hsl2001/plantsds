# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8"]
# ///

import argparse
from array import array
import csv
from collections import defaultdict
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile

VALID_BASE_RUN = re.compile(rb"[ACGT]+")


def fasta_label(path):
    name = path.name
    for suffix in (".gz", ".bgz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for suffix in (".fa", ".fna", ".fasta", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def load_fastas(selected_dir):
    list_path = next(
        (
            path
            for path in (
                selected_dir / "angio_wgd_genomes.files",
                selected_dir / "genomes.files",
            )
            if path.is_file()
        ),
        None,
    )
    if list_path is None:
        raise FileNotFoundError(
            f"FASTA list not found in {selected_dir}: expected "
            "angio_wgd_genomes.files or genomes.files"
        )

    fastas = {}
    for line in list_path.read_text(encoding="utf-8").splitlines():
        value = line.strip()
        if not value:
            continue
        path = Path(value)
        if not path.is_absolute():
            path = Path.cwd() / path
        if not path.is_file():
            raise FileNotFoundError(f"FASTA not found: {path}")
        path = path.resolve()
        label = fasta_label(path)
        if label in fastas:
            previous = fastas[label]
            if path.samefile(previous):
                print(f"[input] Ignoring repeated FASTA path: {path}", file=sys.stderr)
                continue
            raise ValueError(
                f"Ambiguous FASTA label {label!r}: {previous} and {path}. "
                "BED chrom names cannot distinguish these inputs."
            )
        fastas[label] = path
    if not fastas:
        raise ValueError(f"No FASTA paths in {list_path}")
    return fastas


def find_arabidopsis_labels(selected_dir, fastas):
    summary_path = next(
        (
            path
            for path in (
                selected_dir / "angio_wgd_genomes.tsv",
                selected_dir / "genomes.tsv",
            )
            if path.is_file()
        ),
        None,
    )
    if summary_path is None:
        raise FileNotFoundError(
            f"Species manifest not found in {selected_dir}: expected "
            "angio_wgd_genomes.tsv or genomes.tsv"
        )
    with summary_path.open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("species", "").strip().casefold() == "arabidopsis thaliana"
        ]
    if not rows:
        raise ValueError(f"No Arabidopsis thaliana rows in {summary_path}")
    labels = []
    for row in rows:
        label = fasta_label(Path(row["used_fasta"]))
        if label not in fastas:
            raise ValueError(
                f"Arabidopsis FASTA from {summary_path} is absent from the FASTA list: {label}"
            )
        labels.append(label)
    return tuple(dict.fromkeys(labels))


def split_chrom(chrom, labels):
    prefix = ""
    match = None
    for part in chrom.split("-")[:-1]:
        prefix = part if not prefix else f"{prefix}-{part}"
        if prefix in labels:
            match = (prefix, chrom[len(prefix) + 1 :])
    if match is None or not match[1]:
        raise ValueError(f"BED chrom does not match a selected FASTA label: {chrom}")
    return match


def build_fai(fasta_path, fai_path):
    records = []
    names = set()
    current = None

    def finish_record():
        if current is None:
            return
        if current[2] == 0:
            raise ValueError(f"FASTA record has no sequence: {current[0]} in {fasta_path}")
        records.append((current[0], current[2], current[1], current[3], current[4]))

    with fasta_path.open("rb") as fasta:
        for line in fasta:
            if line.startswith(b">"):
                finish_record()
                fields = line[1:].split(None, 1)
                if not fields:
                    raise ValueError(f"Empty FASTA header in {fasta_path}")
                name = fields[0].decode("utf-8", errors="replace")
                if name in names:
                    raise ValueError(f"Duplicate FASTA record {name!r} in {fasta_path}")
                names.add(name)
                current = [name, fasta.tell(), 0, None, None, None, None]
                continue

            if current is None:
                if line.strip():
                    raise ValueError(f"Sequence data before FASTA header in {fasta_path}")
                continue
            bases = line.rstrip(b"\r\n")
            if not bases:
                raise ValueError(f"Blank sequence line in {fasta_path}")
            width = len(line)
            if current[3] is None:
                current[3] = len(bases)
                current[4] = width
            elif current[5] != current[3] or current[6] != current[4]:
                raise ValueError(
                    f"Irregular FASTA line wrapping in {fasta_path}; reformat it or create a valid .fai"
                )
            current[2] += len(bases)
            current[5] = len(bases)
            current[6] = width
    finish_record()

    temporary = fai_path.with_name(fai_path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as output:
        for record in records:
            output.write("\t".join(map(str, record)) + "\n")
    os.replace(temporary, fai_path)


def load_fai(fasta_path):
    fai_path = Path(str(fasta_path) + ".fai")
    if not fai_path.is_file() or fai_path.stat().st_mtime_ns < fasta_path.stat().st_mtime_ns:
        print(f"[index] Building {fai_path}", file=sys.stderr, flush=True)
        build_fai(fasta_path, fai_path)

    index = {}
    with fai_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed FASTA index: {fai_path}:{line_number}")
            name = fields[0]
            try:
                index[name] = tuple(map(int, fields[1:5]))
            except ValueError as error:
                raise ValueError(f"Malformed FASTA index: {fai_path}:{line_number}") from error
    return index


def interval_composition(fasta, record_index, start, end, record_name, fasta_path):
    length, offset, line_bases, line_width = record_index
    if start < 0 or end <= start or end > length:
        raise ValueError(
            f"BED interval outside {fasta_path}:{record_name}: {start}-{end}, length {length}"
        )

    first_line = start // line_bases
    last_base = end - 1
    last_line = last_base // line_bases
    byte_start = offset + first_line * line_width + start % line_bases
    byte_end = offset + last_line * line_width + last_base % line_bases + 1
    fasta.seek(byte_start)
    sequence = fasta.read(byte_end - byte_start).translate(None, b"\r\n").upper()
    if len(sequence) != end - start:
        raise ValueError(f"FASTA index does not match sequence data: {fasta_path}:{record_name}")

    c_count = sequence.count(b"C")
    g_count = sequence.count(b"G")
    valid_bases = sum(sequence.count(base) for base in b"ACGT")
    if not valid_bases:
        return 0, 0, 0, 0.0

    valid_pairs = sum(
        len(run) - 1 for run in VALID_BASE_RUN.findall(sequence) if len(run) > 1
    )
    cpg_observed = sequence.count(b"CG")
    cpg_expected = valid_pairs * c_count * g_count / (valid_bases * valid_bases)
    return c_count + g_count, valid_bases, cpg_observed, cpg_expected


def index_bed(bed_path, fastas, database_path):
    labels = set(fastas)
    connection = sqlite3.connect(database_path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-65536")
    connection.execute(
        "CREATE TABLE loci (genome TEXT, seq TEXT, start INTEGER, end INTEGER, "
        "cluster TEXT, status TEXT, UNIQUE(genome, seq, start, end, cluster, status))"
    )

    count = 0
    batch = []
    connection.execute("BEGIN")
    with bed_path.open(encoding="utf-8") as bed:
        for line_number, line in enumerate(bed, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"{bed_path}:{line_number}: expected at least 4 BED columns")
            genome, seq = split_chrom(fields[0], labels)
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{bed_path}:{line_number}: invalid BED coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{bed_path}:{line_number}: invalid BED interval {start}-{end}")
            status = fields[4] if len(fields) >= 5 else ""
            batch.append((genome, seq, start, end, fields[3], status))
            count += 1
            if len(batch) >= 100_000:
                connection.executemany(
                    "INSERT OR IGNORE INTO loci VALUES (?, ?, ?, ?, ?, ?)", batch
                )
                batch.clear()
    if batch:
        connection.executemany(
            "INSERT OR IGNORE INTO loci VALUES (?, ?, ?, ?, ?, ?)", batch
        )
    connection.commit()
    if count == 0:
        connection.close()
        raise ValueError(f"No BED intervals found in {bed_path}")
    connection.execute("CREATE INDEX loci_order ON loci(genome, seq, start, end, cluster)")
    connection.commit()
    unique_count = connection.execute("SELECT COUNT(*) FROM loci").fetchone()[0]
    return connection, count, unique_count


def collect_cluster_features(connection, genomes, fastas, trajectory_genomes, output, stats):
    connection.execute(
        "CREATE TABLE cluster_composition ("
        "cluster TEXT PRIMARY KEY, gc_bases INTEGER, valid_bases INTEGER, "
        "cpg_observed INTEGER, cpg_expected REAL, segment_count INTEGER, "
        "genome_count INTEGER)"
    )

    for genome_index, genome in enumerate(genomes, start=1):
        fasta_path = fastas[genome]
        fai = load_fai(fasta_path)
        totals = defaultdict(lambda: [0, 0, 0, 0.0, 0])
        query = connection.execute(
            "SELECT seq, start, end, cluster FROM loci "
            "WHERE genome = ? ORDER BY seq, start, end, cluster",
            (genome,),
        )
        print(
            f"[composition] {genome_index}/{len(genomes)} {genome}",
            file=sys.stderr,
            flush=True,
        )
        with fasta_path.open("rb") as fasta:
            for seq, start, end, cluster in query:
                if seq not in fai:
                    raise ValueError(f"Sequence {genome}-{seq} is absent from {fasta_path}")
                gc_bases, valid_bases, cpg_observed, cpg_expected = interval_composition(
                    fasta, fai[seq], start, end, seq, fasta_path
                )
                stats["segments_processed"] += 1
                if not valid_bases:
                    stats["segments_without_valid_bases"] += 1
                cluster_total = totals[cluster]
                cluster_total[0] += gc_bases
                cluster_total[1] += valid_bases
                cluster_total[2] += cpg_observed
                cluster_total[3] += cpg_expected
                cluster_total[4] += 1

        connection.executemany(
            "INSERT INTO cluster_composition VALUES (?, ?, ?, ?, ?, ?, 1) "
            "ON CONFLICT(cluster) DO UPDATE SET "
            "gc_bases = cluster_composition.gc_bases + excluded.gc_bases, "
            "valid_bases = cluster_composition.valid_bases + excluded.valid_bases, "
            "cpg_observed = cluster_composition.cpg_observed + excluded.cpg_observed, "
            "cpg_expected = cluster_composition.cpg_expected + excluded.cpg_expected, "
            "segment_count = cluster_composition.segment_count + excluded.segment_count, "
            "genome_count = cluster_composition.genome_count + 1",
            ((cluster, *values) for cluster, values in totals.items()),
        )
        connection.commit()

    arabidopsis_clusters = {
        genome: {
            row[0]
            for row in connection.execute(
                "SELECT DISTINCT cluster FROM loci WHERE genome = ?", (genome,)
            )
        }
        for genome in trajectory_genomes
    }
    x_values = array("f")
    y_values = array("f")
    cluster_ids = []
    query = connection.execute(
        "SELECT cluster, gc_bases, valid_bases, cpg_observed, cpg_expected, "
        "segment_count, genome_count FROM cluster_composition "
        "ORDER BY CAST(cluster AS INTEGER)"
    )
    for cluster, gc_bases, valid_bases, cpg_observed, cpg_expected, segment_count, genome_count in query:
        stats["clusters_total"] += 1
        gc_percent = 100.0 * gc_bases / valid_bases if valid_bases else None
        cpg_oe = cpg_observed / cpg_expected if cpg_expected > 0 else None
        accessions = [
            genome for genome in trajectory_genomes
            if cluster in arabidopsis_clusters[genome]
        ]
        output.write(
            f"{cluster}\t{genome_count}\t{segment_count}\t{valid_bases}\t"
            f"{gc_bases}\t{cpg_observed}\t{cpg_expected:.9f}\t"
            f"{'' if gc_percent is None else f'{gc_percent:.6f}'}\t"
            f"{'' if cpg_oe is None else f'{cpg_oe:.9f}'}\t"
            f"{';'.join(accessions)}\n"
        )
        if gc_percent is None or cpg_oe is None:
            stats["clusters_without_defined_cpg_oe"] += 1
            continue
        x_values.append(gc_percent)
        y_values.append(cpg_oe)
        cluster_ids.append(cluster)
        stats["clusters_plotted"] += 1

    return x_values, y_values, cluster_ids, arabidopsis_clusters


def make_plots(
    output_prefix,
    x_values,
    y_values,
    cluster_ids,
    arabidopsis_clusters,
    genome_count,
    trajectory_genomes,
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    x = np.frombuffer(x_values, dtype=np.float32)
    y = np.frombuffer(y_values, dtype=np.float32)
    max_cpg_oe = max(float(np.max(y)), 1.0)
    outputs = []
    for include_trajectory, suffix in (
        (False, "_background"),
        (True, "_with_arabidopsis"),
    ):
        figure, axis = plt.subplots(figsize=(8, 7), layout="constrained")
        bins = axis.hexbin(
            x,
            y,
            gridsize=140,
            extent=(0, 75, 0, max_cpg_oe),
            mincnt=1,
            bins="log",
            cmap="viridis",
            linewidths=0,
        )
        if include_trajectory:
            from matplotlib.lines import Line2D

            trajectory_colors = ("#d81b8a", "#f28e2b")
            markers = ("o", "^", "s", "D")
            handles = []
            for index, genome in enumerate(trajectory_genomes):
                color = trajectory_colors[index % len(trajectory_colors)]
                indices = [
                    point_index
                    for point_index, cluster in enumerate(cluster_ids)
                    if cluster in arabidopsis_clusters[genome]
                ]
                if indices:
                    axis.scatter(
                        x[indices],
                        y[indices],
                        s=13,
                        marker=markers[index % len(markers)],
                        color=color,
                        edgecolors="white",
                        linewidths=0.3,
                        alpha=0.9,
                        zorder=5,
                    )
                label = genome.rsplit("_", 1)[-1].removesuffix(".nuclear")
                handles.append(
                    Line2D(
                        [0], [0], color=color, marker=markers[index % len(markers)],
                        markersize=5, linewidth=0, label=label
                    )
                )
            axis.legend(
                handles=handles,
                loc="upper right",
                framealpha=0.9,
            )
        axis.axhline(0, color="#555555", linewidth=0.8, alpha=0.7)
        axis.set(
            xlim=(0, 75),
            ylim=(0, max_cpg_oe),
            xlabel="Cluster GC (%)",
            ylabel="CpG O/E (observed / expected)",
            title=(
                "Cluster composition + Arabidopsis thaliana membership"
                if include_trajectory
                else "Selected genome cluster composition"
            ),
        )
        axis.grid(color="#d9e2e1", linewidth=0.5, alpha=0.5)
        colorbar = figure.colorbar(bins, ax=axis)
        colorbar.set_label("Cluster count (log scale)")
        figure.suptitle(
            f"{genome_count} genomes; {len(x):,} clusters; "
            "cluster statistics pooled across all member segments"
        )
        output_path = Path(str(output_prefix) + suffix + ".png")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=180)
        plt.close(figure)
        outputs.append(output_path)
    return outputs


def run(arguments):
    selected_dir = arguments.selected_dir.resolve()
    bed_path = arguments.bed.resolve()
    output_prefix = arguments.output_prefix.resolve()
    if not bed_path.is_file():
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    if not selected_dir.is_dir():
        raise FileNotFoundError(f"Selected directory not found: {selected_dir}")
    fastas = load_fastas(selected_dir)
    trajectory_genomes = find_arabidopsis_labels(selected_dir, fastas)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    scratch = os.environ.get("TMPDIR") or None
    with tempfile.TemporaryDirectory(prefix="gc_transition_", dir=scratch) as temporary_dir:
        database_path = Path(temporary_dir) / "intervals.sqlite"
        connection, bed_count, unique_locus_count = index_bed(
            bed_path, fastas, database_path
        )
        genomes = [row[0] for row in connection.execute("SELECT DISTINCT genome FROM loci ORDER BY genome")]
        missing = [genome for genome in genomes if genome not in fastas]
        if missing:
            raise ValueError(f"BED genomes missing from FASTA list: {missing[:5]}")
        missing_trajectory = [genome for genome in trajectory_genomes if genome not in genomes]
        if missing_trajectory:
            raise ValueError(
                f"Arabidopsis thaliana has no intervals in {bed_path}: {missing_trajectory}"
            )

        output_path = Path(str(output_prefix) + ".tsv")
        temporary_output = output_path.with_name(output_path.name + ".partial")
        stats = {
            "segments_processed": 0,
            "segments_without_valid_bases": 0,
            "clusters_total": 0,
            "clusters_plotted": 0,
            "clusters_without_defined_cpg_oe": 0,
        }
        with temporary_output.open("w", encoding="utf-8") as output:
            output.write(
                "cluster_id\tgenome_count\tsegment_count\tvalid_bases\tgc_bases\t"
                "cpg_observed\tcpg_expected\tgc_percent\tcpg_oe\t"
                "arabidopsis_accessions\n"
            )
            x_values, y_values, cluster_ids, arabidopsis_clusters = collect_cluster_features(
                connection, genomes, fastas, trajectory_genomes, output, stats
            )
        connection.close()
        if stats["clusters_plotted"] == 0:
            temporary_output.unlink(missing_ok=True)
            raise ValueError("No clusters have both valid GC and defined CpG O/E")
        os.replace(temporary_output, output_path)

    plot_paths = make_plots(
        output_prefix,
        x_values,
        y_values,
        cluster_ids,
        arabidopsis_clusters,
        len(genomes),
        trajectory_genomes,
    )
    clusters_with_arabidopsis = sum(
        any(cluster in arabidopsis_clusters[genome] for genome in trajectory_genomes)
        for cluster in cluster_ids
    )
    print(
        f"BED records={bed_count:,}; unique loci={unique_locus_count:,}; "
        f"segments={stats['segments_processed']:,}; "
        f"clusters={stats['clusters_total']:,}; plotted={stats['clusters_plotted']:,}; "
        f"without_defined_CpG_OE={stats['clusters_without_defined_cpg_oe']:,}; "
        f"segments_without_valid_bases={stats['segments_without_valid_bases']:,}; "
        f"Arabidopsis={','.join(trajectory_genomes)}; "
        f"clusters_with_Arabidopsis="
        f"{clusters_with_arabidopsis:,}",
        file=sys.stderr,
    )
    print(f"Wrote {output_path}, " + ", ".join(map(str, plot_paths)), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Plot pooled cluster GC percentage against CpG observed/expected."
    )
    parser.add_argument("--selected-dir", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        run(arguments)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
