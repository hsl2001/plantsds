# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8"]
# ///

import argparse
from array import array
import csv
from itertools import groupby
import os
from pathlib import Path
import sqlite3
import sys
import tempfile


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
    list_path = selected_dir / "angio_wgd_genomes.files"
    if not list_path.is_file():
        raise FileNotFoundError(f"FASTA list not found: {list_path}")

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


def find_arabidopsis_label(selected_dir, fastas):
    summary_path = selected_dir / "angio_wgd_genomes.tsv"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Species manifest not found: {summary_path}")
    with summary_path.open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("species", "").strip().casefold() == "arabidopsis thaliana"
        ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one Arabidopsis thaliana row in {summary_path}, found {len(rows)}"
        )
    label = fasta_label(Path(rows[0]["used_fasta"]))
    if label not in fastas:
        raise ValueError(f"Arabidopsis FASTA from {summary_path} is absent from the FASTA list: {label}")
    return label


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


def interval_gc(fasta, record_index, start, end, record_name, fasta_path):
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
    sequence = fasta.read(byte_end - byte_start).translate(None, b"\r\n")
    if len(sequence) != end - start:
        raise ValueError(f"FASTA index does not match sequence data: {fasta_path}:{record_name}")

    gc = sequence.count(b"G") + sequence.count(b"g") + sequence.count(b"C") + sequence.count(b"c")
    valid = sum(
        sequence.count(base) + sequence.count(bytes((base + 32,)))
        for base in b"ACGT"
    )
    return (100.0 * gc / valid, valid) if valid else (None, 0)


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


def process_genome(
    connection,
    genome,
    fasta_path,
    max_gap_bp,
    output,
    x_values,
    y_values,
    trajectory_genome,
    trajectory_paths,
    trajectory_x_values,
    trajectory_y_values,
    stats,
):
    fai = load_fai(fasta_path)
    query = connection.execute(
        "SELECT seq, start, end, cluster, status FROM loci "
        "WHERE genome = ? ORDER BY seq, start, end, cluster",
        (genome,),
    )
    with fasta_path.open("rb") as fasta:
        for seq, grouped in groupby(query, key=lambda row: row[0]):
            intervals = list(grouped)
            if seq not in fai:
                raise ValueError(f"Sequence {genome}-{seq} is absent from {fasta_path}")
            gc_intervals = []
            path_x = []
            path_y = []

            def finish_path():
                if path_x:
                    trajectory_paths.append((path_x.copy(), path_y.copy()))
                    path_x.clear()
                    path_y.clear()

            for _, start, end, cluster, status in intervals:
                gc, valid_bases = interval_gc(fasta, fai[seq], start, end, seq, fasta_path)
                gc_intervals.append((start, end, cluster, status, gc, valid_bases))
                stats["intervals"] += 1
                stats["invalid_gc"] += gc is None

            for index in range(1, len(gc_intervals) - 1):
                previous = gc_intervals[index - 1]
                current = gc_intervals[index]
                following = gc_intervals[index + 1]
                if previous[4] is None or current[4] is None or following[4] is None:
                    if genome == trajectory_genome:
                        finish_path()
                    continue
                if previous[1] > current[0] or current[1] > following[0]:
                    stats["overlap_skipped"] += 1
                    if genome == trajectory_genome:
                        finish_path()
                    continue
                previous_gap = current[0] - previous[1]
                following_gap = following[0] - current[1]
                if max_gap_bp and max(previous_gap, following_gap) > max_gap_bp:
                    stats["gap_skipped"] += 1
                    if genome == trajectory_genome:
                        finish_path()
                    continue

                delta_previous = current[4] - previous[4]
                delta_next = following[4] - current[4]
                output.write(
                    f"{genome}\t{seq}\t{current[2]}\t{current[0]}\t{current[1]}\t"
                    f"{current[4]:.6f}\t{previous[4]:.6f}\t{following[4]:.6f}\t"
                    f"{delta_previous:.6f}\t{delta_next:.6f}\t{previous_gap}\t"
                    f"{following_gap}\t{current[3]}\n"
                )
                x_values.append(delta_previous)
                y_values.append(delta_next)
                stats["points"] += 1
                if genome == trajectory_genome:
                    path_x.append(delta_previous)
                    path_y.append(delta_next)
                    trajectory_x_values.append(delta_previous)
                    trajectory_y_values.append(delta_next)
                    stats["trajectory_points"] += 1
            if genome == trajectory_genome:
                finish_path()


def make_plots(
    output_prefix,
    x_values,
    y_values,
    trajectory_paths,
    trajectory_x_values,
    trajectory_y_values,
    max_gap_bp,
    genome_count,
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    x = np.frombuffer(x_values, dtype=np.float32)
    y = np.frombuffer(y_values, dtype=np.float32)
    trajectory_x = np.frombuffer(trajectory_x_values, dtype=np.float32)
    trajectory_y = np.frombuffer(trajectory_y_values, dtype=np.float32)
    gap_note = "unlimited" if not max_gap_bp else f"<= {max_gap_bp:,} bp"
    line_paths = [
        np.column_stack((path_x, path_y))
        for path_x, path_y in trajectory_paths
        if len(path_x) >= 2
    ]
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
            extent=(-100, 100, -100, 100),
            mincnt=1,
            bins="log",
            cmap="viridis",
            linewidths=0,
        )
        if include_trajectory:
            from matplotlib.collections import LineCollection
            from matplotlib.lines import Line2D

            if line_paths:
                axis.add_collection(
                    LineCollection(
                        line_paths,
                        colors="#ffe082",
                        linewidths=0.75,
                        alpha=0.95,
                        zorder=4,
                    )
                )
            axis.scatter(
                trajectory_x,
                trajectory_y,
                s=5,
                color="#fff3b0",
                edgecolors="#263238",
                linewidths=0.15,
                zorder=5,
            )
            axis.legend(
                handles=[
                    Line2D(
                        [0], [0], color="#ffe082", marker="o", markersize=4,
                        linewidth=1, label="Arabidopsis thaliana trajectory"
                    )
                ],
                loc="upper right",
                framealpha=0.9,
            )
        axis.axhline(0, color="#555555", linewidth=0.8, alpha=0.7)
        axis.axvline(0, color="#555555", linewidth=0.8, alpha=0.7)
        axis.set(
            xlim=(-100, 100),
            ylim=(-100, 100),
            xlabel="GC change from previous segment (percentage points)",
            ylabel="GC change to next segment (percentage points)",
            title=(
                "Angiosperm background + Arabidopsis thaliana trajectory"
                if include_trajectory
                else "Angiosperm background"
            ),
        )
        axis.set_aspect("equal", adjustable="box")
        axis.grid(color="#d9e2e1", linewidth=0.5, alpha=0.5)
        colorbar = figure.colorbar(bins, ax=axis)
        colorbar.set_label("Interval occurrences (log count)")
        figure.suptitle(
            f"{genome_count} genomes; {len(x):,} background centers; "
            f"neighbor gaps {gap_note}"
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
    if arguments.max_gap_bp < 0:
        raise ValueError("max-gap-bp must be zero or positive")

    fastas = load_fastas(selected_dir)
    trajectory_genome = find_arabidopsis_label(selected_dir, fastas)
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
            raise ValueError(f"BED genomes missing from {selected_dir}/angio_wgd_genomes.files: {missing[:5]}")
        if trajectory_genome not in genomes:
            raise ValueError(f"Arabidopsis thaliana has no intervals in {bed_path}: {trajectory_genome}")

        output_path = Path(str(output_prefix) + ".tsv")
        temporary_output = output_path.with_name(output_path.name + ".partial")
        x_values = array("f")
        y_values = array("f")
        trajectory_paths = []
        trajectory_x_values = array("f")
        trajectory_y_values = array("f")
        stats = {
            "intervals": 0,
            "invalid_gc": 0,
            "overlap_skipped": 0,
            "gap_skipped": 0,
            "points": 0,
            "trajectory_points": 0,
        }
        with temporary_output.open("w", encoding="utf-8") as output:
            output.write(
                "genome\tseq\tcluster_id\tstart\tend\tgc_current\tgc_previous\tgc_next\t"
                "delta_previous\tdelta_next\tgap_previous_bp\tgap_next_bp\tstatus\n"
            )
            for index, genome in enumerate(genomes, start=1):
                print(f"[gc] {index}/{len(genomes)} {genome}", file=sys.stderr, flush=True)
                process_genome(
                    connection,
                    genome,
                    fastas[genome],
                    arguments.max_gap_bp,
                    output,
                    x_values,
                    y_values,
                    trajectory_genome,
                    trajectory_paths,
                    trajectory_x_values,
                    trajectory_y_values,
                    stats,
                )
        connection.close()
        if stats["points"] == 0 or stats["trajectory_points"] == 0:
            temporary_output.unlink(missing_ok=True)
            raise ValueError(
                "No valid Angiosperm background or Arabidopsis GC transitions found"
            )
        os.replace(temporary_output, output_path)

    plot_paths = make_plots(
        output_prefix,
        x_values,
        y_values,
        trajectory_paths,
        trajectory_x_values,
        trajectory_y_values,
        arguments.max_gap_bp,
        len(genomes),
    )
    print(
        f"BED records={bed_count:,}; unique loci={unique_locus_count:,}; "
        f"processed={stats['intervals']:,}; "
        f"points={stats['points']:,}; invalid_GC={stats['invalid_gc']:,}; "
        f"overlap_skipped={stats['overlap_skipped']:,}; "
        f"gap_skipped={stats['gap_skipped']:,}; "
        f"Arabidopsis={trajectory_genome}; "
        f"trajectory_points={stats['trajectory_points']:,}",
        file=sys.stderr,
    )
    print(f"Wrote {output_path}, " + ", ".join(map(str, plot_paths)), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Build a GC transition table and density plot from SegTrace BED intervals."
    )
    parser.add_argument("--selected-dir", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument(
        "--max-gap-bp",
        type=int,
        default=0,
        help="maximum gap on either side (0 means unlimited; default: 0)",
    )
    arguments = parser.parse_args()
    try:
        run(arguments)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
