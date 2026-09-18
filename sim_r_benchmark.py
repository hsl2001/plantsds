#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy>=1.26", "edlib>=1.2.7"]
# ///
"""Validate SegTrace read mode on Badread-simulated long reads.

The genome defaults intentionally match sim_benchmark.py. Badread read
headers contain the source contig, strand, and source interval, so those
simulation truth fields are used instead of remapping the reads.
"""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import sim_benchmark as sb


DEFAULT_SIMULATION = {
    "species": 3,
    "chromosomes": 3,
    "chrom_length": 1_000_000,
    "fragments": 100,
    "min_fragment_length": 1_000,
    "max_fragment_length": 50_000,
    "min_copies": 2,
    "max_copies": 10,
    "max_snp_rate": 0.10,
    "max_indel_rate": 0.01,
    "filler_chunk": 1_000_000,
}
DEFAULT_QUANTITIES = ("1x", "2x", "4x", "8x", "16x")


@dataclass(frozen=True)
class BadreadTruth:
    source_start: int
    source_end: int
    strand: str
    target: str
    read_length: int


def executable(name: str) -> str:
    path = Path(name)
    if path.exists():
        return str(path.resolve())
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(f"executable not found: {name}")


def run_logged(command: list[str], log_path: Path, stdout_path: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_path.open("w") as log_handle:
            if stdout_path is None:
                subprocess.run(
                    command,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=log_handle,
                )
            else:
                with stdout_path.open("w") as stdout_handle:
                    subprocess.run(
                        command,
                        check=True,
                        stdout=stdout_handle,
                        stderr=log_handle,
                    )
    except subprocess.CalledProcessError as error:
        stderr = log_path.read_text(errors="replace").strip()
        detail = stderr[-4000:] if stderr else "(no stderr output)"
        raise RuntimeError(
            f"command failed with exit status {error.returncode}: "
            f"{shlex.join(command)}\n{detail}"
        ) from error


def build_simulation_args(out_dir: Path, seed: int, force: bool) -> argparse.Namespace:
    values = dict(DEFAULT_SIMULATION)
    values.update(out_dir=str(out_dir), seed=seed, force=force)
    return argparse.Namespace(**values)


def parse_badread_truth(path: Path) -> dict[str, BadreadTruth]:
    """Read source coordinates embedded in Badread FASTQ headers."""
    truth: dict[str, BadreadTruth] = {}
    with path.open() as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            handle.readline()
            handle.readline()
            if not header.startswith("@"):
                continue
            fields = header[1:].strip().split()
            if len(fields) < 2:
                continue
            source = fields[1].split(",")
            if len(source) < 3 or source[1] not in {"+strand", "-strand"}:
                continue
            try:
                source_start, source_end = map(int, source[2].split("-", 1))
            except ValueError:
                continue
            truth[fields[0]] = BadreadTruth(
                source_start=source_start,
                source_end=source_end,
                strand="+" if source[1] == "+strand" else "-",
                target=source[0],
                read_length=len(sequence),
            )
    return truth


def lift_prediction(
    prediction: Path,
    truth: dict[str, BadreadTruth],
    output: Path,
    read_label: str,
) -> tuple[int, int]:
    """Project read-coordinate BED intervals using Badread source truth."""
    mapped = 0
    unmapped = 0
    prefix = f"{read_label}-"
    with prediction.open() as source, output.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                read_start = int(fields[1])
                read_end = int(fields[2])
            except ValueError:
                continue
            if not fields[0].startswith(prefix):
                unmapped += 1
                continue
            read_name = fields[0][len(prefix):]
            read_truth = truth.get(read_name)
            if read_truth is None or read_truth.read_length <= 0:
                unmapped += 1
                continue

            left = max(0, min(read_start, read_truth.read_length))
            right = max(0, min(read_end, read_truth.read_length))
            if right <= left:
                unmapped += 1
                continue
            source_length = read_truth.source_end - read_truth.source_start
            ref_left = round(left * source_length / read_truth.read_length)
            ref_right = round(right * source_length / read_truth.read_length)
            if read_truth.strand == "+":
                ref_start = read_truth.source_start + ref_left
                ref_end = read_truth.source_start + ref_right
            else:
                ref_start = read_truth.source_end - ref_right
                ref_end = read_truth.source_end - ref_left
            if ref_end <= ref_start:
                unmapped += 1
                continue
            writer.writerow([read_truth.target, ref_start, ref_end, fields[3] if len(fields) > 3 else ""])
            mapped += 1
    return mapped, unmapped


def run_segtrace(
    args: argparse.Namespace,
    reads: Path,
    output_dir: Path,
    mode: str,
) -> tuple[Path, Path, str]:
    prefix = output_dir / mode / "segtrace"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [
        executable(args.segtrace_bin),
        "-k", str(args.kmer),
        "-s", str(args.scale),
        "-w", str(args.window_size),
        "-c", str(args.min_report_copies),
        "-p", str(args.threads),
        "-o", str(prefix),
    ]
    if args.step_size:
        command.extend(["-t", str(args.step_size)])
    if mode == "segtrace_r":
        command.append("-r")
    command.append(str(reads))
    log_path = prefix.parent / "segtrace.stderr.txt"
    run_logged(command, log_path)
    prediction = prefix.with_suffix(".seg.bed")
    if not prediction.exists():
        raise FileNotFoundError(f"SegTrace did not create {prediction}")
    return prediction, log_path, shlex.join(command)


def extract_coverage(log_path: Path) -> float | None:
    pattern = re.compile(r"haploid coverage ~ ([0-9]+(?:\.[0-9]+)?)x")
    for line in log_path.read_text().splitlines():
        match = pattern.search(line)
        if match:
            return float(match.group(1))
    return None


def parse_quantities(value: str) -> list[str]:
    quantities = [item for item in re.split(r"[,\s:]+", value.strip()) if item]
    if not quantities:
        raise SystemExit("--quantities must contain at least one Badread quantity")
    return quantities


def write_results(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "quantity", "mode", "status", "raw_prediction", "mapped_prediction", "command",
        "haploid_coverage", "mapped_intervals", "unmapped_intervals",
        *sb.CSV_FIELDS,
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark SegTrace -r using Badread reads and sim_benchmark.py defaults."
    )
    parser.add_argument("--out-dir", default="sim_r_benchmark")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--quantities",
        default=",".join(DEFAULT_QUANTITIES),
        help="Badread quantities separated by commas, spaces, or colons",
    )
    parser.add_argument("--read-length", default="15000,13000", help="Badread mean,stdev")
    parser.add_argument("--identity", default="95,99,2.5", help="Badread identity distribution")
    parser.add_argument("--threads", "-p", type=int, default=8)
    parser.add_argument("--segtrace-bin", default="./segtrace")
    parser.add_argument("--badread-bin", default="badread")
    parser.add_argument("--kmer", "-k", type=int, default=17)
    parser.add_argument("--scale", "-s", type=int, default=16)
    parser.add_argument("--window-size", "-w", type=int, default=1024)
    parser.add_argument("--step-size", "-t", type=int, default=0)
    parser.add_argument("--min-report-copies", "-c", type=int, default=2)
    parser.add_argument("--skip-control", action="store_true", help="Skip the no--r control run")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be positive")
    if args.min_report_copies < 2:
        raise SystemExit("--min-report-copies must be >= 2 for one read input")

    out_dir = Path(args.out_dir)
    if out_dir.exists():
        if not args.force:
            raise SystemExit(f"{out_dir} exists; pass --force to replace it")
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    genome_dir = out_dir / "genome"
    simulation_args = build_simulation_args(genome_dir, args.seed, force=True)
    paths, placements = sb.generate_simulation(simulation_args)
    print(f"[INFO] Generated {len(paths.fasta_paths)} genomes and {len(placements)} truth intervals")

    rows: list[dict[str, object]] = []
    modes = ["segtrace_r"] if args.skip_control else ["segtrace_r", "segtrace_no_r"]
    for quantity in parse_quantities(args.quantities):
        quantity_dir = out_dir / quantity.replace("/", "_")
        quantity_dir.mkdir(parents=True, exist_ok=True)
        reads = quantity_dir / "badread.fastq"
        badread_command = [
            executable(args.badread_bin), "simulate",
            "--reference", str(paths.combined_fasta),
            "--quantity", quantity,
            "--length", args.read_length,
            "--identity", args.identity,
            "--error_model", "random",
            "--qscore_model", "ideal",
            "--seed", str(args.seed),
            "--glitches", "0,0,0",
            "--junk_reads", "0",
            "--random_reads", "0",
            "--chimeras", "0",
            "--start_adapter_seq", "",
            "--end_adapter_seq", "",
        ]
        run_logged(badread_command, quantity_dir / "badread.stderr.txt", reads)
        read_truth = parse_badread_truth(reads)
        print(f"[INFO] {quantity}: parsed Badread truth headers: {len(read_truth)}")

        for mode in modes:
            prediction, log_path, command_text = run_segtrace(
                args, reads, quantity_dir, mode
            )
            mapped_prediction = quantity_dir / mode / "segtrace.reference.bed"
            mapped_count, unmapped_count = lift_prediction(
                prediction, read_truth, mapped_prediction, reads.stem
            )
            metrics = sb.evaluate_with_numpy(mapped_prediction, paths.truth_bed)
            row: dict[str, object] = {
                "quantity": quantity,
                "mode": mode,
                "status": "ok",
                "raw_prediction": str(prediction),
                "mapped_prediction": str(mapped_prediction),
                "command": command_text,
                "haploid_coverage": extract_coverage(log_path) or "",
                "mapped_intervals": mapped_count,
                "unmapped_intervals": unmapped_count,
                **metrics,
            }
            rows.append(row)
            print(
                f"[RESULT] {quantity} {mode}: bp_f1={metrics['bp_f1']:.4f} "
                f"frag_f1={metrics['frag_f1']:.4f} "
                f"mapped={mapped_count} unmapped={unmapped_count}"
            )

    write_results(out_dir / "benchmark_r_results.csv", rows)
    print(f"[INFO] Truth BED: {paths.truth_bed}")
    print(f"[INFO] Results CSV: {out_dir / 'benchmark_r_results.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
