#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy>=1.26"]
# ///
"""
Simulation benchmark for duplicated-segment callers and whole-genome mappers.

Protocol:
- Generate genome FASTA files with 10 chromosomes of 100 Mb per species by default.
- Generate 100 random 1 kb to 50 kb source fragments.
- Draw 2 to 10 copies per fragment.
- Inject 0 to 10% SNPs and 0 to 1% INDELs into each copy.
- Place all copies into non-overlapping genomic loci and write a gold-standard BED.
- Run enabled tools with default-like parameters and evaluate BP/Frag accuracy with numpy.
"""

from __future__ import annotations

import argparse
import csv
import platform
import random
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np


BASES = b"ACGT"
BYTE_TO_BASE = bytes(BASES[value & 3] for value in range(256))
CSV_FIELDS = [
    "tool", "status", "enabled", "command", "prediction_bed", "error",
    "time_perf_seconds", "max_rss_kb",
    "pred_bp", "truth_bp", "intersect_bp", "bp_recall", "bp_precision", "bp_f1",
    "truth_fragments", "pred_fragments", "frag_tp", "frag_fp", "frag_fn", "frag_recall", "frag_precision", "frag_f1",
]


@dataclass(frozen=True)
class CopyPlacement:
    species: str
    chrom: str
    start: int
    end: int
    event_id: str
    copy_id: str
    snp_rate: float
    indel_rate: float
    sequence: bytes


@dataclass(frozen=True)
class SimulationPaths:
    fasta_paths: list[Path]
    combined_fasta: Path
    truth_bed: Path


@dataclass
class Profile:
    wall_seconds: float = 0.0
    max_rss_kb: int = 0
    perf: dict[str, float] = field(default_factory=dict)
    command: str = ""


@dataclass
class ToolResult:
    tool: str
    status: str
    enabled: bool = True
    prediction_bed: Path | None = None
    profile: Profile = field(default_factory=Profile)
    metrics: dict[str, float] = field(default_factory=dict)
    error: str = ""


@dataclass(frozen=True)
class ToolSpec:
    name: str
    no_flag: str
    cli_name: str
    description: str
    runner: Callable[[argparse.Namespace, SimulationPaths, Path], ToolResult]


class FastaWriter:
    def __init__(self, handle, width: int = 80):
        self.handle = handle
        self.width = width
        self.column = 0
        self.records: list[tuple[str, int, int, int, int]] = []
        self.current_name = ""
        self.current_length = 0
        self.current_offset = 0

    def header(self, name: str) -> None:
        self.finish_record()
        if self.column:
            self.handle.write(b"\n")
        self.handle.write(f">{name}\n".encode())
        self.current_name = name
        self.current_length = 0
        self.current_offset = self.handle.tell()
        self.column = 0

    def sequence(self, seq: bytes) -> None:
        self.current_length += len(seq)
        offset = 0
        while offset < len(seq):
            take = min(self.width - self.column, len(seq) - offset)
            self.handle.write(seq[offset:offset + take])
            self.column += take
            offset += take
            if self.column == self.width:
                self.handle.write(b"\n")
                self.column = 0

    def finish_record(self) -> None:
        if self.column:
            self.handle.write(b"\n")
            self.column = 0
        if self.current_name:
            self.records.append((self.current_name, self.current_length, self.current_offset, self.width, self.width + 1))
            self.current_name = ""

    def write_fai(self, fasta_path: Path) -> None:
        with fasta_path.with_suffix(fasta_path.suffix + ".fai").open("w") as handle:
            for record in self.records:
                handle.write("\t".join(str(value) for value in record) + "\n")
            offset += take
            if self.column == self.width:
                self.handle.write(b"\n")
                self.column = 0

    def finish_record(self) -> None:
        if self.column:
            self.handle.write(b"\n")
            self.column = 0
        if self.current_name:
            self.records.append((self.current_name, self.current_length, self.current_offset, self.width, self.width + 1))
            self.current_name = ""

    def write_fai(self, fasta_path: Path) -> None:
        with fasta_path.with_suffix(fasta_path.suffix + ".fai").open("w") as handle:
            for record in self.records:
                handle.write("\t".join(str(value) for value in record) + "\n")


def random_dna(rng: random.Random, length: int) -> bytes:
    return rng.randbytes(length).translate(BYTE_TO_BASE)


def mutate_copy(rng: random.Random, source: bytes, snp_rate: float, indel_rate: float) -> bytes:
    mutated = bytearray()
    for base in source:
        if rng.random() < indel_rate and rng.random() < 0.5:
            continue
        if rng.random() < indel_rate:
            mutated.append(rng.choice(BASES))
        if rng.random() < snp_rate:
            choices = [candidate for candidate in BASES if candidate != base]
            mutated.append(rng.choice(choices))
        else:
            mutated.append(base)
    return bytes(mutated)


def overlaps(existing: list[tuple[int, int]], start: int, end: int) -> bool:
    return any(max(start, old_start) < min(end, old_end) for old_start, old_end in existing)


def choose_locus(
    rng: random.Random,
    used: dict[tuple[str, str], list[tuple[int, int]]],
    species_names: list[str],
    chrom_names: list[str],
    chrom_length: int,
    copy_length: int,
    max_attempts: int = 10000,
) -> tuple[str, str, int, int]:
    if copy_length >= chrom_length:
        raise ValueError(f"copy length {copy_length} must be shorter than chromosome length {chrom_length}")

    for _ in range(max_attempts):
        species = rng.choice(species_names)
        chrom = rng.choice(chrom_names)
        start = rng.randint(0, chrom_length - copy_length)
        end = start + copy_length
        occupied = used[(species, chrom)]
        if not overlaps(occupied, start, end):
            occupied.append((start, end))
            return species, chrom, start, end

    raise RuntimeError("failed to place a non-overlapping copy; reduce copy count or increase genome size")


def generate_simulation(args: argparse.Namespace) -> tuple[SimulationPaths, list[CopyPlacement]]:
    rng = random.Random(args.seed)
    out_dir = Path(args.out_dir)
    fasta_dir = out_dir / "fasta"
    truth_bed = out_dir / "truth.bed"
    combined_fasta = out_dir / "combined.fa"

    if out_dir.exists():
        if not args.force:
            raise FileExistsError(f"{out_dir} already exists; pass --force to replace it")
        shutil.rmtree(out_dir)
    fasta_dir.mkdir(parents=True)

    species_names = [f"sp{idx:02d}" for idx in range(1, args.species + 1)]
    chrom_names = [f"chr{idx:02d}" for idx in range(1, args.chromosomes + 1)]
    used = {(species, chrom): [] for species in species_names for chrom in chrom_names}
    placements_by_locus: dict[tuple[str, str], list[CopyPlacement]] = {key: [] for key in used}
    placements: list[CopyPlacement] = []

    for event_idx in range(1, args.fragments + 1):
        event_id = f"frag{event_idx:04d}"
        fragment_length = rng.randint(args.min_fragment_length, args.max_fragment_length)
        source = random_dna(rng, fragment_length)
        copy_count = rng.randint(args.min_copies, args.max_copies)

        for copy_idx in range(1, copy_count + 1):
            snp_rate = rng.uniform(0.0, args.max_snp_rate)
            indel_rate = rng.uniform(0.0, args.max_indel_rate)
            sequence = mutate_copy(rng, source, snp_rate, indel_rate)
            if not sequence:
                sequence = source[:1]
            species, chrom, start, end = choose_locus(
                rng, used, species_names, chrom_names, args.chrom_length, len(sequence)
            )
            placement = CopyPlacement(
                species=species,
                chrom=chrom,
                start=start,
                end=end,
                event_id=event_id,
                copy_id=f"copy{copy_idx:02d}",
                snp_rate=snp_rate,
                indel_rate=indel_rate,
                sequence=sequence,
            )
            placements.append(placement)
            placements_by_locus[(species, chrom)].append(placement)

    fasta_paths: list[Path] = []
    filler_chunk = max(1, args.filler_chunk)
    with combined_fasta.open("wb") as combined_handle:
        combined_writer = FastaWriter(combined_handle)
        for species in species_names:
            fasta_path = fasta_dir / f"{species}.fa"
            fasta_paths.append(fasta_path)
            with fasta_path.open("wb") as species_handle:
                species_writer = FastaWriter(species_handle)
                for chrom in chrom_names:
                    species_writer.header(chrom)
                    combined_writer.header(f"{species}-{chrom}")
                    cursor = 0
                    chrom_placements = sorted(placements_by_locus[(species, chrom)], key=lambda item: item.start)
                    for placement in chrom_placements:
                        write_random_filler_many([species_writer, combined_writer], rng, placement.start - cursor, filler_chunk)
                        write_sequence_many([species_writer, combined_writer], placement.sequence)
                        cursor = placement.end
                    write_random_filler_many([species_writer, combined_writer], rng, args.chrom_length - cursor, filler_chunk)
                    species_writer.finish_record()
                    combined_writer.finish_record()
                    species_writer.write_fai(fasta_path)
                combined_writer.write_fai(combined_fasta)

    with truth_bed.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for placement in sorted(placements, key=lambda item: (item.species, item.chrom, item.start)):
            seq_name = f"{placement.species}-{placement.chrom}"
            writer.writerow([
                seq_name,
                placement.start,
                placement.end,
                placement.event_id,
                placement.copy_id,
                f"snp={placement.snp_rate:.6f};indel={placement.indel_rate:.6f}",
            ])

    return SimulationPaths(fasta_paths=fasta_paths, combined_fasta=combined_fasta, truth_bed=truth_bed), placements


def write_random_filler_many(writers: list[FastaWriter], rng: random.Random, length: int, chunk_size: int) -> None:
    remaining = length
    while remaining > 0:
        take = min(chunk_size, remaining)
        seq = random_dna(rng, take)
        write_sequence_many(writers, seq)
        remaining -= take


def write_sequence_many(writers: list[FastaWriter], seq: bytes) -> None:
    for writer in writers:
        writer.sequence(seq)


def parse_bed_like(path: Path) -> list[tuple[str, int, int]]:
    if not path.exists():
        return []
    intervals = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            add_interval(fields, 0, intervals)
            if len(fields) >= 6:
                add_interval(fields, 3, intervals)
    return sorted(set(intervals))


def add_interval(fields: list[str], offset: int, intervals: list[tuple[str, int, int]]) -> None:
    try:
        chrom = fields[offset]
        start = int(fields[offset + 1])
        end = int(fields[offset + 2])
    except (IndexError, ValueError):
        return
    if end > start:
        intervals.append((chrom, start, end))


def write_bed(intervals: list[tuple[str, int, int]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for chrom, start, end in sorted(set(intervals)):
            writer.writerow([chrom, start, end])


def group_intervals(intervals: list[tuple[str, int, int]]) -> dict[str, np.ndarray]:
    grouped: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in intervals:
        grouped.setdefault(chrom, []).append((start, end))
    return {
        chrom: np.asarray(sorted(values), dtype=np.int64)
        for chrom, values in grouped.items()
        if values
    }


def merge_array(intervals: np.ndarray) -> np.ndarray:
    if intervals.size == 0:
        return np.empty((0, 2), dtype=np.int64)

    ordered = intervals[np.argsort(intervals[:, 0], kind="mergesort")]
    merged = [ordered[0].copy()]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append(np.array([start, end], dtype=np.int64))
    return np.vstack(merged)


def merged_by_chrom(intervals: list[tuple[str, int, int]]) -> dict[str, np.ndarray]:
    return {chrom: merge_array(array) for chrom, array in group_intervals(intervals).items()}


def bp_length(grouped: dict[str, np.ndarray]) -> int:
    return int(sum(np.sum(array[:, 1] - array[:, 0]) for array in grouped.values() if array.size))


def bp_intersection(left: dict[str, np.ndarray], right: dict[str, np.ndarray]) -> int:
    total = 0
    for chrom in left.keys() & right.keys():
        left_array = left[chrom]
        right_array = right[chrom]
        i = j = 0
        while i < len(left_array) and j < len(right_array):
            overlap = min(left_array[i, 1], right_array[j, 1]) - max(left_array[i, 0], right_array[j, 0])
            if overlap > 0:
                total += int(overlap)
            if left_array[i, 1] < right_array[j, 1]:
                i += 1
            else:
                j += 1
    return total


def reciprocal_match_count(query: dict[str, np.ndarray], target: dict[str, np.ndarray], fraction: float = 0.5) -> int:
    matches = 0
    for chrom in query.keys() & target.keys():
        target_array = target[chrom]
        target_starts = target_array[:, 0]
        target_ends = target_array[:, 1]
        for q_start, q_end in query[chrom]:
            q_length = q_end - q_start
            if q_length <= 0:
                continue
            idx = np.searchsorted(target_starts, q_end, side="left") - 1
            while idx >= 0 and target_ends[idx] > q_start:
                t_start, t_end = target_array[idx]
                t_length = t_end - t_start
                overlap = min(q_end, t_end) - max(q_start, t_start)
                if overlap > 0 and overlap / q_length >= fraction and overlap / t_length >= fraction:
                    matches += 1
                    break
                idx -= 1
    return matches


def evaluate_with_numpy(pred_bed: Path, truth_bed: Path) -> dict[str, float]:
    pred_intervals = parse_bed_like(pred_bed)
    truth_intervals = parse_bed_like(truth_bed)
    pred = group_intervals(pred_intervals)
    truth = group_intervals(truth_intervals)
    pred_merged = merged_by_chrom(pred_intervals)
    truth_merged = merged_by_chrom(truth_intervals)

    pred_bp = bp_length(pred_merged)
    truth_bp = bp_length(truth_merged)
    intersect_bp = bp_intersection(pred_merged, truth_merged)
    bp_recall = intersect_bp / truth_bp if truth_bp else 0.0
    bp_precision = intersect_bp / pred_bp if pred_bp else 0.0
    bp_f1 = f1_score(bp_recall, bp_precision)

    true_positive = reciprocal_match_count(truth, pred, fraction=0.5)
    matched_pred = reciprocal_match_count(pred, truth, fraction=0.5)
    truth_count = len(truth_intervals)
    pred_count = len(pred_intervals)
    frag_recall = true_positive / truth_count if truth_count else 0.0
    frag_precision = matched_pred / pred_count if pred_count else 0.0
    frag_f1 = f1_score(frag_recall, frag_precision)

    return {
        "pred_bp": pred_bp,
        "truth_bp": truth_bp,
        "intersect_bp": intersect_bp,
        "bp_recall": bp_recall,
        "bp_precision": bp_precision,
        "bp_f1": bp_f1,
        "truth_fragments": truth_count,
        "pred_fragments": pred_count,
        "frag_tp": true_positive,
        "frag_fp": pred_count - matched_pred,
        "frag_fn": truth_count - true_positive,
        "frag_recall": frag_recall,
        "frag_precision": frag_precision,
        "frag_f1": frag_f1,
    }


def f1_score(recall: float, precision: float) -> float:
    return 2.0 * recall * precision / (recall + precision) if recall + precision else 0.0


def resolve_executable(name: str) -> Path | None:
    path = Path(name)
    if path.exists():
        return path.resolve()
    found = shutil.which(name)
    return Path(found) if found else None


def profile_command(
    command: list[str],
    work_dir: Path,
    label: str,
    stdout_path: Path | None = None,
) -> Profile:
    profile_dir = work_dir / "profiles"
    profile_dir.mkdir(exist_ok=True)
    time_path = profile_dir / f"{label}.time.txt"
    stderr_path = profile_dir / f"{label}.stderr.txt"

    time_bin = Path("/usr/bin/time")
    time_flag = "-l" if platform.system() == "Darwin" else "-v"
    wrapped = [str(time_bin), time_flag, "-o", str(time_path), *command] if time_bin.exists() else command

    started = time.perf_counter()
    with stderr_path.open("w") as stderr_handle:
        if stdout_path is None:
            result = subprocess.run(wrapped, stdout=subprocess.DEVNULL, stderr=stderr_handle, check=False)
        else:
            with stdout_path.open("w") as stdout_handle:
                result = subprocess.run(wrapped, stdout=stdout_handle, stderr=stderr_handle, check=False)
    elapsed = time.perf_counter() - started

    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, wrapped)

    return Profile(
        wall_seconds=elapsed,
        max_rss_kb=parse_max_rss(time_path),
        command=" ".join(command),
    )


def parse_wall_seconds(path: Path, fallback: float) -> float:
    if not path.exists():
        return fallback
    text = path.read_text(errors="replace")
    match = re.search(r"Elapsed \(wall clock\) time .*: (\S+)", text)
    if not match:
        return fallback
    return parse_elapsed(match.group(1))


def parse_elapsed(value: str) -> float:
    parts = value.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(value)


def parse_max_rss(path: Path) -> int:
    if not path.exists():
        return 0
    text = path.read_text(errors="replace")
    linux_match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if linux_match:
        return int(linux_match.group(1))

    mac_match = re.search(r"^\s*(\d+)\s+maximum resident set size\s*$", text, re.IGNORECASE | re.MULTILINE)
    return int(mac_match.group(1)) // 1024 if mac_match else 0


def empty_bed(path: Path) -> Path:
    path.write_text("")
    return path


def tool_missing(tool: str, executable: str) -> ToolResult:
    return ToolResult(tool=tool, status="missing", error=f"{executable} not found in PATH")


def tool_failed(tool: str, error: Exception) -> ToolResult:
    return ToolResult(tool=tool, status="failed", error=str(error))


def result_from_prediction(tool: str, prediction_bed: Path, truth_bed: Path, profile: Profile) -> ToolResult:
    return ToolResult(
        tool=tool,
        status="ok",
        prediction_bed=prediction_bed,
        profile=profile,
        metrics=evaluate_with_numpy(prediction_bed, truth_bed),
    )


def run_truth(args: argparse.Namespace, paths: SimulationPaths, work_dir: Path) -> ToolResult:
    del args, work_dir
    return ToolResult(
        tool="Truth",
        status="ok",
        prediction_bed=paths.truth_bed,
        metrics=evaluate_with_numpy(paths.truth_bed, paths.truth_bed),
    )


def run_segtrace(args: argparse.Namespace, paths: SimulationPaths, work_dir: Path) -> ToolResult:
    tool = "SegTrace"
    segtrace_bin = resolve_executable(args.segtrace_bin)
    if segtrace_bin is None:
        return tool_missing(tool, args.segtrace_bin)
    prefix = work_dir / "segtrace"
    prediction = prefix.with_suffix(".dup.bed")
    command = [
        str(segtrace_bin), "-k", str(args.kmer), "-s", str(args.scale), "-w", str(args.window_size),
        "-t", str(args.step_size), "-c", str(args.min_report_copies), "-p", str(args.threads),
        "-o", str(prefix), *[str(path) for path in paths.fasta_paths],
    ]
    try:
        profile = profile_command(command, work_dir, "segtrace")
        if not prediction.exists():
            prediction = empty_bed(work_dir / "segtrace.empty.bed")
        return result_from_prediction(tool, prediction, paths.truth_bed, profile)
    except Exception as exc:
        return tool_failed(tool, exc)

def run_minimap2(args: argparse.Namespace, paths: SimulationPaths, work_dir: Path) -> ToolResult:
    tool = "minimap2"
    minimap_bin = resolve_executable(args.minimap2_bin)
    if minimap_bin is None:
        return tool_missing(tool, args.minimap2_bin)
    paf = work_dir / "minimap2.paf"
    prediction = work_dir / "minimap2.bed"
    command = [str(minimap_bin), "-x", "asm5", "-P", "-N", str(args.max_mappings), "-t", str(args.threads), str(paths.combined_fasta), str(paths.combined_fasta)]
    try:
        profile = profile_command(command, work_dir, "minimap2", stdout_path=paf)
        paf_to_bed(paf, prediction, args.min_call_length)
        return result_from_prediction(tool, prediction, paths.truth_bed, profile)
    except Exception as exc:
        return tool_failed(tool, exc)


def run_nucmer(args: argparse.Namespace, paths: SimulationPaths, work_dir: Path) -> ToolResult:
    tool = "MUMmer/nucmer"
    nucmer_bin = resolve_executable(args.nucmer_bin)
    show_coords_bin = resolve_executable(args.show_coords_bin)
    if nucmer_bin is None:
        return tool_missing(tool, args.nucmer_bin)
    if show_coords_bin is None:
        return tool_missing(tool, args.show_coords_bin)
    prefix = work_dir / "nucmer"
    coords = work_dir / "nucmer.coords"
    prediction = work_dir / "nucmer.bed"
    command = [str(nucmer_bin), "--maxmatch", "--prefix", str(prefix), str(paths.combined_fasta), str(paths.combined_fasta)]
    try:
        profile = profile_command(command, work_dir, "nucmer")
        with coords.open("w") as handle:
            subprocess.run([str(show_coords_bin), "-THrd", str(prefix) + ".delta"], stdout=handle, check=True)
        coords_to_bed(coords, prediction, args.min_call_length)
        return result_from_prediction(tool, prediction, paths.truth_bed, profile)
    except Exception as exc:
        return tool_failed(tool, exc)


def run_blastn(args: argparse.Namespace, paths: SimulationPaths, work_dir: Path) -> ToolResult:
    tool = "BLASTN"
    makeblastdb_bin = resolve_executable(args.makeblastdb_bin)
    blastn_bin = resolve_executable(args.blastn_bin)

    if makeblastdb_bin is None:
        return tool_missing(tool, args.makeblastdb_bin)
    if blastn_bin is None:
        return tool_missing(tool, args.blastn_bin)
    db_prefix = work_dir / "blastdb" / "combined"
    db_prefix.parent.mkdir(exist_ok=True)
    out = work_dir / "blastn.tsv"
    prediction = work_dir / "blastn.bed"
    make_db = [str(makeblastdb_bin), "-dbtype", "nucl", "-in", str(paths.combined_fasta), "-out", str(db_prefix)]
    blast = [
        str(blastn_bin), "-query", str(paths.combined_fasta), "-db", str(db_prefix),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-num_threads", str(args.threads),
    ]
    try:
        db_profile = profile_command(make_db, work_dir, "blastn_makeblastdb")
        blast_profile = profile_command(blast, work_dir, "blastn", stdout_path=out)
        profile = merge_profiles([db_profile, blast_profile])
        blast_to_bed(out, prediction, args.min_call_length)
        return result_from_prediction(tool, prediction, paths.truth_bed, profile)
    except Exception as exc:
        return tool_failed(tool, exc)


def merge_profiles(profiles: list[Profile]) -> Profile:
    merged_perf: dict[str, float] = {}
    for profile in profiles:
        for key, value in profile.perf.items():
            merged_perf[key] = merged_perf.get(key, 0.0) + value
    return Profile(
        wall_seconds=sum(profile.wall_seconds for profile in profiles),
        max_rss_kb=max((profile.max_rss_kb for profile in profiles), default=0),
        perf=merged_perf,
        command=" && ".join(profile.command for profile in profiles),
    )


def paf_to_bed(paf: Path, bed: Path, min_len: int) -> None:
    intervals: list[tuple[str, int, int]] = []
    if paf.exists():
        with paf.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 9:
                    continue
                query, q_start, q_end = fields[0], int(fields[2]), int(fields[3])
                target, t_start, t_end = fields[5], int(fields[7]), int(fields[8])
                if is_self_diagonal(query, q_start, q_end, target, t_start, t_end):
                    continue
                add_if_long(intervals, query, q_start, q_end, min_len)
                add_if_long(intervals, target, t_start, t_end, min_len)
    write_bed(intervals, bed)


def coords_to_bed(coords: Path, bed: Path, min_len: int) -> None:
    intervals: list[tuple[str, int, int]] = []
    if coords.exists():
        with coords.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 6:
                    continue
                try:
                    ref_start, ref_end, qry_start, qry_end = [int(value) for value in fields[:4]]
                except ValueError:
                    continue
                ref_name, qry_name = fields[-2], fields[-1]
                ref_start, ref_end = sorted((ref_start - 1, ref_end))
                qry_start, qry_end = sorted((qry_start - 1, qry_end))
                if is_self_diagonal(ref_name, ref_start, ref_end, qry_name, qry_start, qry_end):
                    continue
                add_if_long(intervals, ref_name, ref_start, ref_end, min_len)
                add_if_long(intervals, qry_name, qry_start, qry_end, min_len)
    write_bed(intervals, bed)


def blast_to_bed(blast: Path, bed: Path, min_len: int) -> None:
    intervals: list[tuple[str, int, int]] = []
    if blast.exists():
        with blast.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 10:
                    continue
                query, subject = fields[0], fields[1]
                q_start, q_end = sorted((int(fields[6]) - 1, int(fields[7])))
                s_start, s_end = sorted((int(fields[8]) - 1, int(fields[9])))
                if is_self_diagonal(query, q_start, q_end, subject, s_start, s_end):
                    continue
                add_if_long(intervals, query, q_start, q_end, min_len)
                add_if_long(intervals, subject, s_start, s_end, min_len)
    write_bed(intervals, bed)


def is_self_diagonal(left_name: str, left_start: int, left_end: int, right_name: str, right_start: int, right_end: int) -> bool:
    if left_name != right_name:
        return False
    left_len = left_end - left_start
    right_len = right_end - right_start
    if left_len <= 0 or right_len <= 0:
        return False
    overlap = min(left_end, right_end) - max(left_start, right_start)
    return overlap / left_len >= 0.8 and overlap / right_len >= 0.8


def add_if_long(intervals: list[tuple[str, int, int]], chrom: str, start: int, end: int, min_len: int) -> None:
    if end - start >= min_len:
        intervals.append((chrom, start, end))


def disabled_result(tool: str) -> ToolResult:
    return ToolResult(tool=tool, status="disabled", enabled=False)


def result_to_row(result: ToolResult) -> dict[str, object]:
    row: dict[str, object] = {field_name: "" for field_name in CSV_FIELDS}
    row.update({
        "tool": result.tool,
        "status": result.status,
        "enabled": result.enabled,
        "command": result.profile.command,
        "prediction_bed": str(result.prediction_bed) if result.prediction_bed else "",
        "error": result.error,
        "time_perf_seconds": f"{result.profile.wall_seconds:.6f}" if result.profile.wall_seconds else "0",
        "max_rss_kb": result.profile.max_rss_kb,
    })
    row.update(result.metrics)
    return row


def write_results_csv(path: Path, results: list[ToolResult]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for result in results:
            writer.writerow(result_to_row(result))


def write_tool_survey(path: Path, specs: list[ToolSpec]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tool", "role"])
        writer.writerow(["Truth", "synthetic truth BED baseline; expected perfect BP/Frag scores"])
        for spec in specs:
            writer.writerow([spec.name, spec.description])


def print_report(results: list[ToolResult], csv_path: Path) -> None:
    print("=" * 108)
    print("Simulation benchmark summary")
    print("=" * 108)
    print(f"{'Tool':<16} {'Status':<10} {'Time(s)':>9} {'RSS(MB)':>9} {'BP F1':>9} {'Frag F1':>9}  Error")
    print("-" * 108)
    for result in results:
        metrics = result.metrics
        rss_mb = result.profile.max_rss_kb / 1024 if result.profile.max_rss_kb else 0.0
        print(
            f"{result.tool:<16} {result.status:<10} {result.profile.wall_seconds:>9.2f} {rss_mb:>9.1f} "
            f"{metrics.get('bp_f1', 0.0) * 100:>8.2f}% {metrics.get('frag_f1', 0.0) * 100:>8.2f}%  {result.error}"
        )
    print("=" * 108)
    print(f"[INFO] CSV results: {csv_path}")


def tool_specs() -> list[ToolSpec]:
    return [
        ToolSpec("SegTrace", "no_segtrace", "segtrace", "-", run_segtrace),
        ToolSpec("minimap2", "no_minimap2", "minimap2", "-", run_minimap2),
        ToolSpec("MUMmer/nucmer", "no_nucmer", "nucmer", "-", run_nucmer),
        ToolSpec("BLASTN", "no_blastn", "blastn", "-", run_blastn),
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Simulate duplicated genome fragments and benchmark multiple tools.")
    parser.add_argument("--out-dir", default="sim_benchmark", help="Output directory")
    parser.add_argument("--out-csv", default="benchmark_results.csv", help="Results CSV name under --out-dir")
    parser.add_argument("--force", action="store_true", help="Replace an existing output directory")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--species", type=int, default=3, help="Number of species/genomes")
    parser.add_argument("--chromosomes", type=int, default=5, help="Chromosomes per species")
    parser.add_argument("--chrom-length", type=int, default=100_000_000, help="Length of each chromosome")
    parser.add_argument("--fragments", type=int, default=500, help="Number of source fragments")
    parser.add_argument("--min-fragment-length", type=int, default=1_000, help="Minimum source fragment length")
    parser.add_argument("--max-fragment-length", type=int, default=50_000, help="Maximum source fragment length")
    parser.add_argument("--min-copies", type=int, default=2, help="Minimum copies per source fragment")
    parser.add_argument("--max-copies", type=int, default=10, help="Maximum copies per source fragment")
    parser.add_argument("--max-snp-rate", type=float, default=0.10, help="Maximum SNP rate per copy")
    parser.add_argument("--max-indel-rate", type=float, default=0.01, help="Maximum INDEL rate per copy")
    parser.add_argument("--filler-chunk", type=int, default=1_000_000, help="Random filler chunk size")
    parser.add_argument("--skip-tools", action="store_true", help="Only generate FASTA and truth BED")
    parser.add_argument("--min-call-length", type=int, default=1_000, help="Minimum converted interval length for mapper outputs")
    parser.add_argument("--max-mappings", type=int, default=1000, help="Maximum mappings/chains retained per segment where supported")

    parser.add_argument("--threads", "-p", type=int, default=8, help="Tool threads where supported")
    parser.add_argument("--segtrace-bin", default="./segtrace", help="SegTrace executable")
    parser.add_argument("--kmer", "-k", type=int, default=19, help="SegTrace k-mer size")
    parser.add_argument("--scale", "-s", type=int, default=16, help="SegTrace scale factor")
    parser.add_argument("--window-size", "-w", type=int, default=1024, help="SegTrace window size")
    parser.add_argument("--step-size", "-t", type=int, default=0, help="SegTrace step size")
    parser.add_argument("--min-report-copies", "-c", type=int, default=1, help="SegTrace -c value")
    parser.add_argument("--minimap2-bin", default="minimap2", help="minimap2 executable")
    parser.add_argument("--nucmer-bin", default="nucmer", help="nucmer executable")
    parser.add_argument("--show-coords-bin", default="show-coords", help="show-coords executable")
    parser.add_argument("--blastn-bin", default="blastn", help="blastn executable")
    parser.add_argument("--makeblastdb-bin", default="makeblastdb", help="makeblastdb executable")

    parser.add_argument("--no-truth", action="store_true", help="Disable the synthetic truth baseline row")
    for spec in tool_specs():
        parser.add_argument(f"--no-{spec.cli_name}", dest=spec.no_flag, action="store_true", help=f"Disable {spec.name}")
    args = parser.parse_args()

    if args.min_fragment_length > args.max_fragment_length:
        parser.error("--min-fragment-length must be <= --max-fragment-length")
    if args.min_copies > args.max_copies:
        parser.error("--min-copies must be <= --max-copies")
    if args.chrom_length < args.max_fragment_length:
        parser.error("--chrom-length must be >= --max-fragment-length")
    if args.species < 1 or args.chromosomes < 1 or args.fragments < 1:
        parser.error("--species, --chromosomes, and --fragments must be positive")
    return args


def main() -> int:
    args = parse_args()
    paths, placements = generate_simulation(args)
    out_dir = Path(args.out_dir)
    specs = tool_specs()
    write_tool_survey(out_dir / "benchmark_tools.tsv", specs)
    print(f"[INFO] Wrote {len(paths.fasta_paths)} FASTA files, combined FASTA, and {len(placements)} truth intervals to {out_dir}")

    if args.skip_tools:
        print("[INFO] Skipping benchmark tools")
        return 0

    results: list[ToolResult] = []
    if args.no_truth:
        results.append(disabled_result("Truth"))
    else:
        results.append(run_truth(args, paths, out_dir))

    for spec in specs:
        if getattr(args, spec.no_flag):
            results.append(disabled_result(spec.name))
            continue
        print(f"[INFO] Running {spec.name}...")
        results.append(spec.runner(args, paths, out_dir))

    csv_path = out_dir / args.out_csv
    write_results_csv(csv_path, results)
    print_report(results, csv_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
