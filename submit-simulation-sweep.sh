#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-simulation-sweep}"
THREADS="${THREADS:-16}"
MEM="${MEM:-480gb}"
WALLTIME="${WALLTIME:-96:00:00}"
QUEUE_NODE="${QUEUE_NODE:-node02}"
ENV_SETUP="${ENV_SETUP:-}"

export SNP_VALUES="${SNP_VALUES:-1 2 5 10}"
export INDEL_VALUES="${INDEL_VALUES:-0.1 0.2 0.5 1}"
export GENOME_SIZES="${GENOME_SIZES:-1Mb 10Mb 100Mb 1Gb}"

export OUT_DIR="${OUT_DIR:-simulation_sweep}"
export OUT_CSV="${OUT_CSV:-simulation_sweep.csv}"
export BASE_SEED="${BASE_SEED:-42}"
export MIN_COPIES="${MIN_COPIES:-2}"
export MAX_COPIES="${MAX_COPIES:-10}"
export MIN_FRAGMENT_LENGTH="${MIN_FRAGMENT_LENGTH:-}"
export MAX_FRAGMENT_LENGTH="${MAX_FRAGMENT_LENGTH:-}"
export KEEP_CASES="${KEEP_CASES:-0}"
export SEGTRACE_BIN="${SEGTRACE_BIN:-./segtrace}"
export MINIMAP2_BIN="${MINIMAP2_BIN:-minimap2}"
export MIN_CALL_LENGTH="${MIN_CALL_LENGTH:-1}"
export MAX_MAPPINGS="${MAX_MAPPINGS:-1000}"
export WORKDIR THREADS ENV_SETUP

mkdir -p "$LOGDIR"

qsub -N "$JOB_NAME" \
    -l select=1:ncpus=${THREADS} \
  -l walltime=${WALLTIME} \
  -v WORKDIR="${WORKDIR}",THREADS="${THREADS}",OUT_DIR="${OUT_DIR}",OUT_CSV="${OUT_CSV}",ENV_SETUP="${ENV_SETUP}" \
  -V \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
export PATH="$HOME/.local/bin:$PATH"
${ENV_SETUP:+$ENV_SETUP}

if command -v micromamba >/dev/null 2>&1; then
	RUNNER=(micromamba run -n segtrace-compare uv run)
else
	RUNNER=(uv run)
fi

make clean && make

echo "=== Running simulation sweep in ${OUT_DIR} ==="
"${RUNNER[@]}" --with numpy python - <<'PY'
import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path

import sim_benchmark as sb


def values(name: str) -> list[str]:
    parsed = os.environ.get(name, "").split()
    if not parsed:
        raise ValueError(f"{name} must contain at least one value")
    return parsed


def parse_percent(token: str, name: str) -> tuple[str, float]:
    value = float(token)
    if not 0.0 <= value <= 100.0:
        raise ValueError(f"{name} values must be between 0 and 100: {token}")
    return token.replace(".", "p"), value / 100.0


def parse_size(token: str) -> int:
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([kmgt]?i?b?)", token.lower())
    if match is None:
        raise ValueError(f"invalid genome size: {token}")
    units = {
        "": 1,
        "b": 1,
        "k": 1_000,
        "kb": 1_000,
        "ki": 1_024,
        "kib": 1_024,
        "m": 1_000_000,
        "mb": 1_000_000,
        "mi": 1_048_576,
        "mib": 1_048_576,
        "g": 1_000_000_000,
        "gb": 1_000_000_000,
        "gi": 1_073_741_824,
        "gib": 1_073_741_824,
        "t": 1_000_000_000_000,
        "tb": 1_000_000_000_000,
        "ti": 1_099_511_627_776,
        "tib": 1_099_511_627_776,
    }
    number, suffix = match.groups()
    if suffix not in units:
        raise ValueError(f"invalid genome size suffix: {token}")
    return int(float(number) * units[suffix])


def safe_label(token: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "p", token)


def failure_row(tool: str, error: Exception) -> dict[str, object]:
    return sb.result_to_row(sb.tool_failed(tool, error))


def run_segtrace_with_defaults(segtrace_bin: str, paths: sb.SimulationPaths, work_dir: Path, threads: int) -> sb.ToolResult:
    tool = "SegTrace"
    executable = sb.resolve_executable(segtrace_bin)
    if executable is None:
        return sb.tool_missing(tool, segtrace_bin)
    prefix = work_dir / "segtrace"
    prediction = prefix.with_suffix(".seg.bed")
    command = [
        str(executable), "-o", str(prefix), "-p", str(threads),
        *[str(path) for path in paths.fasta_paths],
    ]
    try:
        profile = sb.profile_command(command, work_dir, "segtrace")
        if not prediction.exists():
            prediction = sb.empty_bed(work_dir / "segtrace.empty.bed")
        return sb.result_from_prediction(tool, prediction, paths.truth_bed, profile)
    except Exception as exc:
        return sb.tool_failed(tool, exc)


out_dir = Path(os.environ["OUT_DIR"])
out_dir.mkdir(parents=True, exist_ok=True)
csv_path = out_dir / os.environ["OUT_CSV"]
case_root = out_dir / "cases"
case_root.mkdir(exist_ok=True)

snp_tokens = values("SNP_VALUES")
indel_tokens = values("INDEL_VALUES")
genome_tokens = values("GENOME_SIZES")
snp_values = [parse_percent(token, "SNP_VALUES") for token in snp_tokens]
indel_values = [parse_percent(token, "INDEL_VALUES") for token in indel_tokens]
genome_values = [(token, parse_size(token)) for token in genome_tokens]

min_copies = int(os.environ["MIN_COPIES"])
max_copies = int(os.environ["MAX_COPIES"])
if min_copies < 1 or min_copies > max_copies:
    raise ValueError("MIN_COPIES and MAX_COPIES must be positive and ordered")

base_seed = int(os.environ["BASE_SEED"])
threads = int(os.environ["THREADS"])
keep_cases = os.environ.get("KEEP_CASES", "0") == "1"
explicit_min_length = os.environ.get("MIN_FRAGMENT_LENGTH", "")
explicit_max_length = os.environ.get("MAX_FRAGMENT_LENGTH", "")

metadata_fields = [
    "snp_percent", "indel_percent", "genome_size", "genome_bp", "segments", "seed",
]
fieldnames = metadata_fields + sb.CSV_FIELDS
total = len(snp_values) * len(indel_values) * len(genome_values)
done = 0

with csv_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()

    for snp_label, snp_rate in snp_values:
        for indel_label, indel_rate in indel_values:
            for genome_label, genome_bp in genome_values:
                done += 1
                segments = max(1, genome_bp // 1_000_000)
                label = (
                    f"snp{snp_label}_indel{indel_label}_size{safe_label(genome_label)}"
                    f"_seg{segments}"
                )
                case_dir = case_root / label
                max_fragment_length = genome_bp // max(segments * max_copies * 4, 1)
                max_fragment_length = max(1, max_fragment_length)
                min_fragment_length = max(1, max_fragment_length // 2)
                if explicit_min_length:
                    min_fragment_length = int(explicit_min_length)
                if explicit_max_length:
                    max_fragment_length = int(explicit_max_length)
                if min_fragment_length > max_fragment_length:
                    raise ValueError(
                        f"fragment length range is invalid for {label}: "
                        f"{min_fragment_length}>{max_fragment_length}"
                    )
                seed = base_seed + done - 1
                metadata = {
                    "snp_percent": snp_label,
                    "indel_percent": indel_label,
                    "genome_size": genome_label,
                    "genome_bp": genome_bp,
                    "segments": segments,
                    "seed": seed,
                }
                print(f"[{done}/{total}] generating {label}", flush=True)

                try:
                    simulation_args = argparse.Namespace(
                        out_dir=str(case_dir),
                        force=True,
                        seed=seed,
                        species=1,
                        chromosomes=1,
                        chrom_length=genome_bp,
                        fragments=segments,
                        min_fragment_length=min_fragment_length,
                        max_fragment_length=max_fragment_length,
                        min_copies=min_copies,
                        max_copies=max_copies,
                        max_snp_rate=snp_rate,
                        max_indel_rate=indel_rate,
                        filler_chunk=1_000_000,
                    )
                    paths, placements = sb.generate_simulation(simulation_args)
                    tool_args = argparse.Namespace(
                        threads=threads,
                        minimap2_bin=os.environ["MINIMAP2_BIN"],
                        max_mappings=int(os.environ["MAX_MAPPINGS"]),
                        min_call_length=int(os.environ["MIN_CALL_LENGTH"]),
                    )
                    minimap_dir = case_dir / "minimap2"
                    segtrace_dir = case_dir / "segtrace"
                    minimap_dir.mkdir()
                    segtrace_dir.mkdir()
                    minimap_result = sb.run_minimap2(tool_args, paths, minimap_dir)
                    segtrace_result = run_segtrace_with_defaults(
                        os.environ["SEGTRACE_BIN"], paths, segtrace_dir, threads,
                    )
                    results = [minimap_result, segtrace_result]
                    print(
                        f"[{done}/{total}] {label} placements={len(placements)} "
                        f"minimap2={minimap_result.status} segtrace={segtrace_result.status}",
                        flush=True,
                    )
                except Exception as exc:
                    results = [failure_row("minimap2", exc), failure_row("SegTrace", exc)]
                    print(f"[{done}/{total}] {label} failed: {exc}", file=sys.stderr, flush=True)

                for result in results:
                    row = dict(metadata)
                    row.update(result if isinstance(result, dict) else sb.result_to_row(result))
                    if not keep_cases:
                        row["prediction_bed"] = ""
                    writer.writerow(row)
                handle.flush()

                if not keep_cases:
                    shutil.rmtree(case_dir, ignore_errors=True)

print(f"[INFO] simulation sweep CSV: {csv_path}")
PY

echo "=== Simulation sweep complete: ${OUT_DIR}/${OUT_CSV} ==="
PBS