#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    printf 'Usage: %s pan-arabidopsis.seg.bed CLUSTER_ID\n' "$0" >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
bed_path=$1
cluster_id=$2

if [[ ! -f "$bed_path" && -f "$script_dir/cons/$bed_path" ]]; then
    bed_path="$script_dir/cons/$bed_path"
fi
if [[ ! -f "$bed_path" ]]; then
    printf '[ERROR] BED file not found: %s\n' "$1" >&2
    exit 1
fi

fasta_dir=${FASTA_DIR:-"$script_dir/cons/pan_arabidopsis_fasta"}
bed_stem=${bed_path%.seg.bed}
if [[ "$bed_stem" == "$bed_path" ]]; then
    bed_stem=${bed_path%.bed}
fi
output_path=${OUTPUT_FASTA:-"${bed_stem}.cluster-${cluster_id}.fasta"}

python3 - "$bed_path" "$cluster_id" "$fasta_dir" "$output_path" <<'PY'
import sys
from collections import defaultdict
from pathlib import Path

bed_path = Path(sys.argv[1])
cluster_id = sys.argv[2]
fasta_dir = Path(sys.argv[3])
output_path = Path(sys.argv[4])

segments = defaultdict(list)
with bed_path.open(encoding="utf-8") as bed_file:
    for line_number, line in enumerate(bed_file, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"{bed_path}:{line_number}: expected at least 4 BED columns")
        if fields[3] != cluster_id:
            continue

        ecotype, separator, seqid = fields[0].rpartition("-")
        if not separator:
            raise ValueError(f"{bed_path}:{line_number}: invalid genome-sequence name")
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError as error:
            raise ValueError(f"{bed_path}:{line_number}: invalid BED coordinates") from error
        if start < 0 or end <= start:
            raise ValueError(f"{bed_path}:{line_number}: invalid BED interval")
        segments[ecotype].append((seqid, start, end, line_number))

if not segments:
    raise ValueError(f"No BED records found for cluster {cluster_id!r} in {bed_path}")


def load_records(fasta_path, wanted_names):
    records = {}
    current_name = None
    chunks = []

    def finish_record():
        if current_name in wanted_names:
            if current_name in records:
                raise ValueError(f"{fasta_path}: duplicate FASTA record {current_name!r}")
            records[current_name] = "".join(chunks).upper()

    with fasta_path.open(encoding="utf-8") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                if current_name is not None:
                    finish_record()
                current_name = line[1:].split()[0]
                chunks = []
            elif current_name in wanted_names:
                chunks.append(line.strip())
        if current_name is not None:
            finish_record()
    return records


def find_fasta(ecotype):
    for suffix in (".fasta", ".fa", ".fna"):
        candidate = fasta_dir / f"{ecotype}{suffix}"
        if candidate.is_file():
            return candidate
    return None


output_path.parent.mkdir(parents=True, exist_ok=True)
written = 0
skipped = 0
with output_path.open("w", encoding="utf-8") as output_file:
    for ecotype, ecotype_segments in segments.items():
        fasta_path = find_fasta(ecotype)
        if fasta_path is None:
            skipped += len(ecotype_segments)
            print(
                f"[WARN] FASTA not found for {ecotype}; skipped "
                f"{len(ecotype_segments)} segment(s)",
                file=sys.stderr,
            )
            continue

        records = load_records(fasta_path, {segment[0] for segment in ecotype_segments})
        for segment_number, (seqid, start, end, line_number) in enumerate(
            ecotype_segments, start=1
        ):
            if seqid not in records:
                raise ValueError(
                    f"{fasta_path}: FASTA record {seqid!r} not found "
                    f"(referenced at {bed_path}:{line_number})"
                )
            sequence = records[seqid]
            if end > len(sequence):
                raise ValueError(
                    f"{bed_path}:{line_number}: interval {start}-{end} exceeds "
                    f"{fasta_path}:{seqid} length {len(sequence)}"
                )

            extracted = sequence[start:end]
            output_file.write(
                f">cluster_{cluster_id}|ecotype={ecotype}|seqid={seqid}|"
                f"start={start + 1}|end={end}|segment={segment_number}\n"
            )
            for offset in range(0, len(extracted), 80):
                output_file.write(extracted[offset : offset + 80] + "\n")
            written += 1

print(
    f"Wrote {output_path} ({written} segment sequences; {skipped} skipped)",
    file=sys.stderr,
)
PY

last_env=${LAST_ENV:-benchmark-segtrace}
last_threads=${LAST_THREADS:-1}
output_prefix=${output_path%.*}
maf_output=${LAST_MAF:-"${output_prefix}.last.maf"}
dotplot_output=${LAST_DOTPLOT:-"${output_prefix}.last-dotplot.png"}
dotplot_width=${LAST_DOTPLOT_WIDTH:-4000}
dotplot_height=${LAST_DOTPLOT_HEIGHT:-900}
last_temp_dir=$(mktemp -d)
trap 'rm -rf "$last_temp_dir"' EXIT
longest_fasta="$last_temp_dir/longest.fasta"
others_fasta="$last_temp_dir/others.fasta"
lastdb_prefix="$last_temp_dir/others"

python3 - "$output_path" "$longest_fasta" "$others_fasta" <<'PY'
import sys
import re
from pathlib import Path

input_path, longest_path, others_path = map(Path, sys.argv[1:])
records = []
name = None
chunks = []

with input_path.open(encoding="utf-8") as input_file:
    for line in input_file:
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))

if len(records) < 2:
    raise ValueError("at least two extracted sequences are required for a dotplot")


def natural_key(value):
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", value)
    )


def compact_name(header, length, longest=False, self_copy=False):
    values = {}
    for field in header.split("|"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
    if longest and self_copy:
        prefix = "longest_self"
    elif longest:
        prefix = "longest"
    else:
        prefix = values.get("ecotype", "sequence")
    return (
        f"{prefix}|{values.get('seqid', 'sequence')}|"
        f"segment={values.get('segment', '1')}|length={length}"
    )


records.sort(key=lambda record: (-len(record[1]), natural_key(record[0])))


def print_label(axis, rank, header, length, longest=False, self_copy=False):
    values = {}
    for field in header.split("|"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
    label = compact_name(header, length, longest, self_copy)
    print(
        f"{axis}\t{rank}\t{label}\t{values.get('ecotype', '')}\t"
        f"{values.get('seqid', '')}"
    )


def write_fasta(path, selected, self_copy=False):
    with path.open("w", encoding="utf-8") as output_file:
        for header, sequence, longest in selected:
            output_file.write(
                f">{compact_name(header, len(sequence), longest, self_copy)}\n"
            )
            for offset in range(0, len(sequence), 80):
                output_file.write(sequence[offset : offset + 80] + "\n")


write_fasta(
    longest_path,
    [(records[0][0], records[0][1], True)],
)
write_fasta(
    others_path,
    [(records[0][0], records[0][1], True)]
    + [(header, sequence, False) for header, sequence in records[1:]],
    self_copy=True,
)
print("axis\trank\tlabel\tecotype\tchromosome")
print_label("x", 1, records[0][0], len(records[0][1]), True, True)
for rank, (header, sequence) in enumerate(records[1:], start=2):
    print_label("x", rank, header, len(sequence))
print_label("y", 1, records[0][0], len(records[0][1]), True)
print(
    f"Longest sequence: {records[0][0]} ({len(records[0][1])} bp); "
    f"x-axis sequences including self: {len(records)}",
    file=sys.stderr,
)
PY

mkdir -p "$(dirname -- "$maf_output")" "$(dirname -- "$dotplot_output")"
micromamba run -n "$last_env" lastdb -P "$last_threads" "$lastdb_prefix" "$others_fasta"
micromamba run -n "$last_env" lastal -P "$last_threads" "$lastdb_prefix" \
    "$longest_fasta" > "$maf_output"
micromamba run -n "$last_env" last-dotplot \
    --sort1=2 \
    --sort2=0 \
    --labels1=0 \
    --labels2=0 \
    --rot1=h \
    --rot2=v \
    --max-gap1=1000000000,1000000000 \
    --max-gap2=1000000000,1000000000 \
    --maxseqs=10000 \
    --width="$dotplot_width" \
    --height="$dotplot_height" \
    "$maf_output" "$dotplot_output"
printf 'Wrote %s and %s using %s\n' "$maf_output" "$dotplot_output" "$last_env" >&2