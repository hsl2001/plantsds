#!/usr/bin/env bash
# Compare SegTrace segments with non-diagonal minimap2 self-alignments.
#
# Usage:
#   ./compare_t2t_nip.sh [minimap2.paf] [segtrace.seg.bed] [output-prefix]
#
# The default inputs are t2t_nip.minimap2.paf and t2t_nip.seg.bed.  Results:
#   <prefix>.venn.txt
#   <prefix>.segtrace-specific.tsv
#   <prefix>.minimap2-specific.tsv
#
# Matching can be adjusted without editing this script:
#   MIN_OVERLAP_FRAC=0.5 MIN_OVERLAP_BP=100 MIN_ALIGN_LEN=1000 ./compare_t2t_nip.sh
set -euo pipefail

PAF=${1:-t2t_nip.minimap2.paf}
BED=${2:-t2t_nip.seg.bed}
PREFIX=${3:-t2t_nip.compare}
MIN_OVERLAP_FRAC=${MIN_OVERLAP_FRAC:-0.5}
MIN_OVERLAP_BP=${MIN_OVERLAP_BP:-100}
MIN_ALIGN_LEN=${MIN_ALIGN_LEN:-1000}

if [[ ! -r "$PAF" ]]; then
  printf 'error: minimap2 PAF not found or not readable: %s\n' "$PAF" >&2
  exit 1
fi
if [[ ! -r "$BED" ]]; then
  printf 'error: SegTrace BED not found or not readable: %s\n' "$BED" >&2
  exit 1
fi

python3 - "$PAF" "$BED" "$PREFIX" "$MIN_OVERLAP_FRAC" "$MIN_OVERLAP_BP" "$MIN_ALIGN_LEN" <<'PY'
import bisect
import sys
from collections import defaultdict


paf_path, bed_path, prefix, min_frac_s, min_bp_s, min_align_len_s = sys.argv[1:]
min_frac = float(min_frac_s)
min_bp = int(min_bp_s)
min_align_len = int(min_align_len_s)

if not 0 < min_frac <= 1:
    raise SystemExit("error: MIN_OVERLAP_FRAC must be in (0, 1]")
if min_bp < 1 or min_align_len < 1:
    raise SystemExit("error: MIN_OVERLAP_BP and MIN_ALIGN_LEN must be positive")


def read_bed(path):
    intervals = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 3:
                print(f"warning: skipping malformed BED line {line_number}", file=sys.stderr)
                continue
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            except ValueError:
                print(f"warning: skipping malformed BED line {line_number}", file=sys.stderr)
                continue
            if end <= start:
                continue
            cluster = fields[3] if len(fields) >= 4 else "."
            intervals.append((chrom, start, end, cluster, line_number))
    return intervals


def read_paf(path):
    """Return query and target intervals, excluding exact self-diagonal hits."""
    intervals = []
    skipped_diagonal = 0
    paf_chroms = set()
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                print(f"warning: skipping malformed PAF line {line_number}", file=sys.stderr)
                continue
            try:
                query, qstart, qend = fields[0], int(fields[2]), int(fields[3])
                target, tstart, tend = fields[5], int(fields[7]), int(fields[8])
                block_length = int(fields[10])
            except ValueError:
                print(f"warning: skipping malformed PAF line {line_number}", file=sys.stderr)
                continue
            paf_chroms.update((query, target))
            if block_length < min_align_len or qend <= qstart or tend <= tstart:
                continue

            # minimap2 reports the query-to-itself diagonal for self alignment;
            # retaining it would make the whole genome appear to be a repeat.
            if query == target and qstart == tstart and qend == tend:
                skipped_diagonal += 1
                continue

            intervals.append((query, qstart, qend, f"{line_number}:query"))
            intervals.append((target, tstart, tend, f"{line_number}:target"))
    return intervals, skipped_diagonal, paf_chroms


def normalize_bed_intervals(intervals, paf_chroms):
    """Match prefixed SegTrace names to the corresponding PAF contig names."""
    normalized = []
    changed = 0
    for chrom, start, end, cluster, line_number in intervals:
        match_chrom = chrom
        if chrom not in paf_chroms:
            parts = chrom.split("-")
            for index in range(1, len(parts)):
                candidate = "-".join(parts[index:])
                if candidate in paf_chroms:
                    match_chrom = candidate
                    changed += 1
                    break
        # Keep the original chromosome for the result TSV, but use match_chrom
        # for interval comparisons.
        normalized.append((match_chrom, start, end, cluster, line_number, chrom))
    return normalized, changed


def merge_intervals(intervals):
    grouped = defaultdict(list)
    for chrom, start, end, *_ in intervals:
        grouped[chrom].append((start, end))

    merged = {}
    for chrom, values in grouped.items():
        values.sort()
        result = []
        for start, end in values:
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], end))
            else:
                result.append((start, end))
        merged[chrom] = result
    return merged


def interval_matches(interval, merged_other):
    chrom, start, end = interval[:3]
    candidates = merged_other.get(chrom, ())
    if not candidates:
        return False

    starts = [candidate[0] for candidate in candidates]
    index = max(0, bisect.bisect_right(starts, start) - 1)
    interval_length = end - start
    while index < len(candidates):
        other_start, other_end = candidates[index]
        if other_start >= end:
            break
        overlap = min(end, other_end) - max(start, other_start)
        if overlap >= min_bp and overlap >= min_frac * min(interval_length, other_end - other_start):
            return True
        index += 1
    return False


def write_tsv(path, header, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


segtrace_raw = read_bed(bed_path)
minimap2, skipped_diagonal, paf_chroms = read_paf(paf_path)
segtrace, normalized_bed_count = normalize_bed_intervals(segtrace_raw, paf_chroms)
segtrace_coverage = merge_intervals(segtrace)
minimap2_coverage = merge_intervals(minimap2)

segtrace_matches = [interval_matches(interval, minimap2_coverage) for interval in segtrace]
minimap2_matches = [interval_matches(interval, segtrace_coverage) for interval in minimap2]

segtrace_specific = [interval for interval, matched in zip(segtrace, segtrace_matches) if not matched]
minimap2_specific = [interval for interval, matched in zip(minimap2, minimap2_matches) if not matched]

# A PAF row contributes query and target intervals. Exact duplicate coordinates
# are collapsed so secondary records do not inflate the minimap2-only list.
unique_minimap2 = []
seen = set()
for interval, matched in zip(minimap2, minimap2_matches):
    key = interval[:3]
    if not matched and key not in seen:
        seen.add(key)
        unique_minimap2.append(interval)

segtrace_only_path = f"{prefix}.segtrace-specific.tsv"
minimap2_only_path = f"{prefix}.minimap2-specific.tsv"
venn_path = f"{prefix}.venn.txt"

write_tsv(
    segtrace_only_path,
    ("chrom", "start", "end", "cluster_id", "bed_line"),
    ((original_chrom, start, end, cluster, line_number)
     for _, start, end, cluster, line_number, original_chrom in segtrace_specific),
)
write_tsv(
    minimap2_only_path,
    ("chrom", "start", "end", "paf_record"),
    ((chrom, start, end, record_id) for chrom, start, end, record_id in unique_minimap2),
)

segtrace_both = sum(segtrace_matches)
minimap2_both = sum(minimap2_matches)
segtrace_only_count = len(segtrace) - segtrace_both
minimap2_unique_count = len({interval[:3] for interval in minimap2})
minimap2_only_count = minimap2_unique_count - len({interval[:3] for interval, matched in zip(minimap2, minimap2_matches) if matched})

venn = f"""SegTrace vs minimap2 interval comparison
=========================================
Matching rule: overlap >= {min_frac:g} of the shorter interval and >= {min_bp} bp
minimap2 filter: alignment block length >= {min_align_len} bp; exact same-sequence diagonal hits excluded

                         .-----------------------.
                    .---'                         '---.
                 .-'       SegTrace only: {segtrace_only_count:>8}          '-.
                /                                             \\
               /       Both: {segtrace_both:>8} SegTrace intervals             \\
               \\       Both: {minimap2_both:>8} minimap2 intervals             /
                \\                                             /
                 '-.      minimap2 only: {minimap2_only_count:>8}       .-'
                    '---.                         .---'
                        '-----------------------'

Input intervals:
  SegTrace: {len(segtrace)}
    BED chromosome names normalized to PAF names: {normalized_bed_count}
    minimap2 query/target intervals: {len(minimap2)} ({minimap2_unique_count} unique coordinates)
  exact minimap2 self-diagonal records skipped: {skipped_diagonal}

Lists:
  SegTrace-specific: {segtrace_only_path}
  minimap2-specific: {minimap2_only_path}
"""

with open(venn_path, "w") as handle:
    handle.write(venn)
print(venn, end="")
PY