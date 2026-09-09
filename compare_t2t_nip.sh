#!/usr/bin/env bash
# Run SegTrace and an all-mapping minimap2 self-alignment, then compare them.
#
# Usage:
#   ./compare_t2t_nip.sh [t2t_nip.fasta] [output-prefix]
#
# The default input is t2t_nip.fasta.  With output prefix t2t_nip, the tool
# outputs are t2t_nip.seg.bed and t2t_nip.minimap2.paf.  The merged intervals
# and comparison results are:
#   <prefix>.seg.merged.bed
#   <prefix>.minimap2.merged.bed
#   <prefix>.compare.venn.txt
#   <prefix>.compare.segtrace-specific.tsv
#   <prefix>.compare.minimap2-specific.tsv
#
# Tool and matching parameters can be adjusted without editing this script:
#   THREADS=36 MAX_MAPPINGS=1000000 ./compare_t2t_nip.sh t2t_nip.fasta t2t_nip
set -euo pipefail

FASTA=${1:-t2t_nip.fasta}
PREFIX=${2:-t2t_nip}
PAF="${PREFIX}.minimap2.paf"
BED="${PREFIX}.seg.bed"
SEG_NORMALIZED="${PREFIX}.seg.normalized.bed"
SEG_MERGED="${PREFIX}.seg.merged.bed"
MINIMAP2_INTERVALS="${PREFIX}.minimap2.intervals.bed"
MINIMAP2_MERGED="${PREFIX}.minimap2.merged.bed"
COMPARE_PREFIX="${PREFIX}.compare"
SEGTRACE_BIN=${SEGTRACE_BIN:-./segtrace}
MINIMAP2_BIN=${MINIMAP2_BIN:-minimap2}
BEDTOOLS_BIN=${BEDTOOLS_BIN:-bedtools}
THREADS=${THREADS:-36}
MIN_COPIES=${MIN_COPIES:-2}
MINIMAP2_PRESET=${MINIMAP2_PRESET:-asm20}
MAX_MAPPINGS=${MAX_MAPPINGS:-1000000}
MIN_OVERLAP_FRAC=${MIN_OVERLAP_FRAC:-0.5}
MIN_OVERLAP_BP=${MIN_OVERLAP_BP:-100}
MIN_ALIGN_LEN=${MIN_ALIGN_LEN:-1000}

if [[ ! -r "$FASTA" ]]; then
    printf 'error: FASTA not found or not readable: %s\n' "$FASTA" >&2
    exit 1
fi

if [[ ! -x "$SEGTRACE_BIN" ]]; then
    printf 'error: SegTrace executable not found or not executable: %s\n' "$SEGTRACE_BIN" >&2
    exit 1
fi
if ! command -v "$MINIMAP2_BIN" >/dev/null 2>&1; then
    printf 'error: minimap2 executable not found in PATH: %s\n' "$MINIMAP2_BIN" >&2
    exit 1
fi
if ! command -v "$BEDTOOLS_BIN" >/dev/null 2>&1; then
    printf 'error: bedtools executable not found in PATH: %s\n' "$BEDTOOLS_BIN" >&2
    exit 1
fi

printf '[1/4] Running SegTrace (-c %s)...\n' "$MIN_COPIES" >&2
"$SEGTRACE_BIN" -c "$MIN_COPIES" -p "$THREADS" -o "$PREFIX" "$FASTA"
if [[ ! -s "$BED" ]]; then
    printf 'error: SegTrace did not produce a non-empty BED: %s\n' "$BED" >&2
    exit 1
fi

# -P retains all chains, while -p 0 and a large -N keep low-scoring and
# numerous secondary mappings.
printf '[2/4] Running minimap2 (%s; all mappings)...\n' "$MINIMAP2_PRESET" >&2
"$MINIMAP2_BIN" -x "$MINIMAP2_PRESET" -P -p 0 \
    -N "$MAX_MAPPINGS" -t "$THREADS" "$FASTA" "$FASTA" > "$PAF"
if [[ ! -s "$PAF" ]]; then
    printf 'error: minimap2 did not produce a non-empty PAF: %s\n' "$PAF" >&2
    exit 1
fi

printf '[3/4] Normalizing and merging intervals with bedtools...\n' >&2
# SegTrace prefixes contig names with the FASTA basename (for example,
# t2t_nip-NC_089035.1), while minimap2 uses the FASTA header (NC_089035.1).
# Normalize the BED names against names observed in the PAF before merging.
awk '
    NR == FNR { paf_chrom[$1] = 1; paf_chrom[$6] = 1; next }
    /^#/ { next }
    NF >= 3 {
        chrom = $1
        if (!(chrom in paf_chrom)) {
            n = split(chrom, parts, "-")
            candidate = chrom
            for (i = 2; i <= n; i++) {
                candidate = parts[i]
                for (j = i + 1; j <= n; j++) candidate = candidate "-" parts[j]
                if (candidate in paf_chrom) { chrom = candidate; break }
            }
        }
        printf "%s\t%s\t%s\t%s\n", chrom, $2, $3, (NF >= 4 ? $4 : ".")
    }
' "$PAF" "$BED" > "$SEG_NORMALIZED"

LC_ALL=C sort --parallel="$THREADS" -k1,1 -k2,2n -k3,3n "$SEG_NORMALIZED" \
    | "$BEDTOOLS_BIN" merge -i - -c 4 -o distinct > "$SEG_MERGED"

# Keep both query and target intervals from every non-diagonal PAF alignment.
awk -v min_len="$MIN_ALIGN_LEN" '
    BEGIN { OFS = "\t" }
    NF >= 12 && $11 >= min_len && !($1 == $6 && $3 == $8 && $4 == $9) {
        print $1, $3, $4
        print $6, $8, $9
    }
' "$PAF" > "$MINIMAP2_INTERVALS"
LC_ALL=C sort --parallel="$THREADS" -k1,1 -k2,2n -k3,3n "$MINIMAP2_INTERVALS" \
    | "$BEDTOOLS_BIN" merge -i - > "$MINIMAP2_MERGED"

if [[ ! -s "$SEG_MERGED" || ! -s "$MINIMAP2_MERGED" ]]; then
    printf 'error: bedtools merge produced an empty result\n' >&2
    exit 1
fi

printf '[4/4] Comparing merged SegTrace and minimap2 intervals...\n' >&2
python3 - "$SEG_MERGED" "$MINIMAP2_MERGED" "$COMPARE_PREFIX" "$MIN_OVERLAP_FRAC" "$MIN_OVERLAP_BP" "$MIN_ALIGN_LEN" <<'PY'
import bisect
import sys
from collections import defaultdict


seg_merged_path, minimap2_merged_path, prefix, min_frac_s, min_bp_s, min_align_len_s = sys.argv[1:]
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


def coverage_venn(first, second):
    """Return first-only, shared, and second-only covered bases."""
    first_total = sum(
        end - start for intervals in first.values() for start, end in intervals
    )
    second_total = sum(
        end - start for intervals in second.values() for start, end in intervals
    )
    shared = 0
    for chrom in set(first) | set(second):
        left = first.get(chrom, ())
        right = second.get(chrom, ())
        i = j = 0
        while i < len(left) and j < len(right):
            left_start, left_end = left[i]
            right_start, right_end = right[j]
            overlap_start = max(left_start, right_start)
            overlap_end = min(left_end, right_end)
            if overlap_start < overlap_end:
                shared += overlap_end - overlap_start
            if left_end <= right_end:
                i += 1
            else:
                j += 1
    return first_total - shared, shared, second_total - shared


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


segtrace = read_bed(seg_merged_path)
minimap2 = read_bed(minimap2_merged_path)
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
    ((chrom, start, end, cluster, line_number)
     for chrom, start, end, cluster, line_number in segtrace_specific),
)
write_tsv(
    minimap2_only_path,
    ("chrom", "start", "end", "paf_record"),
    ((chrom, start, end, f"merged_line:{line_number}")
     for chrom, start, end, _, line_number in unique_minimap2),
)

segtrace_both = sum(segtrace_matches)
minimap2_both = sum(minimap2_matches)
segtrace_only_count = len(segtrace) - segtrace_both
minimap2_unique_count = len({interval[:3] for interval in minimap2})
minimap2_only_count = minimap2_unique_count - len({interval[:3] for interval, matched in zip(minimap2, minimap2_matches) if matched})
segtrace_only_bp, both_bp, minimap2_only_bp = coverage_venn(
    segtrace_coverage, minimap2_coverage
)
total_bp = segtrace_only_bp + both_bp + minimap2_only_bp

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

SegTrace vs minimap2 base-pair coverage comparison
==================================================
Exact overlap of merged coordinates; interval matching thresholds do not apply

                         .-----------------------.
                    .---'                         '---.
                 .-'       SegTrace only: {segtrace_only_bp:>12,} bp      '-.
                /                                             \\
               /       Both: {both_bp:>12,} bp                 \\
               \\                                             /
                \\       minimap2 only: {minimap2_only_bp:>12,} bp       /
                 '-.                         .-'
                    '---.                 .---'
                        '-----------------------'

Base-pair totals:
    SegTrace covered: {segtrace_only_bp + both_bp:>12,} bp
    minimap2 covered: {both_bp + minimap2_only_bp:>12,} bp
    union:            {total_bp:>12,} bp

Input intervals:
    SegTrace merged intervals: {len(segtrace)}
    minimap2 merged query/target intervals: {len(minimap2)} ({minimap2_unique_count} unique coordinates)
    Comparison inputs: {seg_merged_path}, {minimap2_merged_path}

Lists:
  SegTrace-specific: {segtrace_only_path}
  minimap2-specific: {minimap2_only_path}
"""

with open(venn_path, "w") as handle:
    handle.write(venn)
print(venn, end="")
PY