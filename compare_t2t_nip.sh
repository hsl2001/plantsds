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
#   <prefix>.compare.evidence.tsv
#   <prefix>.compare.statistics.txt
#
# Tool and matching parameters can be adjusted without editing this script:
#   THREADS=36 MAX_MAPPINGS=1000000 ./compare_t2t_nip.sh t2t_nip.fasta t2t_nip
# Set FORCE_MINIMAP2=1 to regenerate an existing PAF.
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
FORCE_MINIMAP2=${FORCE_MINIMAP2:-0}
MIN_OVERLAP_FRAC=${MIN_OVERLAP_FRAC:-0.5}
MIN_OVERLAP_BP=${MIN_OVERLAP_BP:-100}
MIN_ALIGN_LEN=${MIN_ALIGN_LEN:-1000}
EVIDENCE_MIN_FRAC=${EVIDENCE_MIN_FRAC:-0.1}
EVIDENCE_MIN_BP=${EVIDENCE_MIN_BP:-100}
EVIDENCE_MIN_IDENTITY=${EVIDENCE_MIN_IDENTITY:-0.8}
N_PERMUTATIONS=${N_PERMUTATIONS:-1000}
RANDOM_SEED=${RANDOM_SEED:-42}

if [[ ! -r "$FASTA" ]]; then
    printf 'error: FASTA not found or not readable: %s\n' "$FASTA" >&2
    exit 1
fi

if [[ ! -x "$SEGTRACE_BIN" ]]; then
    printf 'error: SegTrace executable not found or not executable: %s\n' "$SEGTRACE_BIN" >&2
    exit 1
fi
if [[ ! -s "$PAF" || "$FORCE_MINIMAP2" == 1 ]]; then
    if ! command -v "$MINIMAP2_BIN" >/dev/null 2>&1; then
        printf 'error: minimap2 executable not found in PATH: %s\n' "$MINIMAP2_BIN" >&2
        exit 1
    fi
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

if [[ -s "$PAF" && "$FORCE_MINIMAP2" != 1 ]]; then
    printf '[2/4] Using existing minimap2 PAF: %s\n' "$PAF" >&2
else
    # -P retains all chains, while -p 0 and a large -N keep low-scoring and
    # numerous secondary mappings.
    printf '[2/4] Running minimap2 (%s; all mappings)...\n' "$MINIMAP2_PRESET" >&2
    "$MINIMAP2_BIN" -x "$MINIMAP2_PRESET" --secondary=yes -P -p 0 \
        -N "$MAX_MAPPINGS" -t "$THREADS" "$FASTA" "$FASTA" > "$PAF"
fi
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
python3 - "$PAF" "$SEG_MERGED" "$MINIMAP2_MERGED" "$COMPARE_PREFIX" \
    "$MIN_OVERLAP_FRAC" "$MIN_OVERLAP_BP" "$MIN_ALIGN_LEN" \
    "$EVIDENCE_MIN_FRAC" "$EVIDENCE_MIN_BP" "$EVIDENCE_MIN_IDENTITY" \
    "$N_PERMUTATIONS" "$RANDOM_SEED" <<'PY'
import bisect
import math
import random
import sys
from collections import defaultdict


(
    paf_path, seg_merged_path, minimap2_merged_path, prefix,
    min_frac_s, min_bp_s, min_align_len_s,
    evidence_frac_s, evidence_bp_s, evidence_identity_s,
    permutations_s, seed_s,
) = sys.argv[1:]
min_frac = float(min_frac_s)
min_bp = int(min_bp_s)
min_align_len = int(min_align_len_s)
evidence_frac = float(evidence_frac_s)
evidence_bp = int(evidence_bp_s)
evidence_identity = float(evidence_identity_s)
n_permutations = int(permutations_s)
random_seed = int(seed_s)

if not 0 < min_frac <= 1:
    raise SystemExit("error: MIN_OVERLAP_FRAC must be in (0, 1]")
if min_bp < 1 or min_align_len < 1:
    raise SystemExit("error: MIN_OVERLAP_BP and MIN_ALIGN_LEN must be positive")
if not 0 < evidence_frac <= 1 or evidence_bp < 1 or not 0 < evidence_identity <= 1:
    raise SystemExit("error: invalid relaxed PAF evidence thresholds")
if n_permutations < 1:
    raise SystemExit("error: N_PERMUTATIONS must be positive")


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


def index_intervals(merged):
    return {
        chrom: ([(start, end) for start, end in values], [start for start, _ in values])
        for chrom, values in merged.items()
    }


def interval_matches(interval, indexed_other):
    chrom, start, end = interval[:3]
    indexed = indexed_other.get(chrom)
    if not indexed:
        return False

    candidates, starts = indexed
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


def read_paf_evidence(path):
    endpoints = defaultdict(list)
    chrom_lengths = {}
    skipped_diagonal = 0
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                query, qlen = fields[0], int(fields[1])
                qstart, qend = int(fields[2]), int(fields[3])
                target, tlen = fields[5], int(fields[6])
                tstart, tend = int(fields[7]), int(fields[8])
                matches, block_length = int(fields[9]), int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue
            if block_length <= 0 or qend <= qstart or tend <= tstart:
                continue
            chrom_lengths[query] = max(chrom_lengths.get(query, 0), qlen)
            chrom_lengths[target] = max(chrom_lengths.get(target, 0), tlen)
            if query == target and qstart == tstart and qend == tend:
                skipped_diagonal += 1
                continue
            identity = matches / block_length
            evidence = (identity, block_length, mapq, f"{line_number}:query", target)
            endpoints[query].append((qstart, qend, evidence))
            evidence = (identity, block_length, mapq, f"{line_number}:target", query)
            endpoints[target].append((tstart, tend, evidence))
    return endpoints, chrom_lengths, skipped_diagonal


def best_paf_support(interval, endpoints):
    chrom, start, end = interval[:3]
    interval_length = end - start
    best = None
    for hit_start, hit_end, evidence in endpoints.get(chrom, ()):
        overlap = min(end, hit_end) - max(start, hit_start)
        if overlap <= 0:
            continue
        overlap_fraction = overlap / interval_length
        identity, block_length, mapq, record_id, other_chrom = evidence
        candidate = (overlap_fraction, identity, overlap, block_length, mapq, record_id, other_chrom)
        if best is None or candidate[:5] > best[:5]:
            best = candidate
    if best is None:
        return (0, 0.0, 0.0, 0, 0, ".", ".", False)
    overlap_fraction, identity, overlap, block_length, mapq, record_id, other_chrom = best
    supported = (
        overlap >= evidence_bp
        and overlap_fraction >= evidence_frac
        and identity >= evidence_identity
    )
    return (overlap, overlap_fraction, identity, block_length, mapq,
            record_id, other_chrom, supported)


def wilson_interval(successes, total):
    if total == 0:
        return (0.0, 0.0)
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = (proportion + z * z / (2 * total)) / denominator
    margin = z * math.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return (max(0.0, centre - margin), min(1.0, centre + margin))


def random_interval(interval, chrom_lengths, rng):
    chrom, start, end = interval[:3]
    length = end - start
    chrom_length = chrom_lengths.get(chrom, end)
    max_start = max(0, chrom_length - length)
    random_start = rng.randint(0, max_start) if max_start else 0
    return (chrom, random_start, random_start + length, ".", 0)


def permutation_overlap_test(intervals, target_index, chrom_lengths, rng):
    observed = sum(interval_matches(interval, target_index) for interval in intervals)
    null_counts = []
    for _ in range(n_permutations):
        null_counts.append(sum(
            interval_matches(random_interval(interval, chrom_lengths, rng), target_index)
            for interval in intervals
        ))
    mean = sum(null_counts) / n_permutations
    variance = sum((count - mean) ** 2 for count in null_counts) / n_permutations
    p_high = (1 + sum(count >= observed for count in null_counts)) / (n_permutations + 1)
    p_low_specific = (1 + sum(
        len(intervals) - count <= len(intervals) - observed for count in null_counts
    )) / (n_permutations + 1)
    return observed, mean, math.sqrt(variance), p_high, p_low_specific


def write_tsv(path, header, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


segtrace = read_bed(seg_merged_path)
minimap2 = read_bed(minimap2_merged_path)
segtrace_coverage = merge_intervals(segtrace)
minimap2_coverage = merge_intervals(minimap2)
segtrace_index = index_intervals(segtrace_coverage)
minimap2_index = index_intervals(minimap2_coverage)
endpoints, chrom_lengths, skipped_diagonal = read_paf_evidence(paf_path)

segtrace_matches = [interval_matches(interval, minimap2_index) for interval in segtrace]
minimap2_matches = [interval_matches(interval, segtrace_index) for interval in minimap2]

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
evidence_path = f"{prefix}.evidence.tsv"
statistics_path = f"{prefix}.statistics.txt"

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

segtrace_evidence = [best_paf_support(interval, endpoints) for interval in segtrace_specific]
minimap2_evidence = [best_paf_support(interval, endpoints) for interval in minimap2_specific]
all_chrom_lengths = dict(chrom_lengths)
for interval in segtrace + minimap2:
    all_chrom_lengths[interval[0]] = max(all_chrom_lengths.get(interval[0], 0), interval[2])

rng = random.Random(random_seed)
seg_observed, seg_null_mean, seg_null_sd, seg_p_high, seg_p_low_specific = permutation_overlap_test(
    segtrace, minimap2_index, all_chrom_lengths, rng
)
min_observed, min_null_mean, min_null_sd, min_p_high, min_p_low_specific = permutation_overlap_test(
    minimap2, segtrace_index, all_chrom_lengths, rng
)

seg_supported = sum(evidence[-1] for evidence in segtrace_evidence)
min_supported = sum(evidence[-1] for evidence in minimap2_evidence)
seg_ci_low, seg_ci_high = wilson_interval(seg_supported, len(segtrace_specific))
min_ci_low, min_ci_high = wilson_interval(min_supported, len(minimap2_specific))

write_tsv(
    evidence_path,
    ("set", "chrom", "start", "end", "length", "cluster_id", "callset_line",
     "raw_paf_support", "overlap_bp", "overlap_fraction", "identity",
     "alignment_block_bp", "mapq", "other_chrom", "paf_record", "evidence_label"),
    (
        (label, chrom, start, end, end - start, cluster, line_number,
         "yes" if evidence[-1] else "no", evidence[0], f"{evidence[1]:.6f}",
         f"{evidence[2]:.6f}", evidence[3], evidence[4], evidence[6], evidence[5],
         "strong_raw_paf_support" if evidence[-1] else "unresolved")
        for label, intervals, evidence_rows in (
            ("segtrace-specific", segtrace_specific, segtrace_evidence),
            ("minimap2-specific", minimap2_specific, minimap2_evidence),
        )
        for (chrom, start, end, cluster, line_number), evidence in zip(intervals, evidence_rows)
    ),
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
    SegTrace merged intervals: {len(segtrace)}
    minimap2 merged query/target intervals: {len(minimap2)} ({minimap2_unique_count} unique coordinates)
    Comparison inputs: {seg_merged_path}, {minimap2_merged_path}

Statistical evidence:
    raw non-diagonal PAF records used: {sum(len(values) for values in endpoints.values()) // 2}
    exact self-diagonal PAF records skipped: {skipped_diagonal}
    SegTrace-specific strong raw-PAF support: {seg_supported}/{len(segtrace_specific)} ({seg_supported / len(segtrace_specific) if segtrace_specific else 0:.3f}; Wilson 95% CI {seg_ci_low:.3f}-{seg_ci_high:.3f})
    minimap2-specific strong raw-PAF support: {min_supported}/{len(minimap2_specific)} ({min_supported / len(minimap2_specific) if minimap2_specific else 0:.3f}; Wilson 95% CI {min_ci_low:.3f}-{min_ci_high:.3f})
    SegTrace overlap count: observed {seg_observed}, random mean {seg_null_mean:.1f} +/- {seg_null_sd:.1f}, permutation p(overlap >= observed)={seg_p_high:.4g}
    minimap2 overlap count: observed {min_observed}, random mean {min_null_mean:.1f} +/- {min_null_sd:.1f}, permutation p(overlap >= observed)={min_p_high:.4g}
    SegTrace-specific count: observed {segtrace_only_count}, random mean {len(segtrace) - seg_null_mean:.1f}, permutation p(specific <= observed)={seg_p_low_specific:.4g}
    minimap2-specific count: observed {minimap2_only_count}, random mean {len(minimap2) - min_null_mean:.1f}, permutation p(specific <= observed)={min_p_low_specific:.4g}
    Randomization: {n_permutations} placements, same chromosome and interval length, seed {random_seed}

Interpretation limits:
    Strong raw-PAF support argues against a call being an arbitrary interval, but PAF is not an independent truth set.
    Permutation p-values test non-random overlap, not biological truth and not FP/FN status by themselves.
    Definitive FP/FN claims require an independent truth set or an independent validation method.

Lists:
  SegTrace-specific: {segtrace_only_path}
  minimap2-specific: {minimap2_only_path}
    Per-interval evidence: {evidence_path}
"""

with open(venn_path, "w") as handle:
    handle.write(venn)

with open(statistics_path, "w") as handle:
    handle.write(venn)
print(venn, end="")
PY