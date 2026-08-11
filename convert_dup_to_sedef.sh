#!/usr/bin/env bash
# ==============================================================================
# convert_dup_to_sedef.sh
# Converts segtrace dup.bed format into SEDEF SD BED format (9-column BED with pair info)
# Usage: ./convert_dup_to_sedef.sh input.dup.bed [output_SD.bed]
# ==============================================================================

set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input.dup.bed> [output_SD.bed]"
    exit 1
fi

INPUT_FILE="$1"
if [ "$#" -ge 2 ]; then
    OUTPUT_FILE="$2"
else
    OUTPUT_FILE="${INPUT_FILE%.dup.bed}_SD.bed"
    if [ "$OUTPUT_FILE" = "$INPUT_FILE" ]; then
        OUTPUT_FILE="${INPUT_FILE}.sedef.bed"
    fi
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

python3 - "$INPUT_FILE" "$OUTPUT_FILE" << 'EOF'
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

clusters = {}

with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 4:
            continue
        chrom, start, end, cluster_id = parts[0], parts[1], parts[2], parts[3]
        clusters.setdefault(cluster_id, []).append((chrom, start, end))

with open(output_file, 'w') as out:
    for cid, regions in clusters.items():
        if len(regions) < 2:
            continue
        for i in range(len(regions)):
            for j in range(len(regions)):
                if i == j:
                    continue
                c1, s1, e1 = regions[i]
                c2, s2, e2 = regions[j]
                pair_info = f"{c2}:{s2}-{e2}"
                out.write(f"{c1}\t{s1}\t{e1}\t{pair_info}\t0\t+\t{s1}\t{e1}\t139,139,139\n")

EOF

echo "Successfully converted '$INPUT_FILE' -> '$OUTPUT_FILE'"
