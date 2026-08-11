#!/usr/bin/env bash
# ==============================================================================
# convert_dup_to_sedef.sh
# Converts segtrace dup.bed format (e.g., t2t-chm13_sd.dup.bed) into SEDEF SD BED format
# using accession mapping (ACC_MAP) to convert RefSeq IDs to standard chr names.
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

ACC_MAP = {
    'NC_060925.1': 'chr1', 'NC_060926.1': 'chr2', 'NC_060927.1': 'chr3',
    'NC_060928.1': 'chr4', 'NC_060929.1': 'chr5', 'NC_060930.1': 'chr6',
    'NC_060931.1': 'chr7', 'NC_060932.1': 'chr8', 'NC_060933.1': 'chr9',
    'NC_060934.1': 'chr10', 'NC_060935.1': 'chr11', 'NC_060936.1': 'chr12',
    'NC_060937.1': 'chr13', 'NC_060938.1': 'chr14', 'NC_060939.1': 'chr15',
    'NC_060940.1': 'chr16', 'NC_060941.1': 'chr17', 'NC_060942.1': 'chr18',
    'NC_060943.1': 'chr19', 'NC_060944.1': 'chr20', 'NC_060945.1': 'chr21',
    'NC_060946.1': 'chr22', 'NC_060947.1': 'chrX', 'NC_060948.1': 'chrY',
    'NC_012920.1': 'chrM'
}

def format_chrom(c):
    if '-' in c and not c.startswith('chr'):
        c = c.split('-')[-1]
    return ACC_MAP.get(c, c)

clusters = {}

with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 4:
            continue
        chrom, start, end, cluster_id = format_chrom(parts[0]), parts[1], parts[2], parts[3]
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
