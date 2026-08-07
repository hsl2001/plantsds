#!/usr/bin/env python3
"""
Segmental Duplication (SD) Comparison Script
Compares Segtrace output (.dup.bed) and SEDEF output (.bed) using bedtools after chromosome normalization and merging.
Also outputs a renamed BED file with standard chromosome names (chr1, chr2, ..., chrX, chrY).
"""

import sys
import os
import argparse
import subprocess
import numpy as np
import pandas as pd

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

def parse_chrom(c):
    if '-' in c:
        c = c.split('-')[-1]
    return ACC_MAP.get(c, c)

def create_renamed_bed(in_path, out_path):
    """Creates a new BED file preserving all columns, but with chromosome names converted to standard chr format."""
    count = 0
    with open(in_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if not line.strip():
                continue
            if line.startswith('#'):
                fout.write(line)
                continue
            parts = line.strip().split('\t')
            parts[0] = parse_chrom(parts[0])
            fout.write('\t'.join(parts) + '\n')
            count += 1
    return count

def normalize_segtrace(in_path, out_path):
    count = 0
    with open(in_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            c = parse_chrom(parts[0])
            s, e = parts[1], parts[2]
            fout.write(f"{c}\t{s}\t{e}\n")
            count += 1
    return count

def normalize_sedef(in_path, out_path):
    count = 0
    with open(in_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            c1, s1, e1 = parse_chrom(parts[0]), parts[1], parts[2]
            fout.write(f"{c1}\t{s1}\t{e1}\n")
            count += 1
            if len(parts) >= 12:
                c2, s2, e2 = parse_chrom(parts[9]), parts[10], parts[11]
                fout.write(f"{c2}\t{s2}\t{e2}\n")
                count += 1
    return count

def load_bed_stats(filepath):
    intervals = []
    chrom_bp = {}
    if not os.path.exists(filepath):
        return intervals, chrom_bp
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            c, s, e = parts[0], int(parts[1]), int(parts[2])
            l = e - s
            intervals.append((c, s, e, l))
            chrom_bp[c] = chrom_bp.get(c, 0) + l
    return intervals, chrom_bp

def calc_n50(lengths):
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    acc = 0
    for l in lengths:
        acc += l
        if acc >= total / 2.0:
            return l
    return 0

def run_comparison(segtrace_bed, sedef_bed, out_renamed=None, work_dir="cmp_out"):
    os.makedirs(work_dir, exist_ok=True)

    # 1. Create renamed BED file with standard chromosome names
    if not out_renamed:
        if segtrace_bed.endswith(".dup.bed"):
            out_renamed = segtrace_bed.replace(".dup.bed", ".chr.dup.bed")
        else:
            out_renamed = segtrace_bed + ".chr.bed"

    renamed_count = create_renamed_bed(segtrace_bed, out_renamed)
    print(f"[SUCCESS] Created renamed BED file ({renamed_count:,} records): {out_renamed}")
    
    st_norm = os.path.join(work_dir, "segtrace_norm.bed")
    st_sort = os.path.join(work_dir, "segtrace_sorted.bed")
    st_merged = os.path.join(work_dir, "segtrace_merged.bed")

    sd_norm = os.path.join(work_dir, "sedef_norm.bed")
    sd_sort = os.path.join(work_dir, "sedef_sorted.bed")
    sd_merged = os.path.join(work_dir, "sedef_merged.bed")

    intersect_file = os.path.join(work_dir, "intersect_merged.bed")
    st_unique_file = os.path.join(work_dir, "segtrace_unique_merged.bed")
    sd_unique_file = os.path.join(work_dir, "sedef_unique_merged.bed")

    print(f"[INFO] Normalizing Segtrace BED ({segtrace_bed})...")
    st_raw_count = normalize_segtrace(segtrace_bed, st_norm)
    
    print(f"[INFO] Normalizing SEDEF BED ({sedef_bed})...")
    sd_raw_count = normalize_sedef(sedef_bed, sd_norm)

    print("[INFO] Running bedtools sort & merge for Segtrace...")
    subprocess.run(f"bedtools sort -i {st_norm} > {st_sort}", shell=True, check=True)
    subprocess.run(f"bedtools merge -i {st_sort} > {st_merged}", shell=True, check=True)

    print("[INFO] Running bedtools sort & merge for SEDEF...")
    subprocess.run(f"bedtools sort -i {sd_norm} > {sd_sort}", shell=True, check=True)
    subprocess.run(f"bedtools merge -i {sd_sort} > {sd_merged}", shell=True, check=True)

    print("[INFO] Computing bedtools intersect, subtract, and jaccard...")
    subprocess.run(f"bedtools intersect -a {st_merged} -b {sd_merged} > {intersect_file}", shell=True, check=True)
    subprocess.run(f"bedtools subtract -a {st_merged} -b {sd_merged} > {st_unique_file}", shell=True, check=True)
    subprocess.run(f"bedtools subtract -a {sd_merged} -b {st_merged} > {sd_unique_file}", shell=True, check=True)

    jaccard_out = subprocess.check_output(f"bedtools jaccard -a {st_merged} -b {sd_merged}", shell=True).decode()

    # Parse metrics
    st_ivs, st_cbp = load_bed_stats(st_merged)
    sd_ivs, sd_cbp = load_bed_stats(sd_merged)
    is_ivs, is_cbp = load_bed_stats(intersect_file)
    st_u_ivs, st_u_cbp = load_bed_stats(st_unique_file)
    sd_u_ivs, sd_u_cbp = load_bed_stats(sd_unique_file)

    st_bp = sum(x[3] for x in st_ivs)
    sd_bp = sum(x[3] for x in sd_ivs)
    is_bp = sum(x[3] for x in is_ivs)
    st_u_bp = sum(x[3] for x in st_u_ivs)
    sd_u_bp = sum(x[3] for x in sd_u_ivs)

    recall = is_bp / sd_bp if sd_bp > 0 else 0.0
    precision = is_bp / st_bp if st_bp > 0 else 0.0
    f1 = 2 * recall * precision / (recall + precision) if (recall + precision) > 0 else 0.0
    jaccard = is_bp / (st_bp + sd_bp - is_bp) if (st_bp + sd_bp - is_bp) > 0 else 0.0

    print("\n==========================================================")
    print(" SEGMENTAL DUPLICATION COMPARISON REPORT (BEDTOOLS MERGE)")
    print("==========================================================")
    print(f"Renamed Output BED File: {out_renamed}")
    print(f"Segtrace Raw Records:    {st_raw_count:,}")
    print(f"SEDEF Raw Regions:      {sd_raw_count:,}")
    print("----------------------------------------------------------")
    print(f"Segtrace Merged Footprint: {st_bp:,} bp ({st_bp/1e6:.2f} Mb)")
    print(f"SEDEF Merged Footprint:    {sd_bp:,} bp ({sd_bp/1e6:.2f} Mb)")
    print(f"Intersection (Overlap):    {is_bp:,} bp ({is_bp/1e6:.2f} Mb)")
    print(f"Segtrace Unique Footprint: {st_u_bp:,} bp ({st_u_bp/1e6:.2f} Mb)")
    print(f"SEDEF Unique Footprint:    {sd_u_bp:,} bp ({sd_u_bp/1e6:.2f} Mb)")
    print("----------------------------------------------------------")
    print(f"Recall (vs SEDEF):         {recall*100:.2f}%")
    print(f"Precision (vs SEDEF):      {precision*100:.2f}%")
    print(f"F1-Score:                  {f1*100:.2f}%")
    print(f"Jaccard Index:             {jaccard:.6f}")
    print("==========================================================")

    # Print bedtools jaccard output
    print("\n[BEDTOOLS JACCARD RAW OUTPUT]")
    print(jaccard_out.strip())

    # Print Interval Length Statistics
    print("\n==========================================================")
    print(" INTERVAL LENGTH STATISTICS (AFTER BEDTOOLS MERGE)")
    print("==========================================================")
    print(f"{'Dataset':<20} | {'Count':<8} | {'Min (bp)':<9} | {'Max (bp)':<10} | {'Mean (bp)':<10} | {'Median (bp)':<11} | {'N50 (bp)':<10}")
    print("-" * 90)
    for name, ivs in [("Segtrace Merged", st_ivs), ("SEDEF Merged", sd_ivs), ("Intersection", is_ivs)]:
        lens = [x[3] for x in ivs]
        if lens:
            arr = np.array(lens)
            print(f"{name:<20} | {len(arr):<8,} | {arr.min():<9,} | {arr.max():<10,} | {arr.mean():<10.1f} | {np.median(arr):<11.1f} | {calc_n50(lens):<10,}")

    # Per-chromosome breakdown table
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    print("\n==========================================================")
    print(" PER-CHROMOSOME FOOTPRINT BREAKDOWN (Mb)")
    print("==========================================================")
    print(f"{'Chrom':<7} | {'Segtrace Merged':<15} | {'SEDEF Merged':<12} | {'Intersection':<12} | {'Segtrace Unique':<15} | {'SEDEF Unique':<12}")
    print("-" * 90)
    chrom_rows = []
    for c in chroms:
        st_c = st_cbp.get(c, 0) / 1e6
        sd_c = sd_cbp.get(c, 0) / 1e6
        is_c = is_cbp.get(c, 0) / 1e6
        st_u_c = st_u_cbp.get(c, 0) / 1e6
        sd_u_c = sd_u_cbp.get(c, 0) / 1e6
        print(f"{c:<7} | {st_c:15.2f} | {sd_c:12.2f} | {is_c:12.2f} | {st_u_c:15.2f} | {sd_u_c:12.2f}")
        chrom_rows.append({
            'Chrom': c,
            'Segtrace_Merged_Mb': st_c,
            'SEDEF_Merged_Mb': sd_c,
            'Intersection_Mb': is_c,
            'Segtrace_Unique_Mb': st_u_c,
            'SEDEF_Unique_Mb': sd_u_c
        })

    df_chrom = pd.DataFrame(chrom_rows)
    df_chrom.to_csv(os.path.join(work_dir, "chrom_breakdown.csv"), index=False)
    print(f"\n[INFO] Per-chromosome summary saved to {os.path.join(work_dir, 'chrom_breakdown.csv')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Segtrace and SEDEF BED files and save renamed BED with standard chromosome names.")
    parser.add_argument("--segtrace", default="t2t-chm13_sd.dup.bed", help="Path to Segtrace dup.bed file")
    parser.add_argument("--sedef", default="sedef-human-t2tchm13.bed", help="Path to SEDEF bed file")
    parser.add_argument("--out-renamed", default=None, help="Output path for renamed Segtrace BED file (default: <segtrace_base>.chr.dup.bed)")
    parser.add_argument("--work-dir", default="cmp_out", help="Directory to store intermediate bed files")
    args = parser.parse_args()

    if not os.path.exists(args.segtrace):
        print(f"[ERROR] Segtrace file '{args.segtrace}' not found.")
        sys.exit(1)
    if not os.path.exists(args.sedef):
        print(f"[ERROR] SEDEF file '{args.sedef}' not found.")
        sys.exit(1)

    run_comparison(args.segtrace, args.sedef, out_renamed=args.out_renamed, work_dir=args.work_dir)
