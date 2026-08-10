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

    # Fragment overlap files (Merged)
    st_m_hit_sd = os.path.join(work_dir, "segtrace_merged_hit_sedef.bed")
    sd_m_hit_st = os.path.join(work_dir, "sedef_merged_hit_segtrace.bed")

    # Fragment overlap files (Raw)
    st_raw_hit_sd = os.path.join(work_dir, "segtrace_raw_hit_sedef.bed")
    sd_raw_hit_st = os.path.join(work_dir, "sedef_raw_hit_segtrace.bed")

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

    print("[INFO] Computing BP-level intersections and subtractions...")
    subprocess.run(f"bedtools intersect -a {st_merged} -b {sd_merged} > {intersect_file}", shell=True, check=True)
    subprocess.run(f"bedtools subtract -a {st_merged} -b {sd_merged} > {st_unique_file}", shell=True, check=True)
    subprocess.run(f"bedtools subtract -a {sd_merged} -b {st_merged} > {sd_unique_file}", shell=True, check=True)

    print("[INFO] Computing Fragment-level overlaps (Merged & Raw)...")
    # Merged Fragment Overlap
    subprocess.run(f"bedtools intersect -a {st_merged} -b {sd_merged} -u > {st_m_hit_sd}", shell=True, check=True)
    subprocess.run(f"bedtools intersect -a {sd_merged} -b {st_merged} -u > {sd_m_hit_st}", shell=True, check=True)
    
    # Raw Fragment Overlap
    subprocess.run(f"bedtools intersect -a {st_norm} -b {sd_norm} -u > {st_raw_hit_sd}", shell=True, check=True)
    subprocess.run(f"bedtools intersect -a {sd_norm} -b {st_norm} -u > {sd_raw_hit_st}", shell=True, check=True)

    jaccard_out = subprocess.check_output(f"bedtools jaccard -a {st_merged} -b {sd_merged}", shell=True).decode()

    # Load interval data
    st_raw_ivs, _ = load_bed_stats(st_norm)
    sd_raw_ivs, _ = load_bed_stats(sd_norm)

    st_m_ivs, st_cbp = load_bed_stats(st_merged)
    sd_m_ivs, sd_cbp = load_bed_stats(sd_merged)
    is_ivs, is_cbp = load_bed_stats(intersect_file)
    st_u_ivs, st_u_cbp = load_bed_stats(st_unique_file)
    sd_u_ivs, sd_u_cbp = load_bed_stats(sd_unique_file)

    st_m_hit_ivs, _ = load_bed_stats(st_m_hit_sd)
    sd_m_hit_ivs, _ = load_bed_stats(sd_m_hit_st)

    st_raw_hit_ivs, _ = load_bed_stats(st_raw_hit_sd)
    sd_raw_hit_ivs, _ = load_bed_stats(sd_raw_hit_st)

    # 1. Base-Pair (BP) Level Metrics
    st_bp = sum(x[3] for x in st_m_ivs)
    sd_bp = sum(x[3] for x in sd_m_ivs)
    is_bp = sum(x[3] for x in is_ivs)
    st_u_bp = sum(x[3] for x in st_u_ivs)
    sd_u_bp = sum(x[3] for x in sd_u_ivs)

    bp_recall = is_bp / sd_bp if sd_bp > 0 else 0.0
    bp_precision = is_bp / st_bp if st_bp > 0 else 0.0
    bp_f1 = 2 * bp_recall * bp_precision / (bp_recall + bp_precision) if (bp_recall + bp_precision) > 0 else 0.0
    bp_jaccard = is_bp / (st_bp + sd_bp - is_bp) if (st_bp + sd_bp - is_bp) > 0 else 0.0

    # 2. Merged Fragment (Frag) Level Metrics
    st_m_count = len(st_m_ivs)
    sd_m_count = len(sd_m_ivs)
    st_m_hit_count = len(st_m_hit_ivs)
    sd_m_hit_count = len(sd_m_hit_ivs)

    m_frag_precision = st_m_hit_count / st_m_count if st_m_count > 0 else 0.0
    m_frag_recall = sd_m_hit_count / sd_m_count if sd_m_count > 0 else 0.0
    m_frag_f1 = 2 * m_frag_recall * m_frag_precision / (m_frag_recall + m_frag_precision) if (m_frag_recall + m_frag_precision) > 0 else 0.0

    # 3. Raw Fragment (Frag) Level Metrics
    st_r_count = len(st_raw_ivs)
    sd_r_count = len(sd_raw_ivs)
    st_r_hit_count = len(st_raw_hit_ivs)
    sd_r_hit_count = len(sd_raw_hit_ivs)

    r_frag_precision = st_r_hit_count / st_r_count if st_r_count > 0 else 0.0
    r_frag_recall = sd_r_hit_count / sd_r_count if sd_r_count > 0 else 0.0
    r_frag_f1 = 2 * r_frag_recall * r_frag_precision / (r_frag_recall + r_frag_precision) if (r_frag_recall + r_frag_precision) > 0 else 0.0

    # Output Comprehensive Evaluation Report
    print("\n=================================================================================")
    print("      SEGMENTAL DUPLICATION COMPARISON REPORT: BP & FRAGMENT LEVEL METRICS")
    print("=================================================================================")
    print(f"Renamed Segtrace BED:    {out_renamed}")
    print(f"Working Directory:       {work_dir}")
    print("---------------------------------------------------------------------------------")
    print(" 1. BASE-PAIR (BP) LEVEL METRICS (Genomic Footprint)")
    print("---------------------------------------------------------------------------------")
    print(f"  Segtrace Merged Footprint:  {st_bp:12,} bp ({st_bp/1e6:8.2f} Mb)")
    print(f"  SEDEF Merged Footprint:     {sd_bp:12,} bp ({sd_bp/1e6:8.2f} Mb)")
    print(f"  Overlap (Intersection) BP:  {is_bp:12,} bp ({is_bp/1e6:8.2f} Mb)")
    print(f"  Segtrace Unique Footprint:  {st_u_bp:12,} bp ({st_u_bp/1e6:8.2f} Mb)")
    print(f"  SEDEF Unique Footprint:     {sd_u_bp:12,} bp ({sd_u_bp/1e6:8.2f} Mb)")
    print(f"  -------------------------------------------------------------------------------")
    print(f"  BP Sensitivity / Recall:    {bp_recall*100:8.2f}%")
    print(f"  BP Precision:               {bp_precision*100:8.2f}%")
    print(f"  BP F1-Score:                {bp_f1*100:8.2f}%")
    print(f"  BP Jaccard Index:           {bp_jaccard:10.6f}")

    print("\n---------------------------------------------------------------------------------")
    print(" 2. FRAGMENT (FRAG) LEVEL METRICS: MERGED INTERVALS")
    print("---------------------------------------------------------------------------------")
    print(f"  Segtrace Merged Frags:      {st_m_count:12,} total frags")
    print(f"  SEDEF Merged Frags:         {sd_m_count:12,} total frags")
    print(f"  Segtrace Frags Overlap SEDEF: {st_m_hit_count:10,} / {st_m_count:,} ({m_frag_precision*100:6.2f}% Precision)")
    print(f"  SEDEF Frags Overlap Segtrace: {sd_m_hit_count:10,} / {sd_m_count:,} ({m_frag_recall*100:6.2f}% Recall)")
    print(f"  Merged Frag F1-Score:       {m_frag_f1*100:8.2f}%")

    print("\n---------------------------------------------------------------------------------")
    print(" 3. FRAGMENT (FRAG) LEVEL METRICS: RAW INTERVALS")
    print("---------------------------------------------------------------------------------")
    print(f"  Segtrace Raw Records:       {st_r_count:12,} raw records")
    print(f"  SEDEF Raw Records:          {sd_r_count:12,} raw records")
    print(f"  Segtrace Raw Hit SEDEF:     {st_r_hit_count:10,} / {st_r_count:,} ({r_frag_precision*100:6.2f}% Precision)")
    print(f"  SEDEF Raw Hit Segtrace:     {sd_r_hit_count:10,} / {sd_r_count:,} ({r_frag_recall*100:6.2f}% Recall)")
    print(f"  Raw Frag F1-Score:          {r_frag_f1*100:8.2f}%")
    print("=================================================================================")

    # Print bedtools jaccard output
    print("\n[BEDTOOLS JACCARD RAW OUTPUT]")
    print(jaccard_out.strip())

    # Print Detailed Fragment Length Statistics Table
    print("\n=================================================================================")
    print(" FRAGMENT LENGTH & COUNT STATISTICS")
    print("=================================================================================")
    print(f"{'Dataset':<24} | {'Count':<8} | {'Min (bp)':<9} | {'Max (bp)':<10} | {'Mean (bp)':<10} | {'Median (bp)':<11} | {'N50 (bp)':<10}")
    print("-" * 95)
    stat_datasets = [
        ("Segtrace Raw", st_raw_ivs),
        ("SEDEF Raw", sd_raw_ivs),
        ("Segtrace Merged", st_m_ivs),
        ("SEDEF Merged", sd_m_ivs),
        ("Intersection (Overlap)", is_ivs),
        ("Segtrace Unique", st_u_ivs),
        ("SEDEF Unique", sd_u_ivs)
    ]
    for name, ivs in stat_datasets:
        lens = [x[3] for x in ivs]
        if lens:
            arr = np.array(lens)
            print(f"{name:<24} | {len(arr):<8,} | {arr.min():<9,} | {arr.max():<10,} | {arr.mean():<10.1f} | {np.median(arr):<11.1f} | {calc_n50(lens):<10,}")
        else:
            print(f"{name:<24} | {0:<8} | {'N/A':<9} | {'N/A':<10} | {'N/A':<10} | {'N/A':<11} | {'N/A':<10}")

    # Per-chromosome breakdown table
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    print("\n=================================================================================")
    print(" PER-CHROMOSOME FOOTPRINT BREAKDOWN (Mb)")
    print("=================================================================================")
    print(f"{'Chrom':<7} | {'Segtrace Merged':<15} | {'SEDEF Merged':<12} | {'Intersection':<12} | {'Segtrace Unique':<15} | {'SEDEF Unique':<12}")
    print("-" * 95)
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
    summary_path = os.path.join(work_dir, "chrom_breakdown.csv")
    df_chrom.to_csv(summary_path, index=False)

    # Save summary metrics to CSV as well
    metrics_summary = pd.DataFrame([{
        'ST_Raw_Records': st_r_count,
        'SD_Raw_Records': sd_r_count,
        'ST_Merged_BP': st_bp,
        'SD_Merged_BP': sd_bp,
        'Intersection_BP': is_bp,
        'BP_Recall': bp_recall,
        'BP_Precision': bp_precision,
        'BP_F1': bp_f1,
        'BP_Jaccard': bp_jaccard,
        'ST_Merged_Frags': st_m_count,
        'SD_Merged_Frags': sd_m_count,
        'ST_Merged_Frag_Hit': st_m_hit_count,
        'SD_Merged_Frag_Hit': sd_m_hit_count,
        'Merged_Frag_Precision': m_frag_precision,
        'Merged_Frag_Recall': m_frag_recall,
        'Merged_Frag_F1': m_frag_f1,
        'ST_Raw_Frag_Hit': st_r_hit_count,
        'SD_Raw_Frag_Hit': sd_r_hit_count,
        'Raw_Frag_Precision': r_frag_precision,
        'Raw_Frag_Recall': r_frag_recall,
        'Raw_Frag_F1': r_frag_f1
    }])
    metrics_path = os.path.join(work_dir, "evaluation_metrics.csv")
    metrics_summary.to_csv(metrics_path, index=False)

    print(f"\n[SUCCESS] Evaluation complete!")
    print(f" - Per-chromosome breakdown: {summary_path}")
    print(f" - Summary metrics CSV:      {metrics_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Segtrace and SEDEF BED files at both BP and Fragment levels.")
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
