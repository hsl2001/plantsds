#!/usr/bin/env python3
"""
Segmental Duplication (SD) Comparison Script
Compares Segtrace output (.dup.bed) and SEDEF output (.bed).
Evaluates strictly at:
  1) Base-Pair (BP) level (Genomic Footprint Sensitivity, Precision, F1, Jaccard)
  2) Reciprocal Overlap Fragment (Frag) level (evaluate.py algorithm at >=50% and >=10% threshold)
All temporary intermediate files are cleaned up automatically.
"""

import sys
import os
import argparse
import subprocess
import shutil
import numpy as np

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
    """Creates a new BED file preserving all columns with standardized chromosome names."""
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

def load_bed_bp(filepath):
    total_bp = 0
    if not os.path.exists(filepath):
        return 0
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            s, e = int(parts[1]), int(parts[2])
            total_bp += max(0, e - s)
    return total_bp

def load_segtrace_pairs(in_path):
    """Loads paired SD regions from Segtrace .dup.bed file based on cluster_id."""
    clusters = {}
    if not os.path.exists(in_path):
        return []
    with open(in_path) as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 4:
                c, s, e, cid = parse_chrom(parts[0]), int(parts[1]), int(parts[2]), parts[3]
                if cid not in clusters:
                    clusters[cid] = []
                clusters[cid].append((c, s, e))
    
    pairs = []
    for cid, locs in clusters.items():
        for i in range(len(locs)):
            for j in range(i + 1, len(locs)):
                pairs.append((locs[i], locs[j]))
    return pairs

def load_sedef_pairs(in_path):
    """Loads paired SD regions from SEDEF .bed file (columns 1-3 vs columns 10-12)."""
    pairs = []
    if not os.path.exists(in_path):
        return pairs
    with open(in_path) as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 12:
                c1, s1, e1 = parse_chrom(parts[0]), int(parts[1]), int(parts[2])
                c2, s2, e2 = parse_chrom(parts[9]), int(parts[10]), int(parts[11])
                pairs.append(((c1, s1, e1), (c2, s2, e2)))
    return pairs

def evaluate_frag_pairs(ref_pairs, target_pairs, threshold=0.5):
    """
    Evaluates pair-level fragment metrics using reciprocal overlap threshold (from evaluate.py).
    ref_pairs: list of ((c1, s1, e1), (c2, s2, e2))
    target_pairs: list of ((c1, s1, e1), (c2, s2, e2))
    """
    if not ref_pairs or not target_pairs:
        return 0.0, 0.0, 0.0, 0, 0

    # Group pairs by chromosome pair to speed up matching
    ref_by_chrom = {}
    for (c1, s1, e1), (c2, s2, e2) in ref_pairs:
        key = (c1, c2) if c1 <= c2 else (c2, c1)
        if key not in ref_by_chrom:
            ref_by_chrom[key] = []
        if c1 <= c2:
            ref_by_chrom[key].append(((s1, e1), (s2, e2)))
        else:
            ref_by_chrom[key].append(((s2, e2), (s1, e1)))

    target_by_chrom = {}
    for (c1, s1, e1), (c2, s2, e2) in target_pairs:
        key = (c1, c2) if c1 <= c2 else (c2, c1)
        if key not in target_by_chrom:
            target_by_chrom[key] = []
        if c1 <= c2:
            target_by_chrom[key].append(((s1, e1), (s2, e2)))
        else:
            target_by_chrom[key].append(((s2, e2), (s1, e1)))

    matched_ref = set()
    matched_target = set()
    ref_id = 0
    
    for key, r_list in ref_by_chrom.items():
        if key not in target_by_chrom:
            ref_id += len(r_list)
            continue
        
        t_list = target_by_chrom[key]
        
        PX1 = np.array([p[0][0] for p in t_list], dtype=np.float64)
        PY1 = np.array([p[0][1] for p in t_list], dtype=np.float64)
        PX2 = np.array([p[1][0] for p in t_list], dtype=np.float64)
        PY2 = np.array([p[1][1] for p in t_list], dtype=np.float64)
        LP1 = PY1 - PX1
        LP2 = PY2 - PX2

        for ((s1, e1), (s2, e2)) in r_list:
            Lt1 = float(e1 - s1)
            Lt2 = float(e2 - s2)
            if Lt1 <= 0 or Lt2 <= 0:
                ref_id += 1
                continue

            # Direct orientation overlaps
            o1 = np.maximum(0.0, np.minimum(e1, PY1) - np.maximum(s1, PX1))
            o2 = np.maximum(0.0, np.minimum(e2, PY2) - np.maximum(s2, PX2))

            dir_match = (
                (o1 / Lt1 >= threshold) & (o1 / LP1 >= threshold) &
                (o2 / Lt2 >= threshold) & (o2 / LP2 >= threshold)
            )

            # Cross orientation overlaps
            o12 = np.maximum(0.0, np.minimum(e1, PY2) - np.maximum(s1, PX2))
            o21 = np.maximum(0.0, np.minimum(e2, PY1) - np.maximum(s2, PX1))

            cross_match = (
                (o12 / Lt1 >= threshold) & (o12 / LP2 >= threshold) &
                (o21 / Lt2 >= threshold) & (o21 / LP1 >= threshold)
            )

            matched_idx = np.where(dir_match | cross_match)[0]
            if len(matched_idx) > 0:
                matched_ref.add(ref_id)
                for idx in matched_idx:
                    matched_target.add((key, idx))
            ref_id += 1

    total_ref = len(ref_pairs)
    total_target = len(target_pairs)
    
    Sn = len(matched_ref) / total_ref if total_ref > 0 else 0.0
    Pr = len(matched_target) / total_target if total_target > 0 else 0.0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0.0
    return Sn, Pr, f1, len(matched_ref), len(matched_target)

def run_comparison(segtrace_bed, sedef_bed, out_renamed=None, work_dir="_cmp_tmp", keep_temp=False):
    os.makedirs(work_dir, exist_ok=True)

    try:
        # Create renamed BED file if requested
        if out_renamed:
            create_renamed_bed(segtrace_bed, out_renamed)

        st_norm = os.path.join(work_dir, "segtrace_norm.bed")
        st_sort = os.path.join(work_dir, "segtrace_sorted.bed")
        st_merged = os.path.join(work_dir, "segtrace_merged.bed")

        sd_norm = os.path.join(work_dir, "sedef_norm.bed")
        sd_sort = os.path.join(work_dir, "sedef_sorted.bed")
        sd_merged = os.path.join(work_dir, "sedef_merged.bed")

        intersect_file = os.path.join(work_dir, "intersect_merged.bed")
        st_unique_file = os.path.join(work_dir, "segtrace_unique_merged.bed")
        sd_unique_file = os.path.join(work_dir, "sedef_unique_merged.bed")

        normalize_segtrace(segtrace_bed, st_norm)
        normalize_sedef(sedef_bed, sd_norm)

        # bedtools sort & merge
        subprocess.run(f"bedtools sort -i {st_norm} > {st_sort}", shell=True, check=True)
        subprocess.run(f"bedtools merge -i {st_sort} > {st_merged}", shell=True, check=True)

        subprocess.run(f"bedtools sort -i {sd_norm} > {sd_sort}", shell=True, check=True)
        subprocess.run(f"bedtools merge -i {sd_sort} > {sd_merged}", shell=True, check=True)

        # BP-level intersections and subtractions
        subprocess.run(f"bedtools intersect -a {st_merged} -b {sd_merged} > {intersect_file}", shell=True, check=True)
        subprocess.run(f"bedtools subtract -a {st_merged} -b {sd_merged} > {st_unique_file}", shell=True, check=True)
        subprocess.run(f"bedtools subtract -a {sd_merged} -b {st_merged} > {sd_unique_file}", shell=True, check=True)

        # 1. Base-Pair (BP) Level Metrics
        st_bp = load_bed_bp(st_merged)
        sd_bp = load_bed_bp(sd_merged)
        is_bp = load_bed_bp(intersect_file)
        st_u_bp = load_bed_bp(st_unique_file)
        sd_u_bp = load_bed_bp(sd_unique_file)

        bp_recall = is_bp / sd_bp if sd_bp > 0 else 0.0
        bp_precision = is_bp / st_bp if st_bp > 0 else 0.0
        bp_f1 = 2 * bp_recall * bp_precision / (bp_recall + bp_precision) if (bp_recall + bp_precision) > 0 else 0.0
        bp_jaccard = is_bp / (st_bp + sd_bp - is_bp) if (st_bp + sd_bp - is_bp) > 0 else 0.0

        # 2. Fragment (Frag) Pair-Level Evaluation (evaluate.py algorithm at 50% reciprocal overlap)
        st_pairs = load_segtrace_pairs(segtrace_bed)
        sd_pairs = load_sedef_pairs(sedef_bed)

        pair_sn_50, pair_pr_50, pair_f1_50, m_ref_50, m_tgt_50 = evaluate_frag_pairs(sd_pairs, st_pairs, threshold=0.5)

        # Output Summary Report
        print("=================================================================================")
        print("            SEGMENTAL DUPLICATION COMPARISON REPORT (BP & FRAG LEVEL)")
        print("=================================================================================")
        print(f"Segtrace Input:  {segtrace_bed}")
        print(f"SEDEF Input:     {sedef_bed}")
        print("---------------------------------------------------------------------------------")
        print(" 1. BASE-PAIR (BP) LEVEL EVALUATION (Genomic Footprint)")
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
        print(" 2. RECIPROCAL OVERLAP FRAGMENT (FRAG) LEVEL EVALUATION (50% Threshold)")
        print("---------------------------------------------------------------------------------")
        print(f"  Total Segtrace SD Pairs:    {len(st_pairs):12,}")
        print(f"  Total SEDEF SD Pairs:       {len(sd_pairs):12,}")
        print(f"  Frag Sensitivity / Recall:  {pair_sn_50*100:8.2f}% ({m_ref_50:,} / {len(sd_pairs):,} SEDEF pairs hit)")
        print(f"  Frag Precision:             {pair_pr_50*100:8.2f}% ({m_tgt_50:,} / {len(st_pairs):,} Segtrace pairs hit)")
        print(f"  Frag F1-Score:              {pair_f1_50*100:8.2f}%")
        print("=================================================================================")

    finally:
        # Clean up temporary working directory unless requested to keep
        if not keep_temp:
            shutil.rmtree(work_dir, ignore_errors=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Segtrace and SEDEF BED files strictly at BP and Reciprocal Overlap Frag levels.")
    parser.add_argument("--segtrace", default="t2t-chm13_sd.dup.bed", help="Path to Segtrace dup.bed file")
    parser.add_argument("--sedef", default="sedef-human-t2tchm13.bed", help="Path to SEDEF bed file")
    parser.add_argument("--out-renamed", default=None, help="Optional output path to save renamed Segtrace BED file")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary intermediate files")
    args = parser.parse_args()

    if not os.path.exists(args.segtrace):
        print(f"[ERROR] Segtrace file '{args.segtrace}' not found.")
        sys.exit(1)
    if not os.path.exists(args.sedef):
        print(f"[ERROR] SEDEF file '{args.sedef}' not found.")
        sys.exit(1)

    run_comparison(args.segtrace, args.sedef, out_renamed=args.out_renamed, keep_temp=args.keep_temp)
