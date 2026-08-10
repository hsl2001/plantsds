#!/usr/bin/env python3
"""
cmp_human.py - Comparison CLI & pipeline for real human / pangenome datasets (Segtrace vs SEDEF / BISER).

Includes:
  1) Human RefSeq chromosome accession mappings (NC_060925.1 -> chr1, etc.).
  2) BED file chromosome normalization for Human datasets.
  3) Segtrace, SEDEF, and BISER pair parsing with human chromosome mapping.
  4) Pipeline runner using bedtools (genomic footprint) and bisect (frag pair) evaluation.

Usage:
  python3 cmp_human.py --segtrace results/t2t-chm13_sd.dup.bed --sedef sedef-human-t2tchm13.bed
"""

import sys
import os
import argparse
import subprocess
import shutil

from cmp_core import load_bed_bp, compute_bp_metrics, evaluate_frag_pairs_fast

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
    if '-' in c and not c.startswith('chr'):
        c = c.split('-')[-1]
    return ACC_MAP.get(c, c)

def create_renamed_bed(in_path, out_path):
    """Creates a new BED file with standard chromosome names."""
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

import itertools

def load_segtrace_pairs(in_path):
    """
    Loads paired SD regions from Segtrace .dup.bed file based on cluster_id.
    Filters same subcluster and self-overlapping regions.
    """
    clusters = {}
    if not os.path.exists(in_path):
        return []
    with open(in_path) as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 5:
                c, s, e, cid, subid = parse_chrom(parts[0]), int(parts[1]), int(parts[2]), parts[3], parts[4]
                if cid not in clusters:
                    clusters[cid] = []
                clusters[cid].append((c, s, e, subid))
    
    pairs = []
    for cid, regions in clusters.items():
        if len(regions) < 2:
            continue
        for (ra_c, ra_s, ra_e, ra_sub), (rb_c, rb_s, rb_e, rb_sub) in itertools.combinations(regions, 2):
            if ra_sub == rb_sub and ra_sub != "0":
                continue
            if ra_c == rb_c and ra_s < rb_e and rb_s < ra_e:
                continue
            pairs.append(((ra_c, ra_s, ra_e), (rb_c, rb_s, rb_e)))
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
                if c1 == c2 and max(s1, s2) < min(e1, e2):
                    continue
                pairs.append(((c1, s1, e1), (c2, s2, e2)))
    return pairs

def load_biser_pairs(in_path):
    """Loads paired SD regions from BISER .bed file."""
    pairs = []
    if not os.path.exists(in_path):
        return pairs
    with open(in_path) as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 6:
                c1, s1, e1 = parse_chrom(parts[0]), int(parts[1]), int(parts[2])
                c2, s2, e2 = parse_chrom(parts[3]), int(parts[4]), int(parts[5])
                if c1 == c2 and max(s1, s2) < min(e1, e2):
                    continue
                pairs.append(((c1, s1, e1), (c2, s2, e2)))
    return pairs

def print_report(segtrace_name, sedef_name, st_bp, sd_bp, is_bp, st_u_bp, sd_u_bp,
                 bp_recall, bp_precision, bp_f1, bp_jaccard,
                 st_pairs_count, sd_pairs_count, pair_sn, pair_pr, pair_f1, m_ref, m_tgt):
    """Prints standard human SD comparison report to stdout."""
    print("=================================================================================")
    print("            SEGMENTAL DUPLICATION COMPARISON REPORT (BP & FRAG LEVEL)")
    print("=================================================================================")
    print(f"Segtrace Input:  {segtrace_name}")
    print(f"SEDEF Input:     {sedef_name}")
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
    print(f"  Total Segtrace SD Pairs:    {st_pairs_count:12,}")
    print(f"  Total SEDEF SD Pairs:       {sd_pairs_count:12,}")
    print(f"  Frag Sensitivity / Recall:  {pair_sn*100:8.2f}% ({m_ref:,} / {sd_pairs_count:,} SEDEF pairs hit)")
    print(f"  Frag Precision:             {pair_pr*100:8.2f}% ({m_tgt:,} / {st_pairs_count:,} Segtrace pairs hit)")
    print(f"  Frag F1-Score:              {pair_f1*100:8.2f}%")
    print("=================================================================================")

def run_human_comparison(segtrace_bed, sedef_bed, out_renamed=None, work_dir="_cmp_tmp", keep_temp=False):
    """Human SD comparison pipeline combining bedtools BP evaluation and bisect Frag evaluation."""
    os.makedirs(work_dir, exist_ok=True)

    try:
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

        bp_recall, bp_precision, bp_f1, bp_jaccard = compute_bp_metrics(st_bp, sd_bp, is_bp)

        # 2. Fragment (Frag) Pair-Level Evaluation
        st_pairs = load_segtrace_pairs(segtrace_bed)
        sd_pairs = load_sedef_pairs(sedef_bed)

        pair_sn, pair_pr, pair_f1, m_ref, m_tgt = evaluate_frag_pairs_fast(sd_pairs, st_pairs, threshold=0.5)

        # Print report
        print_report(segtrace_bed, sedef_bed, st_bp, sd_bp, is_bp, st_u_bp, sd_u_bp,
                     bp_recall, bp_precision, bp_f1, bp_jaccard,
                     len(st_pairs), len(sd_pairs), pair_sn, pair_pr, pair_f1, m_ref, m_tgt)

    finally:
        if not keep_temp:
            shutil.rmtree(work_dir, ignore_errors=True)

# Alias for backward compatibility
run_sd_core_comparison = run_human_comparison

def main():
    parser = argparse.ArgumentParser(description="Compare Segtrace and SEDEF BED files on Human/Real genome datasets.")
    parser.add_argument("--segtrace", default="results/t2t-chm13_sd.dup.bed", help="Path to Segtrace dup.bed file")
    parser.add_argument("--sedef", default="sedef-human-t2tchm13.bed", help="Path to SEDEF bed file")
    parser.add_argument("--out-renamed", default=None, help="Optional output path to save renamed Segtrace BED file")
    parser.add_argument("--work-dir", default="_cmp_tmp", help="Temporary working directory")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary intermediate files")
    args = parser.parse_args()

    if not os.path.exists(args.segtrace):
        print(f"[ERROR] Segtrace file '{args.segtrace}' not found.")
        sys.exit(1)
    if not os.path.exists(args.sedef):
        print(f"[ERROR] SEDEF file '{args.sedef}' not found.")
        sys.exit(1)

    run_human_comparison(args.segtrace, args.sedef, out_renamed=args.out_renamed, work_dir=args.work_dir, keep_temp=args.keep_temp)

if __name__ == "__main__":
    main()
