#!/usr/bin/env python3
"""
cmp_human.py - Comparison CLI for real human / pangenome datasets (Segtrace vs SEDEF / CHM13).
Evaluates Base-Pair (BP) footprint and Fragment-level Reciprocal 50% Overlap (bedtools -f 0.5 -r: TP, FP, FN).
"""

import sys
import os
import argparse
from cmp_core import parse_bed_intervals, calc_bp_metrics, eval_reciprocal_overlap

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

def load_normalized_intervals(path, exclude_chrm=True):
    raw_ints = parse_bed_intervals(path)
    norm_ints = []
    for c, s, e in raw_ints:
        c_norm = parse_chrom(c)
        if not (exclude_chrm and c_norm == 'chrM'):
            norm_ints.append((c_norm, s, e))
    return list(set(norm_ints))

def print_report(segtrace_name, sedef_name, bp_m, frag_m):
    """Prints standard human SD comparison report with TP, FP, FN."""
    print("=================================================================================")
    print("            SEGMENTAL DUPLICATION COMPARISON REPORT")
    print("=================================================================================")
    print(f"Segtrace Input:  {segtrace_name}")
    print(f"SEDEF Input:     {sedef_name}")
    print("---------------------------------------------------------------------------------")
    print(" BASE-PAIR (BP) LEVEL EVALUATION (Genomic Footprint)")
    print("---------------------------------------------------------------------------------")
    print(f"  Segtrace Merged Footprint:  {bp_m['pred_bp']:12,} bp ({bp_m['pred_bp']/1e6:8.2f} Mb)")
    print(f"  SEDEF Merged Footprint:     {bp_m['ref_bp']:12,} bp ({bp_m['ref_bp']/1e6:8.2f} Mb)")
    print(f"  Overlap (Intersection) BP:  {bp_m['is_bp']:12,} bp ({bp_m['is_bp']/1e6:8.2f} Mb)")
    print(f"  Segtrace Unique Footprint:  {bp_m['pred_unique_bp']:12,} bp ({bp_m['pred_unique_bp']/1e6:8.2f} Mb)")
    print(f"  SEDEF Unique Footprint:     {bp_m['ref_unique_bp']:12,} bp ({bp_m['ref_unique_bp']/1e6:8.2f} Mb)")
    print("  -------------------------------------------------------------------------------")
    print(f"  BP Sensitivity / Recall:    {bp_m['recall']*100:8.2f}%")
    print(f"  BP Precision:               {bp_m['precision']*100:8.2f}%")
    print(f"  BP F1-Score:                {bp_m['f1']*100:8.2f}%")
    print(f"  BP Jaccard Index:           {bp_m['jaccard']:10.6f}")
    print("---------------------------------------------------------------------------------")
    print(" FRAGMENT (FRAG) LEVEL EVALUATION")
    print("---------------------------------------------------------------------------------")
    print(f"  True Positives (TP):        {frag_m['tp']:12,}")
    print(f"  False Positives (FP):       {frag_m['fp']:12,}")
    print(f"  False Negatives (FN):       {frag_m['fn']:12,}")
    print(f"  Total Reference Fragments:  {frag_m['total_ref']:12,}")
    print(f"  Total Predicted Fragments:  {frag_m['total_pred']:12,}")
    print("  -------------------------------------------------------------------------------")
    print(f"  FRAG Sensitivity / Recall:  {frag_m['recall']*100:8.2f}%")
    print(f"  FRAG Precision:             {frag_m['precision']*100:8.2f}%")
    print(f"  FRAG F1-Score:              {frag_m['f1']*100:8.2f}%")
    print("=================================================================================")

def cmp_human_run(segtrace_bed, sedef_bed, exclude_chrm=True):
    st_ints = load_normalized_intervals(segtrace_bed, exclude_chrm=exclude_chrm)
    sd_ints = load_normalized_intervals(sedef_bed, exclude_chrm=exclude_chrm)

    bp_m = calc_bp_metrics(st_ints, sd_ints)
    frag_m = eval_reciprocal_overlap(st_ints, sd_ints)

    print_report(segtrace_bed, sedef_bed, bp_m, frag_m)
    return bp_m, frag_m

def main():
    parser = argparse.ArgumentParser(description="Compare Segtrace and SEDEF BED files on Human/Real genome datasets.")
    parser.add_argument("--segtrace", default="t2t-chm13_sd.dup.bed", help="Path to Segtrace dup.bed file")
    parser.add_argument("--sedef", default="data/chm13v2.0_SD.bed", help="Path to SEDEF/CHM13 bed file")
    parser.add_argument("--include-chrm", action="store_true", help="Include mitochondrial chromosome (chrM) in evaluation")
    args = parser.parse_args()

    if not os.path.exists(args.segtrace):
        print(f"[ERROR] Segtrace file '{args.segtrace}' not found.")
        sys.exit(1)
    if not os.path.exists(args.sedef):
        print(f"[ERROR] SEDEF file '{args.sedef}' not found.")
        sys.exit(1)

    cmp_human_run(args.segtrace, args.sedef, exclude_chrm=not args.include_chrm)

if __name__ == "__main__":
    main()
