#!/usr/bin/env python3
"""
cmp_human.py - Comparison CLI & pipeline for real human / pangenome datasets (Segtrace vs SEDEF / CHM13).
"""

import sys
import os

import argparse
import subprocess
import shutil
import gzip

from cmp_core import count_bed_bp, calc_bp_metrics, calc_frag_metrics

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

def chrom_key(c):
    name = c[3:] if c.startswith('chr') else c
    if name.isdigit():
        return (0, int(name), c)
    return (1, name, c)

def open_file(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def load_fragments(filepath):
    frags = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                frags.append([(parts[0], int(parts[1]), int(parts[2]))])
    return frags

def cmp_human_normalize_bed(in_path, out_path):
    """Normalizes chromosome names for Segtrace and SEDEF/CHM13 BED files."""
    count = 0
    with open_file(in_path) as fin, open(out_path, 'w') as fout:
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
            elif len(parts) >= 4 and ':' in parts[3] and '-' in parts[3]:
                c2_str, pos_str = parts[3].split(':', 1)
                s2_str, e2_str = pos_str.split('-', 1)
                fout.write(f"{parse_chrom(c2_str)}\t{s2_str}\t{e2_str}\n")
                count += 1
    return count

def print_report(segtrace_name, sedef_name, st_bp, sd_bp, is_bp, st_u_bp, sd_u_bp,
                 bp_recall, bp_precision, bp_f1, bp_jaccard,
                 frag_recall, frag_precision, frag_f1,
                 per_chrom=None):
    """Prints standard human SD comparison report to stdout."""
    print("=================================================================================")
    print("            SEGMENTAL DUPLICATION COMPARISON REPORT")
    print("=================================================================================")
    print(f"Segtrace Input:  {segtrace_name}")
    print(f"SEDEF Input:     {sedef_name}")
    print("---------------------------------------------------------------------------------")
    print(" BASE-PAIR (BP) LEVEL EVALUATION (Genomic Footprint)")
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
    print("---------------------------------------------------------------------------------")
    print(" FRAGMENT (FRAG) LEVEL EVALUATION (Reciprocal Overlap)")
    print("---------------------------------------------------------------------------------")
    print(f"  FRAG Sensitivity / Recall:  {frag_recall*100:8.2f}%")
    print(f"  FRAG Precision:             {frag_precision*100:8.2f}%")
    print(f"  FRAG F1-Score:              {frag_f1*100:8.2f}%")
    if per_chrom:
        print("---------------------------------------------------------------------------------")
        print(" PER-CHROMOSOME ACCURACY EVALUATION")
        print("---------------------------------------------------------------------------------")
        print(f" {'Chr':<7} | {'Segtrace(BP)':>12} | {'SEDEF(BP)':>12} | {'BP Rec.':>8} | {'BP Prec.':>8} | {'BP F1':>7} | {'Frag Rec.':>9} | {'Frag Prec.':>10} | {'Frag F1':>7}")
        print("---------------------------------------------------------------------------------")
        for r in per_chrom:
            print(f" {r['chrom']:<7} | {r['st_bp']:12,} | {r['sd_bp']:12,} | "
                  f"{r['bp_recall']*100:7.2f}% | {r['bp_precision']*100:7.2f}% | {r['bp_f1']*100:6.2f}% | "
                  f"{r['frag_recall']*100:8.2f}% | {r['frag_precision']*100:9.2f}% | {r['frag_f1']*100:6.2f}%")
    print("=================================================================================")

def get_merged(path):
    by_c = {}
    with open(path) as f:
        for l in f:
            p = l.strip().split()
            if len(p) >= 3:
                by_c.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    merged = {}
    for c, ints in by_c.items():
        ints.sort()
        m = [list(ints[0])]
        for curr in ints[1:]:
            if curr[0] <= m[-1][1]: m[-1][1] = max(m[-1][1], curr[1])
            else: m.append(list(curr))
        merged[c] = m
    return merged

def cmp_human_calc_bp(st_norm, sd_norm, work_dir, use_bedtools=False):
    """Calculates BP footprint overlap using fast in-memory Python calculation or optional bedtools."""
    if use_bedtools and shutil.which("bedtools") is not None:
        st_sort, st_merged = os.path.join(work_dir, "st_sort.bed"), os.path.join(work_dir, "st_merged.bed")
        sd_sort, sd_merged = os.path.join(work_dir, "sd_sort.bed"), os.path.join(work_dir, "sd_merged.bed")
        is_file, st_u_file, sd_u_file = [os.path.join(work_dir, f) for f in ("is.bed", "st_u.bed", "sd_u.bed")]

        subprocess.run(f"bedtools sort -i {st_norm} | bedtools merge -i - > {st_merged}", shell=True, check=True)
        subprocess.run(f"bedtools sort -i {sd_norm} | bedtools merge -i - > {sd_merged}", shell=True, check=True)
        subprocess.run(f"bedtools intersect -a {st_merged} -b {sd_merged} > {is_file}", shell=True, check=True)
        subprocess.run(f"bedtools subtract -a {st_merged} -b {sd_merged} > {st_u_file}", shell=True, check=True)
        subprocess.run(f"bedtools subtract -a {sd_merged} -b {st_merged} > {sd_u_file}", shell=True, check=True)

        return count_bed_bp(st_merged), count_bed_bp(sd_merged), count_bed_bp(is_file), count_bed_bp(st_u_file), count_bed_bp(sd_u_file)

    st_m, sd_m = get_merged(st_norm), get_merged(sd_norm)
    st_bp = sum(e - s for ints in st_m.values() for s, e in ints)
    sd_bp = sum(e - s for ints in sd_m.values() for s, e in ints)
    is_bp = 0
    for c in set(st_m.keys()).intersection(sd_m.keys()):
        l1, l2 = st_m[c], sd_m[c]
        i = j = 0
        while i < len(l1) and j < len(l2):
            overlap = max(0, min(l1[i][1], l2[j][1]) - max(l1[i][0], l2[j][0]))
            is_bp += overlap
            if l1[i][1] < l2[j][1]: i += 1
            else: j += 1
    return st_bp, sd_bp, is_bp, st_bp - is_bp, sd_bp - is_bp

def cmp_human_calc_per_chrom(st_norm, sd_norm, st_frags, sd_frags):
    st_m, sd_m = get_merged(st_norm), get_merged(sd_norm)
    all_chroms = sorted(set(st_m.keys()).union(sd_m.keys()), key=chrom_key)
    
    per_chrom = []
    for c in all_chroms:
        st_ints = st_m.get(c, [])
        sd_ints = sd_m.get(c, [])
        st_c_bp = sum(e - s for s, e in st_ints)
        sd_c_bp = sum(e - s for s, e in sd_ints)
        
        is_c_bp = 0
        if st_ints and sd_ints:
            i = j = 0
            while i < len(st_ints) and j < len(sd_ints):
                overlap = max(0, min(st_ints[i][1], sd_ints[j][1]) - max(st_ints[i][0], sd_ints[j][0]))
                is_c_bp += overlap
                if st_ints[i][1] < sd_ints[j][1]: i += 1
                else: j += 1
                
        bp_rec, bp_prec, bp_f1, bp_jacc = calc_bp_metrics(st_c_bp, sd_c_bp, is_c_bp)
        
        st_frags_c = [p for p in st_frags if p[0][0] == c]
        sd_frags_c = [p for p in sd_frags if p[0][0] == c]
        frag_rec, frag_prec, frag_f1, _, _, _ = calc_frag_metrics(sd_frags_c, st_frags_c)
        
        per_chrom.append({
            'chrom': c,
            'st_bp': st_c_bp,
            'sd_bp': sd_c_bp,
            'is_bp': is_c_bp,
            'bp_recall': bp_rec,
            'bp_precision': bp_prec,
            'bp_f1': bp_f1,
            'bp_jaccard': bp_jacc,
            'frag_recall': frag_rec,
            'frag_precision': frag_prec,
            'frag_f1': frag_f1,
        })
    return per_chrom

def cmp_human_run(segtrace_bed, sedef_bed, work_dir="_cmp_tmp", keep_temp=False, show_per_chrom=True):
    """Human SD comparison pipeline combining BP evaluation and footprint Coverage Frag evaluation."""
    os.makedirs(work_dir, exist_ok=True)
    try:
        st_norm, sd_norm = os.path.join(work_dir, "st_norm.bed"), os.path.join(work_dir, "sd_norm.bed")
        cmp_human_normalize_bed(segtrace_bed, st_norm)
        cmp_human_normalize_bed(sedef_bed, sd_norm)

        st_bp, sd_bp, is_bp, st_u_bp, sd_u_bp = cmp_human_calc_bp(st_norm, sd_norm, work_dir)
        bp_recall, bp_precision, bp_f1, bp_jaccard = calc_bp_metrics(st_bp, sd_bp, is_bp)
        
        st_frags = load_fragments(st_norm)
        sd_frags = load_fragments(sd_norm)
        
        frag_recall, frag_precision, frag_f1, _, _, _ = calc_frag_metrics(sd_frags, st_frags)

        per_chrom = cmp_human_calc_per_chrom(st_norm, sd_norm, st_frags, sd_frags) if show_per_chrom else None

        print_report(segtrace_bed, sedef_bed, st_bp, sd_bp, is_bp, st_u_bp, sd_u_bp,
                     bp_recall, bp_precision, bp_f1, bp_jaccard,
                     frag_recall, frag_precision, frag_f1,
                     per_chrom=per_chrom)
    finally:
        if not keep_temp:
            shutil.rmtree(work_dir, ignore_errors=True)

# Alias for backward compatibility
run_sd_core_comparison = cmp_human_run

def main():
    parser = argparse.ArgumentParser(description="Compare Segtrace and SEDEF BED files on Human/Real genome datasets.")
    parser.add_argument("--segtrace", default="t2t-chm13_sd.dup.bed", help="Path to Segtrace dup.bed file")
    parser.add_argument("--sedef", default="data/chm13v2.0_SD.bed", help="Path to SEDEF/CHM13 bed file")
    parser.add_argument("--work-dir", default="_cmp_tmp", help="Temporary working directory")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary intermediate files")
    parser.add_argument("--no-per-chrom", action="store_true", help="Hide per-chromosome accuracy table")
    args = parser.parse_args()

    if not os.path.exists(args.segtrace):
        print(f"[ERROR] Segtrace file '{args.segtrace}' not found.")
        sys.exit(1)
    if not os.path.exists(args.sedef):
        print(f"[ERROR] SEDEF file '{args.sedef}' not found.")
        sys.exit(1)

    cmp_human_run(args.segtrace, args.sedef, work_dir=args.work_dir, keep_temp=args.keep_temp, show_per_chrom=not args.no_per_chrom)

if __name__ == "__main__":
    main()
