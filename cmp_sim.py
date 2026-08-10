#!/usr/bin/env python3
"""
cmp_sim.py - Simulation benchmark script for Segtrace, SEDEF, and BISER.
Generates simulated genomes with synthetic SDs, runs target tools, and computes BP and Frag metrics.
"""

import random
import subprocess
import os
import time
import sys
import argparse
import numpy as np
import pandas as pd
from cmp_core import evaluate_frag_pairs_fast

def generate_simulated_genome(chrom_sizes, num_dups=100, min_dup_len=1000, max_dup_len=10_000, out_fasta="sim.fa"):
    bases_bytes = np.frombuffer(b'ACGT', dtype=np.uint8)
    genomes = {}
    for chrom, size in chrom_sizes.items():
        genomes[chrom] = bytearray(np.random.choice(bases_bytes, size=size).tobytes())
    
    true_pairs = []
    used_intervals = {chrom: [] for chrom in chrom_sizes}
    chrom_names = list(chrom_sizes.keys())
    
    def is_overlap(chrom, s, e):
        for ts, te in used_intervals[chrom]:
            if max(s, ts) < min(e, te):
                return True
        return False

    print(f"[INFO] Injecting {num_dups} duplications into {out_fasta}...")
    for _ in range(num_dups):
        dup_len = random.randint(min_dup_len, max_dup_len)
        
        while True:
            c1 = random.choice(chrom_names)
            s1 = random.randint(0, chrom_sizes[c1] - dup_len)
            if not is_overlap(c1, s1, s1 + dup_len):
                used_intervals[c1].append((s1, s1 + dup_len))
                break
        
        while True:
            c2 = random.choice(chrom_names)
            s2 = random.randint(0, chrom_sizes[c2] - dup_len)
            if not is_overlap(c2, s2, s2 + dup_len):
                used_intervals[c2].append((s2, s2 + dup_len))
                break
        
        genomes[c2][s2:s2 + dup_len] = genomes[c1][s1:s1 + dup_len]

        div = random.uniform(0.0, 0.1)
        num_muts = np.random.binomial(dup_len, div)
        if num_muts > 0:
            mut_offsets = np.random.choice(dup_len, size=num_muts, replace=False)
            for offset in mut_offsets:
                orig = genomes[c1][s1 + offset]
                if orig == 65: mut = random.choice(b'CGT')
                elif orig == 67: mut = random.choice(b'AGT')
                elif orig == 71: mut = random.choice(b'ACT')
                else: mut = random.choice(b'ACG')
                genomes[c2][s2 + offset] = mut
                
        true_pairs.append(((c1, s1, s1 + dup_len), (c2, s2, s2 + dup_len)))
        
    with open(out_fasta, "wb") as f:
        for chrom, seq in genomes.items():
            f.write(f">{chrom}\n".encode())
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + b"\n")
            
    return true_pairs, out_fasta

def evaluate_sim_bp(true_pairs, pred_pairs):
    if not true_pairs or not pred_pairs:
        return 0.0, 0.0, 0.0

    # Calculate merged BP overlap using temp BED files
    work_dir = "_sim_tmp"
    os.makedirs(work_dir, exist_ok=True)
    t_bed = os.path.join(work_dir, "true.bed")
    p_bed = os.path.join(work_dir, "pred.bed")
    t_m = os.path.join(work_dir, "true_m.bed")
    p_m = os.path.join(work_dir, "pred_m.bed")
    is_bed = os.path.join(work_dir, "is.bed")

    with open(t_bed, 'w') as f:
        for (c1, s1, e1), (c2, s2, e2) in true_pairs:
            f.write(f"{c1}\t{s1}\t{e1}\n{c2}\t{s2}\t{e2}\n")

    with open(p_bed, 'w') as f:
        for (c1, s1, e1), (c2, s2, e2) in pred_pairs:
            f.write(f"{c1}\t{s1}\t{e1}\n{c2}\t{s2}\t{e2}\n")

    subprocess.run(f"bedtools sort -i {t_bed} | bedtools merge -i - > {t_m}", shell=True, check=True)
    subprocess.run(f"bedtools sort -i {p_bed} | bedtools merge -i - > {p_m}", shell=True, check=True)
    subprocess.run(f"bedtools intersect -a {p_m} -b {t_m} > {is_bed}", shell=True, check=True)

    def bp_count(filepath):
        total = 0
        with open(filepath) as f:
            for l in f:
                parts = l.strip().split()
                if len(parts) >= 3:
                    total += int(parts[2]) - int(parts[1])
        return total

    t_bp = bp_count(t_m)
    p_bp = bp_count(p_m)
    i_bp = bp_count(is_bed)

    shutil.rmtree(work_dir, ignore_errors=True)

    rec = i_bp / t_bp if t_bp > 0 else 0.0
    prec = i_bp / p_bp if p_bp > 0 else 0.0
    f1 = 2 * rec * prec / (rec + prec) if (rec + prec) > 0 else 0.0
    return rec, prec, f1

def run_segtrace_sim(fasta_path, true_pairs):
    start_time = time.perf_counter()
    out_prefix = "sim_out"
    subprocess.run(["./segtrace/segtrace", "-k", "15", "-p", "8", fasta_path, "-o", out_prefix],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    exec_time = time.perf_counter() - start_time

    clusters = {}
    dup_file = f"{out_prefix}.dup.bed"
    if os.path.exists(dup_file):
        with open(dup_file) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split()
                if len(parts) < 5: continue
                chrom, start, end, cid, subid = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
                if cid not in clusters:
                    clusters[cid] = []
                clusters[cid].append((chrom, start, end, subid))

    pred_pairs = []
    for cid, regions in clusters.items():
        n = len(regions)
        for i in range(n):
            for j in range(i + 1, n):
                c1, s1, e1, sub1 = regions[i]
                c2, s2, e2, sub2 = regions[j]
                if sub1 == sub2 and sub1 != "0": continue
                if c1 == c2 and max(s1, s2) < min(e1, e2): continue
                pred_pairs.append(((c1, s1, e1), (c2, s2, e2)))

    rec_bp, prec_bp, f1_bp = evaluate_sim_bp(true_pairs, pred_pairs)
    rec_f, prec_f, f1_f, _, _ = evaluate_frag_pairs_fast(true_pairs, pred_pairs, threshold=0.5)
    return rec_bp, prec_bp, f1_bp, rec_f, prec_f, f1_f, exec_time

def main():
    parser = argparse.ArgumentParser(description="Run simulation benchmark comparing SD callers.")
    parser.add_argument("--num-dups", type=int, default=10, help="Number of synthetic SDs to inject")
    parser.add_argument("--genome-size", type=int, default=10_000_000, help="Simulated genome size (bp)")
    args = parser.parse_args()

    chrom_sizes = {'chr1': args.genome_size // 2, 'chr2': args.genome_size // 2}
    true_pairs, fasta_path = generate_simulated_genome(chrom_sizes, num_dups=args.num_dups)

    print(f"\n[BENCHMARK] Evaluating Segtrace on simulated genome ({args.genome_size:,} bp)...")
    r_bp, p_bp, f1_bp, r_f, p_f, f1_f, elapsed = run_segtrace_sim(fasta_path, true_pairs)

    results = [{
        'Tool': 'Segtrace',
        'GenomeSize': args.genome_size,
        'Recall_bp': r_bp,
        'Precision_bp': p_bp,
        'F1-Score_bp': f1_bp,
        'Recall_frag': r_f,
        'Precision_frag': p_f,
        'F1-Score_frag': f1_f,
        'Time(s)': elapsed
    }]

    df = pd.DataFrame(results)
    print("\n=================================================================================")
    print("                      SIMULATION EVALUATION RESULTS")
    print("=================================================================================")
    print(df.to_string(index=False))
    print("=================================================================================")

if __name__ == "__main__":
    main()
