#!/usr/bin/env python3
"""
cmp_sim.py - Simulation benchmark script for Segtrace, SEDEF, and BISER.
"""

import random
import subprocess
import os
import time
import argparse
import numpy as np
from cmp_core import parse_bed_intervals, calc_bp_metrics, eval_reciprocal_overlap

def sim_generate_genome(chrom_sizes, num_dups=100, min_dup_len=1000, max_dup_len=10_000, out_fasta="sim.fa"):
    bases_bytes = np.frombuffer(b'ACGT', dtype=np.uint8)
    genomes = {}
    for chrom, size in chrom_sizes.items():
        genomes[chrom] = bytearray(np.random.choice(bases_bytes, size=size).tobytes())
    
    true_intervals = []
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
                
        true_intervals.extend([(c1, s1, s1 + dup_len), (c2, s2, s2 + dup_len)])
        
    with open(out_fasta, "wb") as f:
        for chrom, seq in genomes.items():
            f.write(f">{chrom}\n".encode())
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + b"\n")
                
    for ext in [".fai", ".sdx"]:
        if os.path.exists(out_fasta + ext):
            os.remove(out_fasta + ext)
            
    return list(set(true_intervals)), out_fasta

def sim_run_segtrace(fasta_path, true_intervals):
    start_time = time.perf_counter()
    out_prefix = "sim_out"
    segtrace_bin = "./segtrace" if os.path.isfile("./segtrace") else "./segtrace/segtrace"
    subprocess.run([segtrace_bin, "-k", "15", "-p", "8", fasta_path, "-o", out_prefix],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    exec_time = time.perf_counter() - start_time

    dup_file = f"{out_prefix}.dup.bed"
    pred_intervals = parse_bed_intervals(dup_file)

    bp_m = calc_bp_metrics(pred_intervals, true_intervals)
    frag_m = eval_reciprocal_overlap(pred_intervals, true_intervals, fraction=0.5)

    return bp_m, frag_m, exec_time

def evaluate_bedpe_sim(true_intervals, filepath):
    pred_intervals = parse_bed_intervals(filepath)
    bp_m = calc_bp_metrics(pred_intervals, true_intervals)
    frag_m = eval_reciprocal_overlap(pred_intervals, true_intervals, fraction=0.5)
    return bp_m, frag_m

def main():
    parser = argparse.ArgumentParser(description="Run simulation benchmark comparing SD callers.")
    parser.add_argument("--num-dups", type=int, default=None, help="Number of synthetic SDs to inject")
    parser.add_argument("--genome-size", type=int, default=10_000_000, help="Simulated genome size (bp)")
    parser.add_argument("--no-sedef", action='store_true', help="Skip SEDEF benchmark")
    args = parser.parse_args()

    g_size = args.genome_size
    print(f"\n=================================================================================")
    print(f"            SIMULATION BENCHMARK REPORT (Genome Size: {g_size:,} bp)")
    print(f"=================================================================================")
    
    chrom_sizes = {'chr1': g_size // 2, 'chr2': g_size - (g_size // 2)}
    num_dups = args.num_dups if args.num_dups else max(10, g_size // 1_000_000)
    
    true_intervals, fasta_path = sim_generate_genome(chrom_sizes, num_dups=num_dups)
    
    bp_m, frag_m, elapsed = sim_run_segtrace(fasta_path, true_intervals)
    
    print("\n[Segtrace Performance]")
    print(f"  Time Elapsed:               {elapsed:.2f} s")
    print(f"  BP Recall / Precision / F1: {bp_m['recall']*100:.2f}% / {bp_m['precision']*100:.2f}% / {bp_m['f1']*100:.2f}%")
    print(f"  FRAG TP / FP / FN:          {frag_m['tp']} / {frag_m['fp']} / {frag_m['fn']}")
    print(f"  FRAG Recall / Prec / F1:    {frag_m['recall']*100:.2f}% / {frag_m['precision']*100:.2f}% / {frag_m['f1']*100:.2f}%")

    if not args.no_sedef:
        try:
            env = os.environ.copy()
            sedef_dir = os.path.abspath("sedef")
            env['PATH'] = f"{sedef_dir}:{env.get('PATH', '')}"
            sedef_out_dir = "sedef_out"
            t0 = time.perf_counter()
            subprocess.run([os.path.join(sedef_dir, "sedef.sh"), "-o", sedef_out_dir, "-f", "-j", "8", fasta_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
            sedef_exec_time = time.perf_counter() - t0
            if os.path.exists(f"{sedef_out_dir}/final.bed"):
                s_bp_m, s_frag_m = evaluate_bedpe_sim(true_intervals, f"{sedef_out_dir}/final.bed")
                print("\n[SEDEF Performance]")
                print(f"  Time Elapsed:               {sedef_exec_time:.2f} s")
                print(f"  BP Recall / Precision / F1: {s_bp_m['recall']*100:.2f}% / {s_bp_m['precision']*100:.2f}% / {s_bp_m['f1']*100:.2f}%")
                print(f"  FRAG TP / FP / FN:          {s_frag_m['tp']} / {s_frag_m['fp']} / {s_frag_m['fn']}")
                print(f"  FRAG Recall / Prec / F1:    {s_frag_m['recall']*100:.2f}% / {s_frag_m['precision']*100:.2f}% / {s_frag_m['f1']*100:.2f}%")
        except Exception as e:
            print(f"[WARNING] SEDEF execution failed: {e}")
            
    print("=================================================================================")

if __name__ == "__main__":
    main()
