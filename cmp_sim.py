#!/usr/bin/env python3
"""
cmp_sim.py - Simulation benchmark script for Segtrace, SEDEF, and BISER.
Generates synthetic SDs, runs callers, evaluates TP/FP/FN/Recall/Precision/F1, and generates visualization plots.
"""

import random
import subprocess
import os
import time
import argparse
import numpy as np
import pandas as pd

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

        div = random.uniform(0.0, 0.0)
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

def plot_sim_results(results_df, output_png="evaluation_plots.png"):
    """Generates 3-panel visualization plots for simulation benchmark."""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

        # 1. BP-level Metrics
        df_bp = results_df[['Tool', 'Recall_bp', 'Precision_bp']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
        df_bp['Metric'] = df_bp['Metric'].map({'Recall_bp': 'Recall', 'Precision_bp': 'Precision'})
        sns.barplot(data=df_bp, x='Tool', y='Score', hue='Metric', ax=ax1, alpha=0.7)
        ax1.set_ylim(0, 1.05)
        ax1.set_title('BP-level Recall & Precision')
        ax1.grid(axis='y', alpha=0.3)

        # 2. Fragment-level Metrics
        df_frag = results_df[['Tool', 'Recall_frag', 'Precision_frag']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
        df_frag['Metric'] = df_frag['Metric'].map({'Recall_frag': 'Recall', 'Precision_frag': 'Precision'})
        sns.barplot(data=df_frag, x='Tool', y='Score', hue='Metric', ax=ax2, alpha=0.7)
        ax2.set_ylim(0, 1.05)
        ax2.set_title('Fragment-level (bedtools -f 0.5 -r) Recall & Precision')
        ax2.grid(axis='y', alpha=0.3)

        # 3. Time vs F1-Score
        palette = {'Segtrace': 'blue', 'SEDEF': 'red', 'BISER': 'green'}
        sns.scatterplot(data=results_df, x='Time(s)', y='F1_frag', hue='Tool', s=150, palette=palette, ax=ax3)
        ax3.set_xlabel('Execution Time (seconds)')
        ax3.set_ylabel('Fragment F1-Score')
        ax3.set_title('F1-Score vs Execution Time')
        ax3.set_ylim(0, 1.05)
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_png, dpi=300)
        plt.close()
        print(f"[INFO] Visualization saved to {output_png}")
    except Exception as e:
        print(f"[WARNING] Could not generate plots: {e}")

def main():
    parser = argparse.ArgumentParser(description="Run simulation benchmark comparing SD callers.")
    parser.add_argument("--num-dups", type=int, default=None, help="Number of synthetic SDs to inject")
    parser.add_argument("--genome-size", type=int, default=100_000_000, help="Simulated genome size (bp)")
    parser.add_argument("--no-sedef", action='store_true', help="Skip SEDEF benchmark")
    parser.add_argument("--plot", action='store_true', help="Generate evaluation_plots.png visualization")
    args = parser.parse_args()

    g_size = args.genome_size
    print(f"\n=================================================================================")
    print(f"            SIMULATION BENCHMARK REPORT (Genome Size: {g_size:,} bp)")
    print(f"=================================================================================")
    
    chrom_sizes = {f'chr{i}': g_size // 5 for i in range(1, 6)}
    num_dups = args.num_dups if args.num_dups else max(10, g_size // 1_000_000)
    
    true_intervals, fasta_path = sim_generate_genome(chrom_sizes, num_dups=num_dups)
    
    all_results = []
    
    bp_m, frag_m, elapsed = sim_run_segtrace(fasta_path, true_intervals)
    all_results.append({
        'Tool': 'Segtrace',
        'GenomeSize': g_size,
        'Recall_bp': bp_m['recall'],
        'Precision_bp': bp_m['precision'],
        'F1_bp': bp_m['f1'],
        'Recall_frag': frag_m['recall'],
        'Precision_frag': frag_m['precision'],
        'F1_frag': frag_m['f1'],
        'TP': frag_m['tp'],
        'FP': frag_m['fp'],
        'FN': frag_m['fn'],
        'Time(s)': elapsed
    })

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
                all_results.append({
                    'Tool': 'SEDEF',
                    'GenomeSize': g_size,
                    'Recall_bp': s_bp_m['recall'],
                    'Precision_bp': s_bp_m['precision'],
                    'F1_bp': s_bp_m['f1'],
                    'Recall_frag': s_frag_m['recall'],
                    'Precision_frag': s_frag_m['precision'],
                    'F1_frag': s_frag_m['f1'],
                    'TP': s_frag_m['tp'],
                    'FP': s_frag_m['fp'],
                    'FN': s_frag_m['fn'],
                    'Time(s)': sedef_exec_time
                })
                print("\n[SEDEF Performance]")
                print(f"  Time Elapsed:               {sedef_exec_time:.2f} s")
                print(f"  BP Recall / Precision / F1: {s_bp_m['recall']*100:.2f}% / {s_bp_m['precision']*100:.2f}% / {s_bp_m['f1']*100:.2f}%")
                print(f"  FRAG TP / FP / FN:          {s_frag_m['tp']} / {s_frag_m['fp']} / {s_frag_m['fn']}")
                print(f"  FRAG Recall / Prec / F1:    {s_frag_m['recall']*100:.2f}% / {s_frag_m['precision']*100:.2f}% / {s_frag_m['f1']*100:.2f}%")
        except Exception as e:
            print(f"[WARNING] SEDEF execution failed: {e}")
            
    print("=================================================================================")

    df_results = pd.DataFrame(all_results)
    df_results.to_csv("evaluation_results.csv", index=False)
    
    if args.plot or True:  # Generate plots by default
        plot_sim_results(df_results, "evaluation_plots.png")

if __name__ == "__main__":
    main()
