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
import matplotlib.pyplot as plt
import seaborn as sns
from cmp_core import calc_frag_metrics

def sim_generate_genome(chrom_sizes, num_dups=100, min_dup_len=1000, max_dup_len=10_000, out_fasta="sim.fa"):
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
                
    for ext in [".fai", ".sdx"]:
        if os.path.exists(out_fasta + ext):
            os.remove(out_fasta + ext)
            
    return true_pairs, out_fasta

def sim_merge_intervals(pairs):
    by_chrom = {}
    for (c1, s1, e1), (c2, s2, e2) in pairs:
        by_chrom.setdefault(c1, []).append((s1, e1))
        by_chrom.setdefault(c2, []).append((s2, e2))
    
    merged_by_chrom = {}
    for chrom, intervals in by_chrom.items():
        if not intervals:
            continue
        intervals.sort(key=lambda x: x[0])
        merged = [list(intervals[0])]
        for curr in intervals[1:]:
            prev = merged[-1]
            if curr[0] <= prev[1]:
                prev[1] = max(prev[1], curr[1])
            else:
                merged.append(list(curr))
        merged_by_chrom[chrom] = merged
    return merged_by_chrom

def sim_calc_bp_metrics(true_pairs, pred_pairs):
    if not true_pairs or not pred_pairs:
        return 0.0, 0.0, 0.0

    true_m = sim_merge_intervals(true_pairs)
    pred_m = sim_merge_intervals(pred_pairs)

    t_bp = sum(e - s for intervals in true_m.values() for s, e in intervals)
    p_bp = sum(e - s for intervals in pred_m.values() for s, e in intervals)

    i_bp = 0
    common_chroms = set(true_m.keys()).intersection(pred_m.keys())
    for chrom in common_chroms:
        list_t = true_m[chrom]
        list_p = pred_m[chrom]
        i = j = 0
        while i < len(list_t) and j < len(list_p):
            t_s, t_e = list_t[i]
            p_s, p_e = list_p[j]
            overlap = max(0, min(t_e, p_e) - max(t_s, p_s))
            i_bp += overlap
            if t_e < p_e:
                i += 1
            else:
                j += 1

    rec = i_bp / t_bp if t_bp > 0 else 0.0
    prec = i_bp / p_bp if p_bp > 0 else 0.0
    f1 = 2 * rec * prec / (rec + prec) if (rec + prec) > 0 else 0.0
    return rec, prec, f1

def sim_run_segtrace(fasta_path, true_pairs):
    start_time = time.perf_counter()
    out_prefix = "sim_out"
    segtrace_bin = "./segtrace" if os.path.isfile("./segtrace") else "./segtrace/segtrace"
    subprocess.run([segtrace_bin, "-k", "15", "-p", "8", fasta_path, "-o", out_prefix],
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

    rec_bp, prec_bp, f1_bp = sim_calc_bp_metrics(true_pairs, pred_pairs)
    rec_f, prec_f, f1_f, _, _, _ = calc_frag_metrics(true_pairs, pred_pairs, threshold=0.5)
    return rec_bp, prec_bp, f1_bp, rec_f, prec_f, f1_f, exec_time

def evaluate_bedpe_sim(true_pairs, filepath):
    if not os.path.exists(filepath): return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    pred_pairs = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split()
            if len(parts) < 6: continue
            try:
                c1, ra_s, ra_e = parts[0], int(parts[1]), int(parts[2])
                c2, rb_s, rb_e = parts[3], int(parts[4]), int(parts[5])
                if "-" in c1: c1 = c1.split("-", 1)[-1]
                if "-" in c2: c2 = c2.split("-", 1)[-1]
                pred_pairs.append(((c1, ra_s, ra_e), (c2, rb_s, rb_e)))
            except ValueError: continue
            
    rec_bp, prec_bp, f1_bp = sim_calc_bp_metrics(true_pairs, pred_pairs)
    rec_f, prec_f, f1_f, _, _, _ = calc_frag_metrics(true_pairs, pred_pairs, threshold=0.5)
    return rec_bp, prec_bp, f1_bp, rec_f, prec_f, f1_f

def main():
    parser = argparse.ArgumentParser(description="Run simulation benchmark comparing SD callers.")
    parser.add_argument("--num-dups", type=int, default=None, help="Number of synthetic SDs to inject")
    parser.add_argument("--genome-size", type=int, default=None, help="Simulated genome size (bp)")
    parser.add_argument("--no-sedef", action='store_true', help="Skip SEDEF benchmark")
    args = parser.parse_args()

    all_results = []
    N = 5
    if args.genome_size:
        genome_sizes_to_test = [args.genome_size] * N
    else:
        genome_sizes_to_test = [10_000_000] * N + [20_000_000] * N + [50_000_000] * N

    for g_size in genome_sizes_to_test:
        print(f"\n======================================")
        print(f"Testing Genome Size: {g_size:,} bp")
        print(f"======================================")
        
        chrom_sizes = {'chr1': g_size // 2, 'chr2': g_size - (g_size // 2)}
        num_dups = args.num_dups if args.num_dups else max(10, g_size // 1_000_000)
        
        true_pairs, fasta_path = sim_generate_genome(chrom_sizes, num_dups=num_dups)
        
        # Segtrace
        r_bp, p_bp, f1_bp, r_f, p_f, f1_f, elapsed = sim_run_segtrace(fasta_path, true_pairs)
        all_results.append({
            'Tool': 'Segtrace',
            'GenomeSize': g_size,
            'Recall_bp': r_bp,
            'Precision_bp': p_bp,
            'F1-Score_bp': f1_bp,
            'Recall_frag': r_f,
            'Precision_frag': p_f,
            'F1-Score_frag': f1_f,
            'Time(s)': elapsed
        })
        
        if not args.no_sedef:
            # SEDEF
            try:
                env = os.environ.copy()
                sedef_dir = os.path.abspath("sedef")
                env['PATH'] = f"{sedef_dir}:{env.get('PATH', '')}"
                sedef_out_dir = "sedef_out"
                t0 = time.perf_counter()
                subprocess.run([os.path.join(sedef_dir, "sedef.sh"), "-o", sedef_out_dir, "-f", "-j", "8", fasta_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
                sedef_exec_time = time.perf_counter() - t0
                if os.path.exists(f"{sedef_out_dir}/final.bed"):
                    s_r_bp, s_p_bp, s_f1_bp, s_r_f, s_p_f, s_f1_f = evaluate_bedpe_sim(true_pairs, f"{sedef_out_dir}/final.bed")
                    all_results.append({
                        'Tool': 'SEDEF',
                        'GenomeSize': g_size,
                        'Recall_bp': s_r_bp,
                        'Precision_bp': s_p_bp,
                        'F1-Score_bp': s_f1_bp,
                        'Recall_frag': s_r_f,
                        'Precision_frag': s_p_f,
                        'F1-Score_frag': s_f1_f,
                        'Time(s)': sedef_exec_time
                    })
            except Exception as e:
                print(f"SEDEF Failed: {e}")
                
            # BISER
            try:
                env = os.environ.copy()
                env['PATH'] = f"{os.path.expanduser('~/.local/bin')}:{env.get('PATH', '')}"
                biser_out_file = "biser_out.bedpe"
                t0 = time.perf_counter()
                subprocess.run(["biser", "-o", biser_out_file, fasta_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
                biser_exec_time = time.perf_counter() - t0
                if os.path.exists(biser_out_file):
                    b_r_bp, b_p_bp, b_f1_bp, b_r_f, b_p_f, b_f1_f = evaluate_bedpe_sim(true_pairs, biser_out_file)
                    all_results.append({
                        'Tool': 'BISER',
                        'GenomeSize': g_size,
                        'Recall_bp': b_r_bp,
                        'Precision_bp': b_p_bp,
                        'F1-Score_bp': b_f1_bp,
                        'Recall_frag': b_r_f,
                        'Precision_frag': b_p_f,
                        'F1-Score_frag': b_f1_f,
                        'Time(s)': biser_exec_time
                    })
            except Exception as e:
                print(f"BISER Failed: {e}")

    df_all = pd.DataFrame(all_results)
    if not df_all.empty:
        df_all.to_csv("evaluation_results.csv", index=False)
        print("\nSaved evaluation_results.csv")
    else:
        print("No results to save.")
        sys.exit(0)

    # Plotting
    try:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

        df_bp = df_all[['Tool', 'Recall_bp', 'Precision_bp']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
        df_bp['Metric'] = df_bp['Metric'].map({'Recall_bp': 'Recall', 'Precision_bp': 'Precision'})
        sns.barplot(data=df_bp, x='Tool', y='Score', hue='Metric', ax=ax1, alpha=0.6, capsize=.1)
        sns.stripplot(data=df_bp, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax1, palette='dark:black', alpha=0.7, size=5, legend=False)
        ax1.set_ylim(0, 1.1)
        ax1.set_title('BP-level Recall & Precision')
        ax1.grid(axis='y', alpha=0.3)

        df_frag = df_all[['Tool', 'Recall_frag', 'Precision_frag']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
        df_frag['Metric'] = df_frag['Metric'].map({'Recall_frag': 'Recall', 'Precision_frag': 'Precision'})
        sns.barplot(data=df_frag, x='Tool', y='Score', hue='Metric', ax=ax2, alpha=0.6, capsize=.1)
        sns.stripplot(data=df_frag, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax2, palette='dark:black', alpha=0.7, size=5, legend=False)
        ax2.set_ylim(0, 1.1)
        ax2.set_title('Fragment-level Recall & Precision')
        ax2.grid(axis='y', alpha=0.3)

        palette = {'Segtrace': 'blue', 'SEDEF': 'red', 'BISER': 'green'}
        markers = {'Segtrace': 'o', 'SEDEF': '*', 'BISER': 's'}
        sns.scatterplot(data=df_all, x='Time(s)', y='F1-Score_frag', hue='Tool', style='Tool', 
                        palette=palette, markers=markers, s=150, alpha=0.8, ax=ax3)
        for tool in df_all['Tool'].unique():
            tool_data = df_all[df_all['Tool'] == tool].sort_values(by='GenomeSize')
            ax3.plot(tool_data['Time(s)'], tool_data['F1-Score_frag'], color=palette.get(tool, 'gray'), alpha=0.4)

        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('F1-Score (Frag)')
        ax3.set_title('F1-Score vs Execution Time')
        ax3.set_ylim(0, 1.1)
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('evaluation_plots.png', dpi=300)
        plt.close()
        print("Saved evaluation_plots.png")
    except Exception as e:
        print(f"Plotting failed: {e}")

if __name__ == "__main__":
    main()
