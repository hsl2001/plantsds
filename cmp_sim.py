#!/usr/bin/env python3
"""
cmp_sim.py - Simulation benchmark for SD callers (Segtrace, SEDEF, BISER).

Features:
- Multi-length scaling evaluation across diverse genome sizes.
- Full SD caller comparison (Segtrace, SEDEF, and BISER).
- Base-pair footprint and Fragment-level (50% reciprocal overlap) evaluation.
- Output summary tables, CSV exports, and visualization plots.
"""

import random
import subprocess
import os
import time
import shutil
import argparse
import sys
import numpy as np
import pandas as pd

from cmp_core import (
    parse_bed_intervals,
    calc_bp_metrics,
    eval_fragment_overlap
)

# Register ./sedef to PATH
sedef_local_dir = os.path.abspath("./sedef")
if os.path.isdir(sedef_local_dir):
    os.environ["PATH"] = f"{sedef_local_dir}:{os.environ.get('PATH', '')}"


def sim_generate_genome(chrom_sizes, num_dups=100, min_dup_len=1000, max_dup_len=10_000,
                        ortholog_rate=0.0, flank_len=500, out_fasta="sim.fa"):
    """Generates synthetic chromosomes with segmental duplications."""
    bases_bytes = np.frombuffer(b'ACGT', dtype=np.uint8)
    genomes = {chrom: bytearray(np.random.choice(bases_bytes, size=size).tobytes()) for chrom, size in chrom_sizes.items()}
    true_intervals = []
    used_intervals = {c: [] for c in chrom_sizes}
    chrom_names = list(chrom_sizes.keys())

    def is_overlap(chrom, s, e):
        return any(max(s, ts) < min(e, te) for ts, te in used_intervals[chrom])

    def pick_free(length, extra=0):
        tot = length + 2 * extra
        for _ in range(1000):
            c = random.choice(chrom_names)
            if chrom_sizes[c] <= tot + 10: continue
            s = random.randint(extra, chrom_sizes[c] - length - extra)
            if not is_overlap(c, s - extra, s + length + extra):
                used_intervals[c].append((s - extra, s + length + extra))
                return c, s
        return None, None

    print(f"[INFO] Injecting {num_dups} SDs into {out_fasta}...")
    for _ in range(num_dups):
        dup_len = random.randint(min_dup_len, max_dup_len)
        if ortholog_rate > 0.0 and random.random() < ortholog_rate and len(chrom_names) >= 3:
            c1, s1 = pick_free(dup_len, flank_len)
            c2, s2 = pick_free(dup_len, flank_len)
            c3, s3 = pick_free(dup_len, 0)
            if not c1 or not c2 or not c3: continue
            genomes[c2][s2 - flank_len:s2 + dup_len + flank_len] = genomes[c1][s1 - flank_len:s1 + dup_len + flank_len]
            genomes[c3][s3:s3 + dup_len] = genomes[c1][s1:s1 + dup_len]
            l1, l2, l3 = (c1, s1, s1 + dup_len), (c2, s2, s2 + dup_len), (c3, s3, s3 + dup_len)
            true_intervals.extend([l1, l2, l3])
        else:
            c1, s1 = pick_free(dup_len, 0)
            c2, s2 = pick_free(dup_len, 0)
            if not c1 or not c2: continue
            genomes[c2][s2:s2 + dup_len] = genomes[c1][s1:s1 + dup_len]
            l1, l2 = (c1, s1, s1 + dup_len), (c2, s2, s2 + dup_len)
            true_intervals.extend([l1, l2])

    with open(out_fasta, "wb") as f:
        for c, seq in genomes.items():
            f.write(f">{c}\n".encode())
            for i in range(0, len(seq), 80): f.write(seq[i:i+80] + b"\n")

    for ext in [".fai", ".sdx"]:
        if os.path.exists(out_fasta + ext): os.remove(out_fasta + ext)

    return list(set(true_intervals)), out_fasta

def evaluate_caller_output(bed_path, true_intervals, exec_time):
    """Computes unified BP and fragment metrics for a caller."""
    pred_intervals = parse_bed_intervals(bed_path)
    bp_m = calc_bp_metrics(pred_intervals, true_intervals)
    frag_m = eval_fragment_overlap(pred_intervals, true_intervals, fraction=0.5)

    return {
        'Recall_bp': bp_m['recall'], 'Precision_bp': bp_m['precision'], 'F1_bp': bp_m['f1'],
        'Recall_frag': frag_m['recall'], 'Precision_frag': frag_m['precision'], 'F1_frag': frag_m['f1'],
        'Time(s)': exec_time
    }

def sim_run_segtrace(fasta_path, true_intervals, threads=8, kmer=17,
                     window_size=1024, step_size=0, scale=16):
    """Runs Segtrace with specified parameters and measures Time & Peak RSS Memory (MB)."""
    out_prefix = "sim_out"
    segtrace_bin = "./segtrace" if os.path.isfile("./segtrace") else "./segtrace/segtrace"
    cmd = [
        segtrace_bin,
        "-k", str(kmer),
        "-w", str(window_size),
        "-t", str(step_size),
        "-s", str(scale),
        "-p", str(threads),
        "-o", out_prefix
    ]
    cmd.append(fasta_path)

    if sys.platform == "darwin":
        profile_cmd = ["/usr/bin/time", "-l", *cmd]
    else:
        profile_cmd = ["/usr/bin/time", "-f", "%M", *cmd]

    t0 = time.perf_counter()
    proc = subprocess.run(profile_cmd, check=True, stdout=subprocess.DEVNULL,
                          stderr=subprocess.PIPE, text=True)
    t_elapsed = time.perf_counter() - t0
    try:
        if sys.platform == "darwin":
            mem_bytes = next(int(line.split()[0]) for line in proc.stderr.splitlines()
                             if "maximum resident set size" in line)
            mem_mb = mem_bytes / (1024.0 * 1024.0)
        else:
            mem_mb = int(proc.stderr.strip().splitlines()[-1]) / 1024.0
    except Exception:
        mem_mb = 0.0

    res = evaluate_caller_output(f"{out_prefix}.dup.bed", true_intervals,
                                 t_elapsed)
    res['Tool'] = 'Segtrace'
    res['kmer'] = kmer
    res['window_size'] = window_size
    res['step_size'] = step_size
    res['scale'] = scale
    res['Memory(MB)'] = mem_mb
    return res

def sim_run_sedef(fasta_path, true_intervals, threads=8):
    """Runs SEDEF with robust path handling."""
    sedef_sh = os.path.abspath("sedef/sedef.sh") if os.path.isfile("sedef/sedef.sh") else shutil.which("sedef.sh")
    if not sedef_sh: return None
    try:
        env = os.environ.copy()
        sedef_dir = os.path.dirname(os.path.abspath(sedef_sh))
        env['PATH'] = f"{sedef_dir}:{env.get('PATH', '')}"
        sedef_out_dir = os.path.abspath("sedef_out")
        shutil.rmtree(sedef_out_dir, ignore_errors=True)

        t0 = time.perf_counter()
        res = subprocess.run([sedef_sh, "-o", sedef_out_dir, "-f", "-j", str(threads), os.path.abspath(fasta_path)],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env, text=True)
        t_elapsed = time.perf_counter() - t0

        if res.returncode != 0:
            err_msg = res.stderr.strip() or res.stdout.strip()
            print(f"[WARNING] SEDEF returned non-zero exit status {res.returncode}: {err_msg.splitlines()[-1] if err_msg else 'Unknown error'}")
            return None

        final_bed = os.path.join(sedef_out_dir, "final.bed")
        if os.path.exists(final_bed):
            res_dict = evaluate_caller_output(final_bed, true_intervals, t_elapsed)
            res_dict['Tool'] = 'SEDEF'
            res_dict['Memory(MB)'] = 0.0
            return res_dict
    except Exception as e:
        print(f"[WARNING] SEDEF execution failed: {e}")
    return None

def sim_run_biser(fasta_path, true_intervals, threads=8):
    """Runs BISER if installed."""
    biser_bin = shutil.which("biser") or (os.path.expanduser("~/.local/bin/biser") if os.path.isfile(os.path.expanduser("~/.local/bin/biser")) else None)
    if not biser_bin: return None
    try:
        env = os.environ.copy()
        env['PATH'] = f"{os.path.dirname(os.path.abspath(biser_bin))}:{env.get('PATH', '')}"
        biser_out = "biser_out.bedpe"
        t0 = time.perf_counter()
        cmd = [biser_bin, "-t", str(threads), "-o", biser_out, fasta_path] if threads > 1 else [biser_bin, "-o", biser_out, fasta_path]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, check=True)
        t_elapsed = time.perf_counter() - t0
        if os.path.exists(biser_out):
            res_dict = evaluate_caller_output(biser_out, true_intervals, t_elapsed)
            res_dict['Tool'] = 'BISER'
            res_dict['Memory(MB)'] = 0.0
            return res_dict
    except Exception as e:
        print(f"[WARNING] BISER execution failed: {e}")
    return None

def plot_sim_results(results_df, output_png="evaluation_plots.png"):
    """Generates BP, fragment, and runtime metric plots."""
    if results_df.empty: return
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(19, 5.5))
        palette = {'Segtrace': '#1f77b4', 'SEDEF': '#d62728', 'BISER': '#2ca02c'}
        markers = {'Segtrace': 'o', 'SEDEF': 's', 'BISER': '^'}

        # 1. BP-level metrics
        df_bp = results_df[['Tool', 'Recall_bp', 'Precision_bp']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
        df_bp['Metric'] = df_bp['Metric'].map({'Recall_bp': 'Recall', 'Precision_bp': 'Precision'})
        sns.barplot(data=df_bp, x='Tool', y='Score', hue='Metric', ax=ax1, alpha=0.7, capsize=0.1)
        sns.stripplot(data=df_bp, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax1, palette='dark:black', alpha=0.7, size=5, legend=False)
        ax1.set_ylim(0, 1.08); ax1.set_title('1. Base-Pair Footprint', fontweight='bold'); ax1.grid(axis='y', alpha=0.3)

        # 2. Fragment-level metrics
        df_frag = results_df[['Tool', 'Recall_frag', 'Precision_frag']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
        df_frag['Metric'] = df_frag['Metric'].map({'Recall_frag': 'Recall', 'Precision_frag': 'Precision'})
        sns.barplot(data=df_frag, x='Tool', y='Score', hue='Metric', ax=ax2, alpha=0.7, capsize=0.1)
        sns.stripplot(data=df_frag, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax2, palette='dark:black', alpha=0.7, size=5, legend=False)
        ax2.set_ylim(0, 1.08); ax2.set_title('2. Fragment Level', fontweight='bold'); ax2.grid(axis='y', alpha=0.3)

        # 3. Execution time vs Fragment F1 scaling
        tools = results_df['Tool'].unique()
        sns.scatterplot(data=results_df, x='Time(s)', y='F1_frag', hue='Tool', style='Tool',
                        palette={t: palette.get(t, '#7f7f7f') for t in tools}, markers={t: markers.get(t, 'o') for t in tools}, s=150, alpha=0.85, ax=ax3)
        for t in tools:
            mean_t = results_df[results_df['Tool'] == t].groupby('GenomeSize')[['Time(s)', 'F1_frag']].mean().reset_index()
            ax3.plot(mean_t['Time(s)'], mean_t['F1_frag'], color=palette.get(t, '#7f7f7f'), linestyle='--', alpha=0.6)
        ax3.set_ylim(0, 1.08); ax3.set_title('3. Time vs Fragment F1 Scaling', fontweight='bold'); ax3.grid(True, alpha=0.3)

        plt.tight_layout(); plt.savefig(output_png, dpi=300); plt.close()
        print(f"[INFO] Visualization plot saved to {output_png}")
    except Exception as e:
        print(f"[WARNING] Plotting failed: {e}")

# ==============================================================================
# PARAMETER SWEEP PLOTTING (FRAG F1, MEMORY, TIME)
# ==============================================================================

def plot_sweep_kmer(df, output_png="plot_kmer.png"):
    """Generates 3-panel plot for k-mer size sweep: Frag F1, Memory, Time."""
    if df.empty: return
    try:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
        df_g = df.groupby('kmer')[['Recall_frag', 'Precision_frag', 'F1_frag', 'Memory(MB)', 'Time(s)']].mean().reset_index()

        # 1. Fragment F1
        ax1.plot(df_g['kmer'], df_g['F1_frag'], marker='o', label='Frag F1', color='#1f77b4', linewidth=2.5)
        ax1.plot(df_g['kmer'], df_g['Recall_frag'], marker='s', label='Frag Recall', color='#2ca02c', linewidth=1.5, linestyle='--')
        ax1.plot(df_g['kmer'], df_g['Precision_frag'], marker='^', label='Frag Precision', color='#ff7f0e', linewidth=1.5, linestyle='--')
        ax1.set_xlabel('k-mer Size (-k)'); ax1.set_ylabel('Score'); ax1.set_ylim(0, 1.05); ax1.legend(); ax1.grid(True, alpha=0.3)

        # 2. Peak Memory (MB)
        ax2.plot(df_g['kmer'], df_g['Memory(MB)'], marker='s', color='#d62728', linewidth=2.5)
        ax2.set_xlabel('k-mer Size (-k)'); ax2.set_ylabel('Peak Memory (MB)'); ax2.grid(True, alpha=0.3)

        # 3. Execution Time (s)
        ax3.plot(df_g['kmer'], df_g['Time(s)'], marker='D', color='#9467bd', linewidth=2.5)
        ax3.set_xlabel('k-mer Size (-k)'); ax3.set_ylabel('Time (seconds)'); ax3.grid(True, alpha=0.3)

        plt.tight_layout(); plt.savefig(output_png, dpi=300); plt.close()
        print(f"[INFO] Saved k-mer plot to {output_png}")
    except Exception as e:
        print(f"[WARNING] Plotting k-mer sweep failed: {e}")

def plot_sweep_scale(df, output_png="plot_scale.png"):
    """Generates 3-panel plot for Scale Factor sweep: Frag F1, Memory, Time."""
    if df.empty: return
    try:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
        df_g = df.groupby('scale')[['Recall_frag', 'Precision_frag', 'F1_frag', 'Memory(MB)', 'Time(s)']].mean().reset_index()

        # 1. Fragment F1
        ax1.plot(df_g['scale'], df_g['F1_frag'], marker='o', label='Frag F1', color='#1f77b4', linewidth=2.5)
        ax1.plot(df_g['scale'], df_g['Recall_frag'], marker='s', label='Frag Recall', color='#2ca02c', linewidth=1.5, linestyle='--')
        ax1.plot(df_g['scale'], df_g['Precision_frag'], marker='^', label='Frag Precision', color='#ff7f0e', linewidth=1.5, linestyle='--')
        ax1.set_xlabel('Scale Factor (-s)'); ax1.set_ylabel('Score'); ax1.set_ylim(0, 1.05); ax1.legend(); ax1.grid(True, alpha=0.3)

        # 2. Peak Memory (MB)
        ax2.plot(df_g['scale'], df_g['Memory(MB)'], marker='s', color='#d62728', linewidth=2.5)
        ax2.set_xlabel('Scale Factor (-s)'); ax2.set_ylabel('Peak Memory (MB)'); ax2.grid(True, alpha=0.3)

        # 3. Execution Time (s)
        ax3.plot(df_g['scale'], df_g['Time(s)'], marker='D', color='#9467bd', linewidth=2.5)
        ax3.set_xlabel('Scale Factor (-s)'); ax3.set_ylabel('Time (seconds)'); ax3.grid(True, alpha=0.3)

        plt.tight_layout(); plt.savefig(output_png, dpi=300); plt.close()
        print(f"[INFO] Saved scale plot to {output_png}")
    except Exception as e:
        print(f"[WARNING] Plotting scale sweep failed: {e}")

def plot_sweep_window(df, output_png="plot_window.png"):
    """Generates 3-panel plot for Window Size sweep: Frag F1, Memory, Time."""
    if df.empty: return
    try:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
        df_g = df.groupby('window_size')[['Recall_frag', 'Precision_frag', 'F1_frag', 'Memory(MB)', 'Time(s)']].mean().reset_index()

        # 1. Fragment F1
        ax1.plot(df_g['window_size'], df_g['F1_frag'], marker='o', label='Frag F1', color='#1f77b4', linewidth=2.5)
        ax1.plot(df_g['window_size'], df_g['Recall_frag'], marker='s', label='Frag Recall', color='#2ca02c', linewidth=1.5, linestyle='--')
        ax1.plot(df_g['window_size'], df_g['Precision_frag'], marker='^', label='Frag Precision', color='#ff7f0e', linewidth=1.5, linestyle='--')
        ax1.set_xlabel('Window Size (-w, bp)'); ax1.set_ylabel('Score'); ax1.set_ylim(0, 1.05); ax1.legend(); ax1.grid(True, alpha=0.3)

        # 2. Peak Memory (MB)
        ax2.plot(df_g['window_size'], df_g['Memory(MB)'], marker='s', color='#d62728', linewidth=2.5)
        ax2.set_xlabel('Window Size (-w, bp)'); ax2.set_ylabel('Peak Memory (MB)'); ax2.grid(True, alpha=0.3)

        # 3. Execution Time (s)
        ax3.plot(df_g['window_size'], df_g['Time(s)'], marker='D', color='#9467bd', linewidth=2.5)
        ax3.set_xlabel('Window Size (-w, bp)'); ax3.set_ylabel('Time (seconds)'); ax3.grid(True, alpha=0.3)

        plt.tight_layout(); plt.savefig(output_png, dpi=300); plt.close()
        print(f"[INFO] Saved window plot to {output_png}")
    except Exception as e:
        print(f"[WARNING] Plotting window sweep failed: {e}")

def plot_sweep_step(df, output_png="plot_step.png"):
    """Generates 3-panel plot for Step Size sweep: Frag F1, Memory, Time."""
    if df.empty: return
    try:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
        df_g = df.groupby('step_size')[['Recall_frag', 'Precision_frag', 'F1_frag', 'Memory(MB)', 'Time(s)']].mean().reset_index()

        # 1. Fragment F1
        ax1.plot(df_g['step_size'], df_g['F1_frag'], marker='o', label='Frag F1', color='#1f77b4', linewidth=2.5)
        ax1.plot(df_g['step_size'], df_g['Recall_frag'], marker='s', label='Frag Recall', color='#2ca02c', linewidth=1.5, linestyle='--')
        ax1.plot(df_g['step_size'], df_g['Precision_frag'], marker='^', label='Frag Precision', color='#ff7f0e', linewidth=1.5, linestyle='--')
        ax1.set_xlabel('Step Size (-t, bp)'); ax1.set_ylabel('Score'); ax1.set_ylim(0, 1.05); ax1.legend(); ax1.grid(True, alpha=0.3)

        # 2. Peak Memory (MB)
        ax2.plot(df_g['step_size'], df_g['Memory(MB)'], marker='s', color='#d62728', linewidth=2.5)
        ax2.set_xlabel('Step Size (-t, bp)'); ax2.set_ylabel('Peak Memory (MB)'); ax2.grid(True, alpha=0.3)

        # 3. Execution Time (s)
        ax3.plot(df_g['step_size'], df_g['Time(s)'], marker='D', color='#9467bd', linewidth=2.5)
        ax3.set_xlabel('Step Size (-t, bp)'); ax3.set_ylabel('Time (seconds)'); ax3.grid(True, alpha=0.3)

        plt.tight_layout(); plt.savefig(output_png, dpi=300); plt.close()
        print(f"[INFO] Saved step plot to {output_png}")
    except Exception as e:
        print(f"[WARNING] Plotting step sweep failed: {e}")

# ==============================================================================
# PARAMETER SWEEP RUNNERS
# ==============================================================================

def run_kmer_sweep(genome_sizes, reps, threads, kmers, window_size, step_size, scale, ortholog_rate, num_dups):
    """Executes sweep for k-mer size measuring Time, Memory, and Frag F1."""
    print("\n" + "=" * 85 + "\n    PARAM SWEEP: K-MER SIZE (Time, Memory, Frag F1)\n" + "=" * 85)
    print(f" k-mers: {kmers}\n Window Size: {window_size} | Step Size: {step_size} | Scale: {scale} | Threads: {threads}\n" + "=" * 85 + "\n")

    results = []
    for g_size in genome_sizes:
        chrom_sizes = {f'chr{i}': g_size // 5 for i in range(1, 6)}
        n_dups = num_dups if num_dups else max(10, g_size // 1_000_000)

        for rep in range(1, reps + 1):
            print(f">>> [Genome Size: {g_size:,} bp | Rep {rep}/{reps}] Injecting {n_dups} SDs...")
            true_intervals, fasta_path = sim_generate_genome(chrom_sizes, num_dups=n_dups, ortholog_rate=ortholog_rate, out_fasta="sim.fa")

            for k in kmers:
                res = sim_run_segtrace(fasta_path, true_intervals, threads=threads, kmer=k, window_size=window_size, step_size=step_size, scale=scale)
                res['GenomeSize'], res['Rep'] = g_size, rep
                results.append(res)
                print(f"  [k={k:2d}] Frag F1: {res['F1_frag']*100:6.2f}% | Mem: {res['Memory(MB)']:6.2f}MB | Time: {res['Time(s)']:5.2f}s")

    df = pd.DataFrame(results)
    df.to_csv("sweep_kmer.csv", index=False)
    print(f"\n[INFO] Saved sweep results to sweep_kmer.csv")
    plot_sweep_kmer(df, "plot_kmer.png")
    return df

def run_scale_sweep(genome_sizes, reps, threads, scales, window_size, step_size, kmer, ortholog_rate, num_dups):
    """Executes sweep for scale factor measuring Time, Memory, and Frag F1."""
    print("\n" + "=" * 85 + "\n    PARAM SWEEP: SCALE FACTOR (Time, Memory, Frag F1)\n" + "=" * 85)
    print(f" Scales: {scales}\n Window Size: {window_size} | Step Size: {step_size} | k-mer: {kmer} | Threads: {threads}\n" + "=" * 85 + "\n")

    results = []
    for g_size in genome_sizes:
        chrom_sizes = {f'chr{i}': g_size // 5 for i in range(1, 6)}
        n_dups = num_dups if num_dups else max(10, g_size // 1_000_000)

        for rep in range(1, reps + 1):
            print(f">>> [Genome Size: {g_size:,} bp | Rep {rep}/{reps}] Injecting {n_dups} SDs...")
            true_intervals, fasta_path = sim_generate_genome(chrom_sizes, num_dups=n_dups, ortholog_rate=ortholog_rate, out_fasta="sim.fa")

            for s in scales:
                res = sim_run_segtrace(fasta_path, true_intervals, threads=threads, kmer=kmer, window_size=window_size, step_size=step_size, scale=s)
                res['GenomeSize'], res['Rep'] = g_size, rep
                results.append(res)
                print(f"  [Scale={s:3d}] Frag F1: {res['F1_frag']*100:6.2f}% | Mem: {res['Memory(MB)']:6.2f}MB | Time: {res['Time(s)']:5.2f}s")

    df = pd.DataFrame(results)
    df.to_csv("sweep_scale.csv", index=False)
    print(f"\n[INFO] Saved sweep results to sweep_scale.csv")
    plot_sweep_scale(df, "plot_scale.png")
    return df

def run_window_sweep(genome_sizes, reps, threads, window_sizes, step_size, kmer, scale, ortholog_rate, num_dups):
    """Executes sweep for window size measuring Time, Memory, and Frag F1."""
    print("\n" + "=" * 85 + "\n    PARAM SWEEP: WINDOW SIZE (Time, Memory, Frag F1)\n" + "=" * 85)
    print(f" Window Sizes: {window_sizes}\n Step Size: {step_size} (Auto) | k-mer: {kmer} | Scale: {scale} | Threads: {threads}\n" + "=" * 85 + "\n")

    results = []
    for g_size in genome_sizes:
        chrom_sizes = {f'chr{i}': g_size // 5 for i in range(1, 6)}
        n_dups = num_dups if num_dups else max(10, g_size // 1_000_000)

        for rep in range(1, reps + 1):
            print(f">>> [Genome Size: {g_size:,} bp | Rep {rep}/{reps}] Injecting {n_dups} SDs...")
            true_intervals, fasta_path = sim_generate_genome(chrom_sizes, num_dups=n_dups, ortholog_rate=ortholog_rate, out_fasta="sim.fa")

            for w in window_sizes:
                res = sim_run_segtrace(fasta_path, true_intervals, threads=threads, kmer=kmer, window_size=w, step_size=step_size, scale=scale)
                res['GenomeSize'], res['Rep'] = g_size, rep
                results.append(res)
                print(f"  [Window={w:4d}] Frag F1: {res['F1_frag']*100:6.2f}% | Mem: {res['Memory(MB)']:6.2f}MB | Time: {res['Time(s)']:5.2f}s")

    df = pd.DataFrame(results)
    df.to_csv("sweep_window.csv", index=False)
    print(f"\n[INFO] Saved sweep results to sweep_window.csv")
    plot_sweep_window(df, "plot_window.png")
    return df

def run_step_sweep(genome_sizes, reps, threads, window_size, step_sizes, kmer, scale, ortholog_rate, num_dups):
    """Executes sweep for step size measuring Time, Memory, and Frag F1."""
    print("\n" + "=" * 85 + "\n    PARAM SWEEP: STEP SIZE (Time, Memory, Frag F1)\n" + "=" * 85)
    print(f" Fixed Window: {window_size}\n Step Sizes: {step_sizes}\n k-mer: {kmer} | Scale: {scale} | Threads: {threads}\n" + "=" * 85 + "\n")

    results = []
    for g_size in genome_sizes:
        chrom_sizes = {f'chr{i}': g_size // 5 for i in range(1, 6)}
        n_dups = num_dups if num_dups else max(10, g_size // 1_000_000)

        for rep in range(1, reps + 1):
            print(f">>> [Genome Size: {g_size:,} bp | Rep {rep}/{reps}] Injecting {n_dups} SDs...")
            true_intervals, fasta_path = sim_generate_genome(chrom_sizes, num_dups=n_dups, ortholog_rate=ortholog_rate, out_fasta="sim.fa")

            for t in step_sizes:
                res = sim_run_segtrace(fasta_path, true_intervals, threads=threads, kmer=kmer, window_size=window_size, step_size=t, scale=scale)
                res['GenomeSize'], res['Rep'] = g_size, rep
                results.append(res)
                print(f"  [Step={t:4d}] Frag F1: {res['F1_frag']*100:6.2f}% | Mem: {res['Memory(MB)']:6.2f}MB | Time: {res['Time(s)']:5.2f}s")

    df = pd.DataFrame(results)
    df.to_csv("sweep_step.csv", index=False)
    print(f"\n[INFO] Saved sweep results to sweep_step.csv")
    plot_sweep_step(df, "plot_step.png")
    return df

def main():
    parser = argparse.ArgumentParser(description="Simulation benchmark & parameter sweep for Segtrace.")
    parser.add_argument("--sweep-param", choices=['none', 'kmer', 'scale', 'window', 'step', 'all'], default='none', help="Parameter sweep mode")
    parser.add_argument("--genome-sizes", type=int, nargs='+', default=None, help="Genome sizes to test (e.g. 10000000 20000000)")
    parser.add_argument("--genome-size", type=int, default=None, help="Single genome size to test")
    parser.add_argument("--reps", "-N", type=int, default=3, help="Replicates per genome size (default: 3)")
    parser.add_argument("--num-dups", type=int, default=None, help="Number of synthetic SDs to inject")
    parser.add_argument("--ortholog-rate", type=float, default=0.0, help="Fraction of shared locus orthologs (default: 0.0)")
    parser.add_argument("--threads", "-p", type=int, default=8, help="CPU threads to use (default: 8)")
    
    # Base defaults
    parser.add_argument("--kmer", "-k", type=int, default=17, help="K-mer size for Segtrace (default: 17)")
    parser.add_argument("--window-size", "-w", type=int, default=1024, help="Window size in bp (default: 1024)")
    parser.add_argument("--step-size", "-t", type=int, default=0, help="Step size in bp (default: 0)")
    parser.add_argument("--scale", "-s", type=int, default=16, help="Scale factor (default: 16)")

    # Sweep target ranges
    parser.add_argument("--kmers", type=int, nargs='+', default=[11, 13, 15, 17, 19, 21, 23, 25], help="k-mer sizes for sweep")
    parser.add_argument("--scales", type=int, nargs='+', default=[4, 8, 16, 32, 64, 128, 256], help="Scale factors for sweep")
    parser.add_argument("--window-sizes", type=int, nargs='+', default=[128, 256, 512, 1024, 2048, 4096, 8192], help="Window sizes for sweep")
    parser.add_argument("--step-sizes", type=int, nargs='+', default=[64, 128, 256, 341, 512, 768, 1024], help="Step sizes for sweep (with W=1024)")

    parser.add_argument("--no-sedef", action='store_true', help="Skip SEDEF benchmark")
    parser.add_argument("--no-biser", action='store_true', help="Skip BISER benchmark")
    parser.add_argument("--out-csv", default="evaluation_results.csv", help="Output results CSV")
    parser.add_argument("--out-plot", default="evaluation_plots.png", help="Output plot PNG")
    parser.add_argument("--no-plot", action='store_true', help="Disable plotting")
    args = parser.parse_args()

    genome_sizes = args.genome_sizes if args.genome_sizes else ([args.genome_size] if args.genome_size else [10_000_000, 20_000_000])
    reps = max(1, args.reps)

    if args.sweep_param != 'none':
        if args.sweep_param in ['kmer', 'all']:
            run_kmer_sweep(genome_sizes, reps, args.threads, args.kmers, args.window_size, args.step_size, args.scale, args.ortholog_rate, args.num_dups)
        if args.sweep_param in ['scale', 'all']:
            run_scale_sweep(genome_sizes, reps, args.threads, args.scales, args.window_size, args.step_size, args.kmer, args.ortholog_rate, args.num_dups)
        if args.sweep_param in ['window', 'all']:
            run_window_sweep(genome_sizes, reps, args.threads, args.window_sizes, args.step_size, args.kmer, args.scale, args.ortholog_rate, args.num_dups)
        if args.sweep_param in ['step', 'all']:
            run_step_sweep(genome_sizes, reps, args.threads, args.window_size, args.step_sizes, args.kmer, args.scale, args.ortholog_rate, args.num_dups)
        return

    # Standard comparative benchmark if sweep_param is 'none'
    biser_avail = not args.no_biser and (shutil.which("biser") is not None or os.path.isfile(os.path.expanduser("~/.local/bin/biser")))
    sedef_avail = not args.no_sedef and (os.path.isfile("sedef/sedef.sh") or shutil.which("sedef.sh") is not None)
    if not biser_avail and not args.no_biser: print("[INFO] BISER binary not detected in PATH. Skipping BISER.")
    if not sedef_avail and not args.no_sedef: print("[INFO] SEDEF binary not detected in 'sedef/sedef.sh'. Skipping SEDEF.")

    print("\n" + "=" * 85 + "\n        SEGMENTAL DUPLICATION CALLER SIMULATION BENCHMARK\n" + "=" * 85)
    print(f" Genome Sizes:   {', '.join(f'{g:,} bp' for g in genome_sizes)}\n Replicates:     {reps} per size\n CPU Threads:    {args.threads}\n" + "=" * 85 + "\n")

    all_results = []
    for g_size in genome_sizes:
        chrom_sizes = {f'chr{i}': g_size // 5 for i in range(1, 6)}
        n_dups = args.num_dups if args.num_dups else max(10, g_size // 1_000_000)

        for rep in range(1, reps + 1):
            print(f">>> [Genome Size: {g_size:,} bp | Rep {rep}/{reps} | SDs: {n_dups}] <<<")
            true_intervals, fasta_path = sim_generate_genome(chrom_sizes, num_dups=n_dups, ortholog_rate=args.ortholog_rate, out_fasta="sim.fa")

            # 1. Segtrace
            st_res = sim_run_segtrace(fasta_path, true_intervals, threads=args.threads, kmer=args.kmer, window_size=args.window_size, step_size=args.step_size, scale=args.scale)
            st_res['GenomeSize'], st_res['Rep'] = g_size, rep
            all_results.append(st_res)
            print(f"  [Segtrace] Time: {st_res['Time(s)']:.2f}s | BP F1: {st_res['F1_bp']*100:.2f}% | Frag F1: {st_res['F1_frag']*100:.2f}%")

            # 2. SEDEF
            if sedef_avail:
                sd_res = sim_run_sedef(fasta_path, true_intervals, threads=args.threads)
                if sd_res:
                    sd_res['GenomeSize'], sd_res['Rep'] = g_size, rep
                    all_results.append(sd_res)
                    print(f"  [SEDEF]    Time: {sd_res['Time(s)']:.2f}s | BP F1: {sd_res['F1_bp']*100:.2f}% | Frag F1: {sd_res['F1_frag']*100:.2f}%")

            # 3. BISER
            if biser_avail:
                bi_res = sim_run_biser(fasta_path, true_intervals, threads=args.threads)
                if bi_res:
                    bi_res['GenomeSize'], bi_res['Rep'] = g_size, rep
                    all_results.append(bi_res)
                    print(f"  [BISER]    Time: {bi_res['Time(s)']:.2f}s | BP F1: {bi_res['F1_bp']*100:.2f}% | Frag F1: {bi_res['F1_frag']*100:.2f}%")
            print("-" * 85)

    df_results = pd.DataFrame(all_results)
    if not df_results.empty:
        df_results.to_csv(args.out_csv, index=False)
        print(f"\n[INFO] Complete results saved to {args.out_csv}")

        summary_cols = ['Tool', 'GenomeSize', 'Recall_bp', 'Precision_bp', 'F1_bp', 'Recall_frag', 'Precision_frag', 'F1_frag', 'Time(s)']
        df_summary = df_results[summary_cols].groupby(['Tool', 'GenomeSize']).mean().reset_index()

        print("\n" + "=" * 115 + "\n                                     AVERAGE BENCHMARK SUMMARY\n" + "=" * 115)
        print(f"{'Tool':<12} {'Size':<12} {'BP F1':<9} {'Frag Rec':<10} {'Frag Prec':<10} {'Frag F1':<10} {'Time(s)':<8}")
        print("-" * 115)
        for _, row in df_summary.iterrows():
            print(f"{row['Tool']:<12} {int(row['GenomeSize']):<12,} {row['F1_bp']*100:>6.2f}%   {row['Recall_frag']*100:>6.2f}%   {row['Precision_frag']*100:>6.2f}%   {row['F1_frag']*100:>6.2f}%   {row['Time(s)']:>6.2f}s")
        print("=" * 115 + "\n")

        if not args.no_plot:
            plot_sim_results(df_results, args.out_plot)

if __name__ == "__main__":
    main()

