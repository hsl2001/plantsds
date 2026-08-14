#!/usr/bin/env python3
"""
cmp_sim.py - Comprehensive Simulation Benchmark for SD Callers (Segtrace, SEDEF, BISER).

Features:
1. Subclustering Comparison (Default):
   - Evaluates Segtrace both with subclustering consideration (Segtrace) and without subclustering
     (Segtrace (No Subcluster)) alongside SEDEF and BISER.
2. Multi-Length Scaling Benchmark:
   - Evaluates performance across diverse genome sizes (e.g. 10Mb, 20Mb, 50Mb, 100Mb) with replicate runs.
3. SEDEF & BISER Benchmark Integration:
   - Runs SEDEF and BISER with graceful fallback if binaries are not present.
4. Comprehensive Metrics & Visualizations:
   - Evaluates Base-Pair (BP) level footprint, Fragment-level (50% reciprocal overlap), and Pair-level metrics.
   - Generates publication-ready 3-panel visualization plots and results CSV.
"""

import random
import subprocess
import os
import sys
import time
import shutil
import argparse
import numpy as np
import pandas as pd

from cmp_core import (
    parse_bed_intervals,
    calc_bp_metrics,
    eval_fragment_overlap,
    load_segtrace_pairs,
    load_bedpe_pairs,
    evaluate_frag_pairs_fast
)

def sim_generate_genome(chrom_sizes, num_dups=100, min_dup_len=1000, max_dup_len=10_000,
                        divergence=0.0, ortholog_rate=0.0, flank_len=500, out_fasta="sim.fa"):
    """
    Generates synthetic chromosomes and injects segmental duplications.
    - Standard Paralogs (rate 1 - ortholog_rate): Duplications to distinct genomic loci with different flanking sequences.
    - Pangenome Shared Locus / Orthologs (rate ortholog_rate): Duplications where the same insertion site
      (with matching flanking sequences) is shared across chromosomes/genomes, plus an additional novel insertion.
    Returns: true_pairs, true_intervals, out_fasta
    """
    bases_bytes = np.frombuffer(b'ACGT', dtype=np.uint8)
    genomes = {}
    for chrom, size in chrom_sizes.items():
        genomes[chrom] = bytearray(np.random.choice(bases_bytes, size=size).tobytes())
    
    true_pairs = []
    true_intervals = []
    used_intervals = {chrom: [] for chrom in chrom_sizes}
    chrom_names = list(chrom_sizes.keys())
    
    def is_overlap(chrom, s, e):
        for ts, te in used_intervals[chrom]:
            if max(s, ts) < min(e, te):
                return True
        return False

    def pick_free_interval(length, extra_flank=0):
        total_len = length + 2 * extra_flank
        for _ in range(1000):
            c = random.choice(chrom_names)
            if chrom_sizes[c] <= total_len + 10:
                continue
            s = random.randint(extra_flank, chrom_sizes[c] - length - extra_flank)
            if not is_overlap(c, s - extra_flank, s + length + extra_flank):
                used_intervals[c].append((s - extra_flank, s + length + extra_flank))
                return c, s
        return None, None

    def mutate_seq(src_seq, div):
        if div <= 0.0:
            return bytearray(src_seq)
        dst_seq = bytearray(src_seq)
        n_muts = np.random.binomial(len(src_seq), div)
        if n_muts > 0:
            offsets = np.random.choice(len(src_seq), size=n_muts, replace=False)
            for off in offsets:
                orig = dst_seq[off]
                if orig == 65: dst_seq[off] = random.choice(b'CGT')
                elif orig == 67: dst_seq[off] = random.choice(b'AGT')
                elif orig == 71: dst_seq[off] = random.choice(b'ACT')
                else: dst_seq[off] = random.choice(b'ACG')
        return dst_seq

    print(f"[INFO] Injecting {num_dups} SD families into {out_fasta}...")

    for _ in range(num_dups):
        dup_len = random.randint(min_dup_len, max_dup_len)
        div = divergence if divergence > 0.0 else random.uniform(0.0, 0.0)

        if ortholog_rate > 0.0 and random.random() < ortholog_rate and len(chrom_names) >= 3:
            # Shared insertion locus / ortholog across chromosomes
            c1, s1 = pick_free_interval(dup_len, extra_flank=flank_len)
            c2, s2 = pick_free_interval(dup_len, extra_flank=flank_len)
            c3, s3 = pick_free_interval(dup_len, extra_flank=0)
            if not c1 or not c2 or not c3:
                continue

            sd_seq = genomes[c1][s1:s1 + dup_len]
            genomes[c2][s2 - flank_len : s2 + dup_len + flank_len] = genomes[c1][s1 - flank_len : s1 + dup_len + flank_len]
            genomes[c3][s3:s3 + dup_len] = mutate_seq(sd_seq, div)

            loc1 = (c1, s1, s1 + dup_len)
            loc2 = (c2, s2, s2 + dup_len)
            loc3 = (c3, s3, s3 + dup_len)

            true_intervals.extend([loc1, loc2, loc3])
            true_pairs.append((loc1, loc3))
            true_pairs.append((loc2, loc3))
        else:
            # Standard distinct locus duplication
            c1, s1 = pick_free_interval(dup_len, extra_flank=0)
            c2, s2 = pick_free_interval(dup_len, extra_flank=0)
            if not c1 or not c2:
                continue

            genomes[c2][s2:s2 + dup_len] = mutate_seq(genomes[c1][s1:s1 + dup_len], div)

            loc1 = (c1, s1, s1 + dup_len)
            loc2 = (c2, s2, s2 + dup_len)
            true_intervals.extend([loc1, loc2])
            true_pairs.append((loc1, loc2))

    with open(out_fasta, "wb") as f:
        for chrom, seq in genomes.items():
            f.write(f">{chrom}\n".encode())
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + b"\n")

    for ext in [".fai", ".sdx"]:
        if os.path.exists(out_fasta + ext):
            os.remove(out_fasta + ext)

    return true_pairs, list(set(true_intervals)), out_fasta

def sim_run_segtrace(fasta_path, true_pairs, true_intervals, threads=8, kmer=15):
    """Runs Segtrace and computes BP, Fragment, and Pair metrics for both subclustered and unsubclustered states."""
    out_prefix = "sim_out"
    segtrace_bin = "./segtrace" if os.path.isfile("./segtrace") else "./segtrace/segtrace"
    
    start_time = time.perf_counter()
    subprocess.run([segtrace_bin, "-k", str(kmer), "-p", str(threads), fasta_path, "-o", out_prefix],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    exec_time = time.perf_counter() - start_time

    dup_bed = f"{out_prefix}.dup.bed"
    pred_intervals = parse_bed_intervals(dup_bed)
    pred_pairs_sub = load_segtrace_pairs(dup_bed, use_subclusters=True)
    pred_pairs_nosub = load_segtrace_pairs(dup_bed, use_subclusters=False)

    bp_m = calc_bp_metrics(pred_intervals, true_intervals)
    frag_m = eval_fragment_overlap(pred_intervals, true_intervals, fraction=0.5)

    pair_sub_m = evaluate_frag_pairs_fast(true_pairs, pred_pairs_sub, threshold=0.5)
    pair_nosub_m = evaluate_frag_pairs_fast(true_pairs, pred_pairs_nosub, threshold=0.5)

    # 1. Segtrace (with Subclustering)
    sub_res = {
        'Tool': 'Segtrace',
        'Recall_bp': bp_m['recall'],
        'Precision_bp': bp_m['precision'],
        'F1_bp': bp_m['f1'],
        'Recall_frag': frag_m['recall'],
        'Precision_frag': frag_m['precision'],
        'F1_frag': frag_m['f1'],
        'TP_frag': frag_m['tp'],
        'FP_frag': frag_m['fp'],
        'FN_frag': frag_m['fn'],
        'Recall_pair': pair_sub_m['recall'],
        'Precision_pair': pair_sub_m['precision'],
        'F1_pair': pair_sub_m['f1'],
        'Pairs_Count': len(pred_pairs_sub),
        'Time(s)': exec_time
    }

    # 2. Segtrace (No Subcluster / Unsubclustered)
    nosub_res = {
        'Tool': 'Segtrace (No Subcluster)',
        'Recall_bp': bp_m['recall'],
        'Precision_bp': bp_m['precision'],
        'F1_bp': bp_m['f1'],
        'Recall_frag': frag_m['recall'],
        'Precision_frag': frag_m['precision'],
        'F1_frag': frag_m['f1'],
        'TP_frag': frag_m['tp'],
        'FP_frag': frag_m['fp'],
        'FN_frag': frag_m['fn'],
        'Recall_pair': pair_nosub_m['recall'],
        'Precision_pair': pair_nosub_m['precision'],
        'F1_pair': pair_nosub_m['f1'],
        'Pairs_Count': len(pred_pairs_nosub),
        'Time(s)': exec_time
    }

    return sub_res, nosub_res

def sim_run_sedef(fasta_path, true_pairs, true_intervals, threads=8):
    """Runs SEDEF benchmark if sedef.sh is found."""
    sedef_sh = None
    if os.path.isfile("sedef/sedef.sh"):
        sedef_sh = os.path.abspath("sedef/sedef.sh")
    elif shutil.which("sedef.sh"):
        sedef_sh = shutil.which("sedef.sh")

    if not sedef_sh:
        return None

    try:
        env = os.environ.copy()
        sedef_dir = os.path.dirname(sedef_sh)
        env['PATH'] = f"{sedef_dir}:{env.get('PATH', '')}"
        sedef_out_dir = "sedef_out"
        os.makedirs(sedef_out_dir, exist_ok=True)

        start_time = time.perf_counter()
        subprocess.run([sedef_sh, "-o", sedef_out_dir, "-f", "-j", str(threads), fasta_path],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, check=True)
        exec_time = time.perf_counter() - start_time

        final_bed = f"{sedef_out_dir}/final.bed"
        if os.path.exists(final_bed):
            pred_intervals = parse_bed_intervals(final_bed)
            pred_pairs = load_bedpe_pairs(final_bed)
            bp_m = calc_bp_metrics(pred_intervals, true_intervals)
            frag_m = eval_fragment_overlap(pred_intervals, true_intervals, fraction=0.5)
            pair_m = evaluate_frag_pairs_fast(true_pairs, pred_pairs, threshold=0.5)

            return {
                'Tool': 'SEDEF',
                'Recall_bp': bp_m['recall'],
                'Precision_bp': bp_m['precision'],
                'F1_bp': bp_m['f1'],
                'Recall_frag': frag_m['recall'],
                'Precision_frag': frag_m['precision'],
                'F1_frag': frag_m['f1'],
                'TP_frag': frag_m['tp'],
                'FP_frag': frag_m['fp'],
                'FN_frag': frag_m['fn'],
                'Recall_pair': pair_m['recall'],
                'Precision_pair': pair_m['precision'],
                'F1_pair': pair_m['f1'],
                'Pairs_Count': len(pred_pairs),
                'Time(s)': exec_time
            }
    except Exception as e:
        print(f"[WARNING] SEDEF execution failed: {e}")
    return None

def sim_run_biser(fasta_path, true_pairs, true_intervals, threads=8):
    """Runs BISER benchmark if biser executable is present."""
    biser_bin = shutil.which("biser")
    if not biser_bin:
        local_biser = os.path.expanduser("~/.local/bin/biser")
        if os.path.isfile(local_biser):
            biser_bin = local_biser

    if not biser_bin:
        return None

    try:
        env = os.environ.copy()
        env['PATH'] = f"{os.path.dirname(biser_bin)}:{env.get('PATH', '')}"
        biser_out_file = "biser_out.bedpe"

        start_time = time.perf_counter()
        cmd = [biser_bin, "-o", biser_out_file, fasta_path]
        if threads > 1:
            cmd = [biser_bin, "-t", str(threads), "-o", biser_out_file, fasta_path]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, check=True)
        exec_time = time.perf_counter() - start_time

        if os.path.exists(biser_out_file):
            pred_intervals = parse_bed_intervals(biser_out_file)
            pred_pairs = load_bedpe_pairs(biser_out_file)
            bp_m = calc_bp_metrics(pred_intervals, true_intervals)
            frag_m = eval_fragment_overlap(pred_intervals, true_intervals, fraction=0.5)
            pair_m = evaluate_frag_pairs_fast(true_pairs, pred_pairs, threshold=0.5)

            return {
                'Tool': 'BISER',
                'Recall_bp': bp_m['recall'],
                'Precision_bp': bp_m['precision'],
                'F1_bp': bp_m['f1'],
                'Recall_frag': frag_m['recall'],
                'Precision_frag': frag_m['precision'],
                'F1_frag': frag_m['f1'],
                'TP_frag': frag_m['tp'],
                'FP_frag': frag_m['fp'],
                'FN_frag': frag_m['fn'],
                'Recall_pair': pair_m['recall'],
                'Precision_pair': pair_m['precision'],
                'F1_pair': pair_m['f1'],
                'Pairs_Count': len(pred_pairs),
                'Time(s)': exec_time
            }
    except Exception as e:
        print(f"[WARNING] BISER execution failed: {e}")
    return None

def plot_sim_results(results_df, output_png="evaluation_plots.png"):
    """Generates 3-panel visualization plots for simulation benchmark."""
    if results_df.empty:
        return
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

        # 1. Base-Pair (BP) Metrics
        df_bp = results_df[['Tool', 'Recall_bp', 'Precision_bp']].melt(
            id_vars='Tool', var_name='Metric', value_name='Score'
        )
        df_bp['Metric'] = df_bp['Metric'].map({'Recall_bp': 'Recall', 'Precision_bp': 'Precision'})
        sns.barplot(data=df_bp, x='Tool', y='Score', hue='Metric', ax=ax1, alpha=0.7, capsize=0.1)
        sns.stripplot(data=df_bp, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax1,
                      palette='dark:black', alpha=0.7, size=6, legend=False)
        ax1.set_ylim(0, 1.08)
        ax1.set_title('Base-Pair Level Footprint (Recall & Precision)', fontsize=13, fontweight='bold')
        ax1.set_ylabel('Score', fontsize=11)
        ax1.grid(axis='y', alpha=0.3)

        # 2. Fragment Level Metrics (50% Reciprocal Overlap)
        df_frag = results_df[['Tool', 'Recall_frag', 'Precision_frag']].melt(
            id_vars='Tool', var_name='Metric', value_name='Score'
        )
        df_frag['Metric'] = df_frag['Metric'].map({'Recall_frag': 'Recall', 'Precision_frag': 'Precision'})
        sns.barplot(data=df_frag, x='Tool', y='Score', hue='Metric', ax=ax2, alpha=0.7, capsize=0.1)
        sns.stripplot(data=df_frag, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax2,
                      palette='dark:black', alpha=0.7, size=6, legend=False)
        ax2.set_ylim(0, 1.08)
        ax2.set_title('Fragment Level (50% Reciprocal Overlap)', fontsize=13, fontweight='bold')
        ax2.set_ylabel('Score', fontsize=11)
        ax2.grid(axis='y', alpha=0.3)

        # 3. Execution Time vs Fragment F1-Score Scaling
        palette = {
            'Segtrace': '#1f77b4',
            'Segtrace (No Subcluster)': '#aec7e8',
            'SEDEF': '#d62728',
            'BISER': '#2ca02c'
        }
        markers = {
            'Segtrace': 'o',
            'Segtrace (No Subcluster)': 'X',
            'SEDEF': 's',
            'BISER': '^'
        }
        
        tools_in_df = results_df['Tool'].unique()
        tool_palette = {t: palette.get(t, '#7f7f7f') for t in tools_in_df}
        tool_markers = {t: markers.get(t, 'o') for t in tools_in_df}

        sns.scatterplot(
            data=results_df, x='Time(s)', y='F1_frag', hue='Tool', style='Tool',
            palette=tool_palette, markers=tool_markers, s=160, alpha=0.85, ax=ax3
        )

        for tool in tools_in_df:
            td = results_df[results_df['Tool'] == tool].sort_values(by='GenomeSize')
            mean_trend = td.groupby('GenomeSize')[['Time(s)', 'F1_frag']].mean().reset_index()
            ax3.plot(mean_trend['Time(s)'], mean_trend['F1_frag'],
                     color=tool_palette[tool], linestyle='--', alpha=0.6, label=f"{tool} trend")

        ax3.set_xlabel('Execution Time (seconds)', fontsize=11)
        ax3.set_ylabel('Fragment F1-Score', fontsize=11)
        ax3.set_title('F1-Score vs Execution Time Scaling', fontsize=13, fontweight='bold')
        ax3.set_ylim(0, 1.08)
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_png, dpi=300)
        plt.close()
        print(f"[INFO] Visualization plot saved to {output_png}")
    except Exception as e:
        print(f"[WARNING] Could not generate plot: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive Simulation Benchmark for Segtrace, SEDEF, and BISER with subclustering."
    )
    parser.add_argument("--genome-sizes", type=int, nargs='+', default=None,
                        help="List of genome sizes to test (e.g. --genome-sizes 10000000 20000000 50000000)")
    parser.add_argument("--genome-size", type=int, default=None,
                        help="Single genome size to test (default fallback: 10_000_000, 20_000_000, 50_000_000)")
    parser.add_argument("--reps", "-N", type=int, default=3,
                        help="Number of replicates per genome size (default: 3)")
    parser.add_argument("--num-dups", type=int, default=None,
                        help="Number of synthetic SDs to inject per replicate (default: auto based on size)")
    parser.add_argument("--divergence", type=float, default=0.0,
                        help="Sequence divergence rate for duplications (default: 0.0)")
    parser.add_argument("--ortholog-rate", type=float, default=0.2,
                        help="Fraction of SDs simulating shared insertion locus / orthologs to test subclustering (default: 0.2)")
    parser.add_argument("--threads", "-p", type=int, default=8,
                        help="Number of CPU threads to use (default: 8)")
    parser.add_argument("--kmer", "-k", type=int, default=15,
                        help="K-mer size for Segtrace (default: 15)")
    parser.add_argument("--no-sedef", action='store_true',
                        help="Skip SEDEF benchmark")
    parser.add_argument("--no-biser", action='store_true',
                        help="Skip BISER benchmark")
    parser.add_argument("--out-csv", default="evaluation_results.csv",
                        help="Output path for results CSV (default: evaluation_results.csv)")
    parser.add_argument("--out-plot", default="evaluation_plots.png",
                        help="Output path for visualization plot (default: evaluation_plots.png)")
    parser.add_argument("--no-plot", action='store_true',
                        help="Disable plot generation")

    args = parser.parse_args()

    # Determine genome sizes to evaluate
    if args.genome_sizes:
        genome_sizes = args.genome_sizes
    elif args.genome_size:
        genome_sizes = [args.genome_size]
    else:
        genome_sizes = [10_000_000, 20_000_000, 50_000_000]

    reps = max(1, args.reps)

    # Check tool availability
    biser_avail = shutil.which("biser") is not None or os.path.isfile(os.path.expanduser("~/.local/bin/biser"))
    sedef_avail = os.path.isfile("sedef/sedef.sh") or shutil.which("sedef.sh") is not None

    if not biser_avail and not args.no_biser:
        print("[INFO] BISER binary not detected in PATH. BISER benchmark will be skipped.")

    if not sedef_avail and not args.no_sedef:
        print("[INFO] SEDEF binary not detected in 'sedef/sedef.sh'. SEDEF benchmark will be skipped.")

    all_results = []

    print("\n" + "=" * 85)
    print("        SEGMENTAL DUPLICATION CALLER SIMULATION BENCHMARK")
    print("=" * 85)
    print(f" Genome Sizes:   {', '.join(f'{g:,} bp' for g in genome_sizes)}")
    print(f" Replicates:     {reps} per size")
    print(f" Ortholog Rate:  {args.ortholog_rate * 100:.1f}% (testing subclustering of shared insertion sites)")
    print(f" CPU Threads:    {args.threads}")
    print("=" * 85 + "\n")

    for g_size in genome_sizes:
        chrom_sizes = {f'chr{i}': g_size // 5 for i in range(1, 6)}
        n_dups = args.num_dups if args.num_dups else max(10, g_size // 1_000_000)

        for rep in range(1, reps + 1):
            print(f">>> [Genome Size: {g_size:,} bp | Rep {rep}/{reps} | SDs: {n_dups}] <<<")
            
            true_pairs, true_intervals, fasta_path = sim_generate_genome(
                chrom_sizes, num_dups=n_dups, divergence=args.divergence,
                ortholog_rate=args.ortholog_rate, out_fasta="sim.fa"
            )

            # 1. Segtrace (both Subclustered and No-Subcluster)
            st_sub_res, st_nosub_res = sim_run_segtrace(
                fasta_path, true_pairs, true_intervals, threads=args.threads, kmer=args.kmer
            )
            st_sub_res['GenomeSize'] = g_size
            st_sub_res['Rep'] = rep
            all_results.append(st_sub_res)

            st_nosub_res['GenomeSize'] = g_size
            st_nosub_res['Rep'] = rep
            all_results.append(st_nosub_res)

            print(f"  [Segtrace (Subclustered)] Time: {st_sub_res['Time(s)']:.2f}s | "
                  f"BP F1: {st_sub_res['F1_bp']*100:.2f}% | "
                  f"Frag F1: {st_sub_res['F1_frag']*100:.2f}% | "
                  f"Pair F1: {st_sub_res['F1_pair']*100:.2f}% (Pairs: {st_sub_res['Pairs_Count']})")
            print(f"  [Segtrace (No Subcluster)]Time: {st_nosub_res['Time(s)']:.2f}s | "
                  f"BP F1: {st_nosub_res['F1_bp']*100:.2f}% | "
                  f"Frag F1: {st_nosub_res['F1_frag']*100:.2f}% | "
                  f"Pair F1: {st_nosub_res['F1_pair']*100:.2f}% (Pairs: {st_nosub_res['Pairs_Count']})")

            # 2. SEDEF
            if not args.no_sedef and sedef_avail:
                sd_res = sim_run_sedef(fasta_path, true_pairs, true_intervals, threads=args.threads)
                if sd_res:
                    sd_res['GenomeSize'] = g_size
                    sd_res['Rep'] = rep
                    all_results.append(sd_res)
                    print(f"  [SEDEF]                   Time: {sd_res['Time(s)']:.2f}s | "
                          f"BP F1: {sd_res['F1_bp']*100:.2f}% | "
                          f"Frag F1: {sd_res['F1_frag']*100:.2f}% | "
                          f"Pair F1: {sd_res['F1_pair']*100:.2f}% (Pairs: {sd_res['Pairs_Count']})")

            # 3. BISER
            if not args.no_biser and biser_avail:
                bi_res = sim_run_biser(fasta_path, true_pairs, true_intervals, threads=args.threads)
                if bi_res:
                    bi_res['GenomeSize'] = g_size
                    bi_res['Rep'] = rep
                    all_results.append(bi_res)
                    print(f"  [BISER]                   Time: {bi_res['Time(s)']:.2f}s | "
                          f"BP F1: {bi_res['F1_bp']*100:.2f}% | "
                          f"Frag F1: {bi_res['F1_frag']*100:.2f}% | "
                          f"Pair F1: {bi_res['F1_pair']*100:.2f}% (Pairs: {bi_res['Pairs_Count']})")

            print("-" * 85)

    df_results = pd.DataFrame(all_results)
    if not df_results.empty:
        df_results.to_csv(args.out_csv, index=False)
        print(f"\n[INFO] Complete evaluation results saved to {args.out_csv}")

        # Summary Table
        summary_cols = ['Tool', 'GenomeSize', 'Recall_bp', 'Precision_bp', 'F1_bp',
                        'Recall_frag', 'Precision_frag', 'F1_frag',
                        'Recall_pair', 'Precision_pair', 'F1_pair', 'Time(s)']
        df_summary = df_results[summary_cols].groupby(['Tool', 'GenomeSize']).mean().reset_index()

        print("\n" + "=" * 115)
        print("                                     AVERAGE BENCHMARK SUMMARY")
        print("=" * 115)
        print(f"{'Tool':<26} {'Size':<12} {'BP F1':<9} {'Frag Rec':<10} {'Frag Prec':<10} {'Frag F1':<10} {'Pair Rec':<10} {'Pair Prec':<10} {'Pair F1':<10} {'Time(s)':<8}")
        print("-" * 115)
        for _, row in df_summary.iterrows():
            print(f"{row['Tool']:<26} {int(row['GenomeSize']):<12,} "
                  f"{row['F1_bp']*100:>6.2f}%   "
                  f"{row['Recall_frag']*100:>6.2f}%   {row['Precision_frag']*100:>6.2f}%   {row['F1_frag']*100:>6.2f}%   "
                  f"{row['Recall_pair']*100:>6.2f}%   {row['Precision_pair']*100:>6.2f}%   {row['F1_pair']*100:>6.2f}%   "
                  f"{row['Time(s)']:>6.2f}s")
        print("=" * 115 + "\n")

        if not args.no_plot:
            plot_sim_results(df_results, args.out_plot)

if __name__ == "__main__":
    main()
