import random
import subprocess
import os
import bisect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
import sys
import argparse
import glob

def generate_simulated_genome(chrom_sizes, num_dups=100, min_dup_len=1000, max_dup_len=10_000):
    bases_bytes = np.frombuffer(b'ACGT', dtype=np.uint8)
    genomes = {}
    for chrom, size in chrom_sizes.items():
        genomes[chrom] = bytearray(np.random.choice(bases_bytes, size=size).tobytes())
    
    true_pairs = []
    used_intervals = {chrom: [] for chrom in chrom_sizes}
    chrom_names = list(chrom_sizes.keys())
    
    def is_overlap(chrom, s, e):
        intervals = used_intervals[chrom]
        idx = bisect.bisect_left(intervals, (s, e))
        if idx > 0 and intervals[idx - 1][1] > s:
            return True
        if idx < len(intervals) and intervals[idx][0] < e:
            return True
        return False

    def add_interval(chrom, s, e):
        bisect.insort_left(used_intervals[chrom], (s, e))

    print(f"Injecting {num_dups} duplications (inter and intra-chromosomal)...")
    for _ in range(num_dups):
        dup_len = random.randint(min_dup_len, max_dup_len)
        
        while True:
            c1 = random.choice(chrom_names)
            s1 = random.randint(0, chrom_sizes[c1] - dup_len)
            if not is_overlap(c1, s1, s1 + dup_len):
                add_interval(c1, s1, s1 + dup_len)
                break
        
        while True:
            c2 = random.choice(chrom_names)
            s2 = random.randint(0, chrom_sizes[c2] - dup_len)
            if not is_overlap(c2, s2, s2 + dup_len):
                add_interval(c2, s2, s2 + dup_len)
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
        
    fasta_path = "sim.fa"
    with open(fasta_path, "wb") as f:
        for chrom, seq in genomes.items():
            f.write(f">{chrom}\n".encode())
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + b"\n")
            
    for ext in [".fai", ".sdx"]:
        if os.path.exists(fasta_path + ext):
            os.remove(fasta_path + ext)
            
    global_offset = 0
    chrom_offsets = {}
    for chrom, size in chrom_sizes.items():
        chrom_offsets[chrom] = global_offset
        global_offset += size + 100_000_000
            
    return true_pairs, fasta_path, chrom_offsets

def to_global_pairs(pairs, chrom_offsets):
    global_pairs = []
    for (c1, s1, e1), (c2, s2, e2) in pairs:
        g1_s = chrom_offsets[c1] + s1
        g1_e = chrom_offsets[c1] + e1
        g2_s = chrom_offsets[c2] + s2
        g2_e = chrom_offsets[c2] + e2
        global_pairs.append(((g1_s, g1_e), (g2_s, g2_e)))
    return global_pairs

def save_bedpe(pairs, filepath):
    norm_pairs = []
    for (c1, s1, e1), (c2, s2, e2) in pairs:
        if c1 > c2 or (c1 == c2 and (s1 > s2 or (s1 == s2 and e1 > e2))):
            norm_pairs.append(((c2, s2, e2), (c1, s1, e1)))
        else:
            norm_pairs.append(((c1, s1, e1), (c2, s2, e2)))
    norm_pairs.sort(key=lambda p: (p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]))
    with open(filepath, "w") as f:
        for (c1, s1, e1), (c2, s2, e2) in norm_pairs:
            f.write(f"{c1}\t{s1}\t{e1}\t{c2}\t{s2}\t{e2}\n")
    return norm_pairs

def merge_intervals(intervals):
    if not intervals:
        return 0
    intervals.sort(key=lambda x: x[0])
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return sum(e - s for s, e in merged)

def evaluate_bp(true_pairs, predicted_pairs):
    if not true_pairs:
        return 0.0, 0.0, 0.0

    norm_true = []
    for (s1, e1), (s2, e2) in true_pairs:
        if s1 > s2 or (s1 == s2 and e1 > e2):
            norm_true.append(((s2, e2), (s1, e1)))
        else:
            norm_true.append(((s1, e1), (s2, e2)))
            
    norm_pred = []
    for (x1, y1), (x2, y2) in predicted_pairs:
        if x1 > x2 or (x1 == x2 and y1 > y2):
            norm_pred.append(((x2, y2), (x1, y1)))
        else:
            norm_pred.append(((x1, y1), (x2, y2)))

    total_true_bp = sum(max(0, e1 - s1) for (s1, e1), (s2, e2) in norm_true)
    total_pred_bp = sum(max(0, y1 - x1) for (x1, y1), (x2, y2) in norm_pred)

    if total_true_bp == 0:
        return 0.0, 0.0, 0.0

    if not norm_pred:
        return 0.0, 0.0, 0.0

    PX1 = np.array([p[0][0] for p in norm_pred], dtype=np.int64)
    PY1 = np.array([p[0][1] for p in norm_pred], dtype=np.int64)
    PX2 = np.array([p[1][0] for p in norm_pred], dtype=np.int64)
    PY2 = np.array([p[1][1] for p in norm_pred], dtype=np.int64)

    # TP for True Pairs (Recall denominator: total_true_bp)
    tp_true_bp = 0
    for (s1, e1), (s2, e2) in norm_true:
        L_t = max(0, e1 - s1)
        if L_t == 0: continue

        cand_mask = (
            ((PX1 < e1) & (PY1 > s1) & (PX2 < s2 + L_t) & (PY2 > s2)) |
            ((PX2 < e1) & (PY2 > s1) & (PX1 < s2 + L_t) & (PY1 > s2))
        )
        if not np.any(cand_mask):
            continue

        px1_c, py1_c, px2_c, py2_c = PX1[cand_mask], PY1[cand_mask], PX2[cand_mask], PY2[cand_mask]

        # Direct orientation
        ks = np.maximum(0, np.maximum(px1_c - s1, px2_c - s2))
        ke = np.minimum(L_t, np.minimum(py1_c - s1, py2_c - s2))
        valid_dir = ks < ke

        # Cross orientation
        ks2 = np.maximum(0, np.maximum(px2_c - s1, px1_c - s2))
        ke2 = np.minimum(L_t, np.minimum(py2_c - s1, py1_c - s2))
        valid_cross = ks2 < ke2

        intervals = []
        if np.any(valid_dir):
            intervals.extend(zip(ks[valid_dir].tolist(), ke[valid_dir].tolist()))
        if np.any(valid_cross):
            intervals.extend(zip(ks2[valid_cross].tolist(), ke2[valid_cross].tolist()))

        tp_true_bp += merge_intervals(intervals)

    # TP for Predicted Pairs (Precision denominator: total_pred_bp)
    tp_pred_bp = 0
    if total_pred_bp > 0:
        TX1 = np.array([t[0][0] for t in norm_true], dtype=np.int64)
        TY1 = np.array([t[0][1] for t in norm_true], dtype=np.int64)
        TX2 = np.array([t[1][0] for t in norm_true], dtype=np.int64)
        TY2 = np.array([t[1][1] for t in norm_true], dtype=np.int64)

        for (x1, y1), (x2, y2) in norm_pred:
            L_p = max(0, y1 - x1)
            if L_p == 0: continue

            cand_mask = (
                ((TX1 < y1) & (TY1 > x1) & (TX2 < x2 + L_p) & (TY2 > x2)) |
                ((TX2 < y1) & (TY2 > x1) & (TX1 < x2 + L_p) & (TY1 > x2))
            )
            if not np.any(cand_mask):
                continue

            tx1_c, ty1_c, tx2_c, ty2_c = TX1[cand_mask], TY1[cand_mask], TX2[cand_mask], TY2[cand_mask]

            # Direct orientation
            ks = np.maximum(0, np.maximum(tx1_c - x1, tx2_c - x2))
            ke = np.minimum(L_p, np.minimum(ty1_c - x1, ty2_c - x2))
            valid_dir = ks < ke

            # Cross orientation
            ks2 = np.maximum(0, np.maximum(tx2_c - x1, tx1_c - x2))
            ke2 = np.minimum(L_p, np.minimum(ty2_c - x1, ty1_c - x2))
            valid_cross = ks2 < ke2

            intervals = []
            if np.any(valid_dir):
                intervals.extend(zip(ks[valid_dir].tolist(), ke[valid_dir].tolist()))
            if np.any(valid_cross):
                intervals.extend(zip(ks2[valid_cross].tolist(), ke2[valid_cross].tolist()))

            tp_pred_bp += merge_intervals(intervals)

    Sn = tp_true_bp / total_true_bp if total_true_bp > 0 else 0.0
    Pr = tp_pred_bp / total_pred_bp if total_pred_bp > 0 else 0.0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0.0
    return Sn, Pr, f1

def evaluate_frag(true_pairs, predicted_pairs, threshold=0.5):
    if not true_pairs or not predicted_pairs:
        return 0.0, 0.0, 0.0

    norm_true = []
    for (s1, e1), (s2, e2) in true_pairs:
        if s1 > s2 or (s1 == s2 and e1 > e2):
            norm_true.append(((s2, e2), (s1, e1)))
        else:
            norm_true.append(((s1, e1), (s2, e2)))
            
    norm_pred = []
    for (x1, y1), (x2, y2) in predicted_pairs:
        if x1 > x2 or (x1 == x2 and y1 > y2):
            norm_pred.append(((x2, y2), (x1, y1)))
        else:
            norm_pred.append(((x1, y1), (x2, y2)))

    N = len(norm_true)
    M = len(norm_pred)
    if N == 0 or M == 0:
        return 0.0, 0.0, 0.0

    PX1 = np.array([p[0][0] for p in norm_pred], dtype=np.float64)
    PY1 = np.array([p[0][1] for p in norm_pred], dtype=np.float64)
    PX2 = np.array([p[1][0] for p in norm_pred], dtype=np.float64)
    PY2 = np.array([p[1][1] for p in norm_pred], dtype=np.float64)
    LP1 = PY1 - PX1
    LP2 = PY2 - PX2

    matched_true = set()
    matched_pred = set()

    for i, ((s1, e1), (s2, e2)) in enumerate(norm_true):
        Lt1 = float(e1 - s1)
        Lt2 = float(e2 - s2)
        if Lt1 == 0 or Lt2 == 0: continue

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

        matched_p_indices = np.where(dir_match | cross_match)[0]
        if len(matched_p_indices) > 0:
            matched_true.add(i)
            matched_pred.update(matched_p_indices.tolist())

    Sn = len(matched_true) / N
    Pr = len(matched_pred) / M
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0.0
    return Sn, Pr, f1

def evaluate(true_pairs, fasta_path, chrom_offsets):
    start_time = time.time()
    segtrace_out = "sim_out"
    subprocess.run(["./segtrace", "-p", "8", fasta_path, "-o", segtrace_out], 
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    exec_time = time.time() - start_time
    
    clusters = {}
    if os.path.exists(f"{segtrace_out}.dup.bed"):
        with open(f"{segtrace_out}.dup.bed") as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split()
                if len(parts) < 5: continue
                chrom, start, end, cluster_id, subcluster_id = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
                if "-" in chrom: chrom = chrom.split("-", 1)[-1]
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append((chrom, start, end, subcluster_id))
                
    predicted_pairs = []
    for cluster_id, regions in clusters.items():
        for i in range(len(regions)):
            for j in range(i + 1, len(regions)):
                ra_c, ra_s, ra_e, ra_sub = regions[i]
                rb_c, rb_s, rb_e, rb_sub = regions[j]
                if ra_sub == rb_sub:
                    continue
                if ra_c == rb_c and max(ra_s, rb_s) < min(ra_e, rb_e):
                    continue
                predicted_pairs.append(((ra_c, ra_s, ra_e), (rb_c, rb_s, rb_e)))

    global_true = to_global_pairs(true_pairs, chrom_offsets)
    global_pred = to_global_pairs(predicted_pairs, chrom_offsets)

    Sn_bp, Pr_bp, f1_bp = evaluate_bp(global_true, global_pred)
    Sn_frag, Pr_frag, f1_frag = evaluate_frag(global_true, global_pred)
    return Sn_bp, Pr_bp, f1_bp, Sn_frag, Pr_frag, f1_frag, exec_time, predicted_pairs

def evaluate_bedpe(true_pairs, filepath, chrom_offsets):
    if not os.path.exists(filepath): return 0, 0, 0, 0, 0, 0
    predicted_pairs = []
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
                predicted_pairs.append(((c1, ra_s, ra_e), (c2, rb_s, rb_e)))
            except ValueError: continue
            
    save_bedpe(predicted_pairs, "sedef_predict.bedpe")
    
    global_true = to_global_pairs(true_pairs, chrom_offsets)
    global_pred = to_global_pairs(predicted_pairs, chrom_offsets)
    
    Sn_bp, Pr_bp, f1_bp = evaluate_bp(global_true, global_pred)
    Sn_frag, Pr_frag, f1_frag = evaluate_frag(global_true, global_pred)
    return Sn_bp, Pr_bp, f1_bp, Sn_frag, Pr_frag, f1_frag


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate segtrace performance.")
    parser.add_argument('--no-sedef', action='store_true', help="Skip SEDEF benchmark")
    args = parser.parse_args()

    all_results = []

    genome_sizes_to_test = [100_000_00, 200_000_00, 500_000_00] * 5
    # genome_sizes_to_test = [100_000_000]

    for g_size in genome_sizes_to_test:
        print(f"\n======================================")
        print(f"Testing Genome Size: {g_size:,} bp")
        print(f"======================================")
            
        chrom_sizes = {"chr1": g_size}
        num_dups = num_dups = g_size // 500_000
        
        true_pairs, fasta_path, chrom_offsets = generate_simulated_genome(chrom_sizes, num_dups=num_dups)
        
        # Segtrace
        Sn_bp, Pr_bp, f1_bp, Sn_frag, Pr_frag, f1_frag, exec_time, pred_pairs = evaluate(true_pairs, fasta_path, chrom_offsets)
        all_results.append({
            'Tool': 'Segtrace',
            'GenomeSize': g_size,
            'Recall_bp': Sn_bp,
            'Precision_bp': Pr_bp,
            'F1-Score_bp': f1_bp,
            'Recall_frag': Sn_frag,
            'Precision_frag': Pr_frag,
            'F1-Score_frag': f1_frag,
            'Time(s)': exec_time
        })
            
        if not args.no_sedef:
            # SEDEF
            try:
                env = os.environ.copy()
                sedef_dir = os.path.abspath("sedef")
                env['PATH'] = f"{sedef_dir}:{env.get('PATH', '')}"
                sedef_out_dir = "sedef_out"
                t0 = time.time()
                subprocess.run([os.path.join(sedef_dir, "sedef.sh"), "-o", sedef_out_dir, "-f", "-j", "8", fasta_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
                sedef_exec_time = time.time() - t0
                if os.path.exists(f"{sedef_out_dir}/final.bed"):
                    sedef_sn_bp, sedef_pr_bp, sedef_f1_bp, sedef_sn_frag, sedef_pr_frag, sedef_f1_frag = evaluate_bedpe(true_pairs, f"{sedef_out_dir}/final.bed", chrom_offsets)
                    all_results.append({
                        'Tool': 'SEDEF',
                        'GenomeSize': g_size,
                        'Recall_bp': sedef_sn_bp,
                        'Precision_bp': sedef_pr_bp,
                        'F1-Score_bp': sedef_f1_bp,
                        'Recall_frag': sedef_sn_frag,
                        'Precision_frag': sedef_pr_frag,
                        'F1-Score_frag': sedef_f1_frag,
                        'Time(s)': sedef_exec_time
                    })
            except Exception as e:
                print(f"SEDEF Failed: {e}")
                
            # BISER
            try:
                env = os.environ.copy()
                env['PATH'] = f"{os.path.expanduser('~/.local/bin')}:{env.get('PATH', '')}"
                biser_out_file = "biser_out.bedpe"
                t0 = time.time()
                subprocess.run(["biser", "-o", biser_out_file, fasta_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
                biser_exec_time = time.time() - t0
                if os.path.exists(biser_out_file):
                    biser_sn_bp, biser_pr_bp, biser_f1_bp, biser_sn_frag, biser_pr_frag, biser_f1_frag = evaluate_bedpe(true_pairs, biser_out_file, chrom_offsets)
                    all_results.append({
                        'Tool': 'BISER',
                        'GenomeSize': g_size,
                        'Recall_bp': biser_sn_bp,
                        'Precision_bp': biser_pr_bp,
                        'F1-Score_bp': biser_f1_bp,
                        'Recall_frag': biser_sn_frag,
                        'Precision_frag': biser_pr_frag,
                        'F1-Score_frag': biser_f1_frag,
                        'Time(s)': biser_exec_time
                    })
            except Exception as e:
                print(f"BISER Failed: {e}")
        
        # Clean up fasta for this iteration
        #for ext in ["", ".fai", ".sdx"]:
        #    if os.path.exists(fasta_path + ext):
        #        os.remove(fasta_path + ext)

    df_all = pd.DataFrame(all_results)
    if not df_all.empty:
        df_all.to_csv("evaluation_results.csv", index=False)
        print("\nSaved evaluation_results.csv")
    else:
        print("No results to save.")
        sys.exit(0)

    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

    # A1) BP Level Bar Plot with Jitter
    df_bp = df_all[['Tool', 'Recall_bp', 'Precision_bp']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
    df_bp['Metric'] = df_bp['Metric'].map({'Recall_bp': 'Recall', 'Precision_bp': 'Precision'})
    
    sns.barplot(data=df_bp, x='Tool', y='Score', hue='Metric', ax=ax1, alpha=0.6, capsize=.1)
    sns.stripplot(data=df_bp, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax1, palette='dark:black', alpha=0.7, size=5, legend=False)
    
    ax1.set_ylim(0, 1.1)
    ax1.set_title('BP-level Recall & Precision')
    ax1.grid(axis='y', alpha=0.3)

    # A2) Fragment Level Bar Plot with Jitter
    df_frag = df_all[['Tool', 'Recall_frag', 'Precision_frag']].melt(id_vars='Tool', var_name='Metric', value_name='Score')
    df_frag['Metric'] = df_frag['Metric'].map({'Recall_frag': 'Recall', 'Precision_frag': 'Precision'})
    
    sns.barplot(data=df_frag, x='Tool', y='Score', hue='Metric', ax=ax2, alpha=0.6, capsize=.1)
    sns.stripplot(data=df_frag, x='Tool', y='Score', hue='Metric', dodge=True, ax=ax2, palette='dark:black', alpha=0.7, size=5, legend=False)
    
    ax2.set_ylim(0, 1.1)
    ax2.set_title('Fragment-level Recall & Precision')
    ax2.grid(axis='y', alpha=0.3)

    # B) F1 vs Time Tradeoff (Scatter Plot)
    palette = {'Segtrace': 'blue', 'SEDEF': 'red', 'BISER': 'green'}
    markers = {'Segtrace': 'o', 'SEDEF': '*', 'BISER': 's'}
    
    sns.scatterplot(data=df_all, x='Time(s)', y='F1-Score_frag', hue='Tool', style='Tool', 
                    palette=palette, markers=markers, s=150, alpha=0.8, ax=ax3)
    
    # Draw lines connecting the points for each tool to show scaling trend
    for tool in df_all['Tool'].unique():
        tool_data = df_all[df_all['Tool'] == tool].sort_values(by='GenomeSize')
        ax3.plot(tool_data['Time(s)'], tool_data['F1-Score_frag'], color=palette[tool], alpha=0.4)

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
