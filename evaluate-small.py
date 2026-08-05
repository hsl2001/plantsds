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

chrom_sizes = {
    "chr1": 30_000_000,
    "chr2": 30_000_000,
}

def generate_simulated_genome(num_dups=1000, min_dup_len=1000, max_dup_len=10_000):
    
    bases_bytes = np.frombuffer(b'ACGT', dtype=np.uint8)
    print("Generating 5 chromosomes in memory ...")
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
        
    fasta_path = f"sim_{os.getpid()}.fa"
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

def evaluate(true_pairs, max_dist, sub_dist, flank_ratio, k_val, s_val, fasta_path, chrom_offsets):
    start_time = time.time()
    segtrace_out = f"sim_out_{os.getpid()}"
    print(f"Running Segtrace with -d {max_dist} -D {sub_dist} -f {flank_ratio} -k {k_val} -s {s_val}...")
    subprocess.run(["./segtrace", "-d", str(max_dist), "-D", str(sub_dist), "-f", str(flank_ratio), "-k", str(k_val), "-s", str(s_val), "-p", "8", "-w", "1000", fasta_path, "-o", segtrace_out], 
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

    results = []
    
    print(f"Starting Benchmark: 100 Iterations with Default Parameters")
    print(f"Segtrace Defaults: -d 0.15, -D 0.20, -f 0.30, -k 21, -s 16")
    print("-" * 100)
    
    for i in range(1, 101):
        print(f"\n[Iteration {i}/100]")
        true_pairs, fasta_path, chrom_offsets = generate_simulated_genome()
        
        # Segtrace parameters
        md, sd, f_val, k_val, s_val = 0.15, 0.2, 0.3, 21, 16
        
        Sn_bp, Pr_bp, f1_bp, Sn_frag, Pr_frag, f1_frag, exec_time, pred_pairs = evaluate(
            true_pairs, md, sd, f_val, k_val, s_val, fasta_path, chrom_offsets
        )
        print(f"Segtrace -> BP F1: {f1_bp:.4f} | Frag F1: {f1_frag:.4f} | Time: {exec_time:.4f}s")
        
        results.append({
            'Iteration': i,
            'Tool': 'Segtrace',
            'Recall_bp': Sn_bp, 'Precision_bp': Pr_bp, 'F1-Score_bp': f1_bp,
            'Recall_frag': Sn_frag, 'Precision_frag': Pr_frag, 'F1-Score_frag': f1_frag,
            'Time(s)': exec_time
        })
        
        if not args.no_sedef:
            sedef_sn_bp, sedef_pr_bp, sedef_f1_bp = 0, 0, 0
            sedef_sn_frag, sedef_pr_frag, sedef_f1_frag = 0, 0, 0
            try:
                env = os.environ.copy()
                sedef_dir = os.path.abspath("sedef")
                env['PATH'] = f"{sedef_dir}:{env.get('PATH', '')}"
                sedef_out_dir = f"sedef_out_{os.getpid()}_{i}"
                t0 = time.time()
                subprocess.run([os.path.join(sedef_dir, "sedef.sh"), "-o", sedef_out_dir, "-f", "-j", "8", fasta_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
                sedef_exec_time = time.time() - t0
                if os.path.exists(f"{sedef_out_dir}/final.bed"):
                    sedef_sn_bp, sedef_pr_bp, sedef_f1_bp, sedef_sn_frag, sedef_pr_frag, sedef_f1_frag = evaluate_bedpe(true_pairs, f"{sedef_out_dir}/final.bed", chrom_offsets)
                print(f"SEDEF    -> BP F1: {sedef_f1_bp:.4f} | Frag F1: {sedef_f1_frag:.4f} | Time: {sedef_exec_time:.4f}s")
                results.append({
                    'Iteration': i,
                    'Tool': 'SEDEF',
                    'Recall_bp': sedef_sn_bp, 'Precision_bp': sedef_pr_bp, 'F1-Score_bp': sedef_f1_bp,
                    'Recall_frag': sedef_sn_frag, 'Precision_frag': sedef_pr_frag, 'F1-Score_frag': sedef_f1_frag,
                    'Time(s)': sedef_exec_time
                })
                subprocess.run(["rm", "-rf", sedef_out_dir])
            except Exception as e:
                print(f"SEDEF -> Failed to run: {e}")
                
        # Clean up fasta to save space
        for ext in ["", ".fai", ".sdx"]:
            if os.path.exists(fasta_path + ext):
                os.remove(fasta_path + ext)
        
        # Save intermediate results in case it crashes
        df = pd.DataFrame(results)
        df.to_csv("benchmark_100_results.csv", index=False)

    print("\n" + "=" * 50)
    print("Benchmark Complete!")
    print("Results saved to benchmark_100_results.csv")
    
    try:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        sns.boxplot(data=df, x='Tool', y='F1-Score_bp', ax=axes[0])
        axes[0].set_title('BP F1-Score')
        sns.boxplot(data=df, x='Tool', y='F1-Score_frag', ax=axes[1])
        axes[1].set_title('Fragment F1-Score')
        sns.boxplot(data=df, x='Tool', y='Time(s)', ax=axes[2])
        axes[2].set_title('Execution Time (s)')
        plt.tight_layout()
        plt.savefig('benchmark_100_summary.png', dpi=300)
        plt.close()
        print("Saved summary plot to benchmark_100_summary.png")
    except Exception as e:
        print(f"Failed to generate plot: {e}")
