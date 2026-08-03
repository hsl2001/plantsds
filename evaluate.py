import random
import subprocess
import os
import bisect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

def generate_simulated_genome(genome_length=100_000_000, num_dups=1000, min_dup_len=1000, max_dup_len=10_000):
    bases_bytes = np.frombuffer(b'ACGT', dtype=np.uint8)
    genome = bytearray(np.random.choice(bases_bytes, size=genome_length).tobytes())
    
    true_pairs = []
    used_intervals = []
    
    def is_overlap(s, e):
        idx = bisect.bisect_left(used_intervals, (s, e))
        if idx > 0 and used_intervals[idx - 1][1] > s:
            return True
        if idx < len(used_intervals) and used_intervals[idx][0] < e:
            return True
        return False

    def add_interval(s, e):
        bisect.insort_left(used_intervals, (s, e))

    for _ in range(num_dups):
        dup_len = random.randint(min_dup_len, max_dup_len)
        while True:
            s1 = random.randint(0, genome_length - dup_len)
            if not is_overlap(s1, s1 + dup_len):
                add_interval(s1, s1 + dup_len)
                break
        
        while True:
            s2 = random.randint(0, genome_length - dup_len)
            if not is_overlap(s2, s2 + dup_len):
                add_interval(s2, s2 + dup_len)
                break
        
        genome[s2:s2 + dup_len] = genome[s1:s1 + dup_len]

        div = random.uniform(0.0, 0.1)
        num_muts = np.random.binomial(dup_len, div)
        if num_muts > 0:
            mut_offsets = np.random.choice(dup_len, size=num_muts, replace=False)
            for offset in mut_offsets:
                orig = genome[s1 + offset]
                if orig == 65: mut = random.choice(b'CGT')
                elif orig == 67: mut = random.choice(b'AGT')
                elif orig == 71: mut = random.choice(b'ACT')
                else: mut = random.choice(b'ACG')
                genome[s2 + offset] = mut
                
        true_pairs.append(((s1, s1 + dup_len), (s2, s2 + dup_len)))
        
    fasta_path = f"sim_{os.getpid()}.fa"
    with open(fasta_path, "wb") as f:
        f.write(b">sim\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + b"\n")
            
    for ext in [".fai", ".sdx"]:
        if os.path.exists(fasta_path + ext):
            os.remove(fasta_path + ext)
            
    return true_pairs, fasta_path

def save_bedpe(pairs, filepath, chrom="sim"):
    norm_pairs = []
    for (s1, e1), (s2, e2) in pairs:
        if s1 > s2 or (s1 == s2 and e1 > e2):
            norm_pairs.append(((s2, e2), (s1, e1)))
        else:
            norm_pairs.append(((s1, e1), (s2, e2)))
    norm_pairs.sort(key=lambda p: (p[0][0], p[0][1], p[1][0], p[1][1]))
    with open(filepath, "w") as f:
        for (s1, e1), (s2, e2) in norm_pairs:
            f.write(f"{chrom}\t{s1}\t{e1}\t{chrom}\t{s2}\t{e2}\n")
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

def evaluate(true_pairs, max_dist, sub_dist, flank_ratio, fasta_path):
    start_time = time.time()
    plantsds_out = f"sim_out_{os.getpid()}"
    subprocess.run(["./plantsds", "-d", str(max_dist), "-D", str(sub_dist), "-f", str(flank_ratio), "-p", "8", "-w", "1000", fasta_path, "-o", plantsds_out], 
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    exec_time = time.time() - start_time
    
    clusters = {}
    if os.path.exists(f"{plantsds_out}.dup.bed"):
        with open(f"{plantsds_out}.dup.bed") as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split()
                if len(parts) < 5: continue
                chrom, start, end, cluster_id, subcluster_id = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append((start, end, subcluster_id))
                
    predicted_pairs = []
    for cluster_id, regions in clusters.items():
        for i in range(len(regions)):
            for j in range(i + 1, len(regions)):
                ra_s, ra_e, ra_sub = regions[i]
                rb_s, rb_e, rb_sub = regions[j]
                predicted_pairs.append(((ra_s, ra_e), (rb_s, rb_e)))

    Sn_bp, Pr_bp, f1_bp = evaluate_bp(true_pairs, predicted_pairs)
    Sn_frag, Pr_frag, f1_frag = evaluate_frag(true_pairs, predicted_pairs)
    return Sn_bp, Pr_bp, f1_bp, Sn_frag, Pr_frag, f1_frag, exec_time, predicted_pairs

def evaluate_bedpe(true_pairs, filepath):
    if not os.path.exists(filepath): return 0, 0, 0, 0, 0, 0
    predicted_pairs = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split()
            if len(parts) < 6: continue
            try:
                ra_s, ra_e = int(parts[1]), int(parts[2])
                rb_s, rb_e = int(parts[4]), int(parts[5])
                predicted_pairs.append(((ra_s, ra_e), (rb_s, rb_e)))
            except ValueError: continue
            
    save_bedpe(predicted_pairs, "sedef_predict.bedpe")
    Sn_bp, Pr_bp, f1_bp = evaluate_bp(true_pairs, predicted_pairs)
    Sn_frag, Pr_frag, f1_frag = evaluate_frag(true_pairs, predicted_pairs)
    return Sn_bp, Pr_bp, f1_bp, Sn_frag, Pr_frag, f1_frag


if __name__ == "__main__":
    true_pairs, fasta_path = generate_simulated_genome()
    print("Genome generated\n")

    save_bedpe(true_pairs, "true.bedpe")
    print("Saved true.bedpe\n")

    d_values = [i / 1000.0 + 0.01 for i in range(201)]
    D_values = [0.1]
    f_values = [0.1]
    
    print(f"{'sub_dist(-D)':>12} | {'max_dist(-d)':>12} | {'flank_ratio(-f)':>15} | {'BP Sn':>10} | {'BP Pr':>10} | {'BP F1':>10} | {'Frag Sn':>10} | {'Frag Pr':>10} | {'Frag F1':>10} | {'Time(s)':>8}")
    print("-" * 135)
    
    results = []
    best_f1_bp = -1.0
    best_pred_pairs = []
    
    for sd in D_values:
        for md in d_values:
            for f_val in f_values:
                Sn_bp, Pr_bp, f1_bp, Sn_frag, Pr_frag, f1_frag, exec_time, pred_pairs = evaluate(true_pairs, md, sd, f_val, fasta_path)
                if f1_bp > best_f1_bp:
                    best_f1_bp = f1_bp
                    best_pred_pairs = pred_pairs
                print(f"{sd:12.2f} | {md:12.2f} | {f_val:15.2f} | {Sn_bp:10.4f} | {Pr_bp:10.4f} | {f1_bp:10.4f} | {Sn_frag:10.4f} | {Pr_frag:10.4f} | {f1_frag:10.4f} | {exec_time:8.4f}")
                results.append({
                    'sub_dist': sd,
                    'max_dist': md,
                    'flank_ratio': f_val,
                    'Recall_bp': Sn_bp,
                    'Precision_bp': Pr_bp,
                    'F1-Score_bp': f1_bp,
                    'Recall_frag': Sn_frag,
                    'Precision_frag': Pr_frag,
                    'F1-Score_frag': f1_frag,
                    'Time(s)': exec_time
                })
            print("-" * 135)

    save_bedpe(best_pred_pairs, "predict.bedpe")
    print("Saved predict.bedpe\n")

    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv("evaluation_results.csv", index=False)
    print("Saved evaluation_results.csv")

    if not df.empty:
        best_row_bp = df.loc[df['F1-Score_bp'].idxmax()]
        print(f"\n[Best Combination (Max BP F1-Score)]")
        print(f"sub_dist(-D)   : {best_row_bp['sub_dist']:.2f}")
        print(f"max_dist(-d)   : {best_row_bp['max_dist']:.2f}")
        print(f"flank_ratio(-f): {best_row_bp['flank_ratio']:.2f}")
        print(f"BP Recall      : {best_row_bp['Recall_bp']:.4f}")
        print(f"BP Precision   : {best_row_bp['Precision_bp']:.4f}")
        print(f"BP F1-Score    : {best_row_bp['F1-Score_bp']:.4f}")
        print(f"Frag Recall    : {best_row_bp['Recall_frag']:.4f}")
        print(f"Frag Precision : {best_row_bp['Precision_frag']:.4f}")
        print(f"Frag F1-Score  : {best_row_bp['F1-Score_frag']:.4f}")
        print(f"Time(s)        : {best_row_bp['Time(s)']:.4f}")

    print("\nRunning SEDEF benchmark...")
    sedef_sn_bp, sedef_pr_bp, sedef_f1_bp = 0, 0, 0
    sedef_sn_frag, sedef_pr_frag, sedef_f1_frag = 0, 0, 0
    try:
        env = os.environ.copy()
        sedef_dir = os.path.abspath("sedef")
        env['PATH'] = f"{sedef_dir}:{env.get('PATH', '')}"
        sedef_out_dir = f"sedef_out_{os.getpid()}"
        t0 = time.time()
        subprocess.run([os.path.join(sedef_dir, "sedef.sh"), "-o", sedef_out_dir, "-f", "-j", "8", fasta_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        sedef_exec_time = time.time() - t0
        if os.path.exists(f"{sedef_out_dir}/final.bed"):
            sedef_sn_bp, sedef_pr_bp, sedef_f1_bp, sedef_sn_frag, sedef_pr_frag, sedef_f1_frag = evaluate_bedpe(true_pairs, f"{sedef_out_dir}/final.bed")
        print(f"SEDEF -> BP Recall: {sedef_sn_bp:.4f}, Pr: {sedef_pr_bp:.4f}, F1: {sedef_f1_bp:.4f} | Frag Recall: {sedef_sn_frag:.4f}, Pr: {sedef_pr_frag:.4f}, F1: {sedef_f1_frag:.4f}, Time: {sedef_exec_time:.4f}s")
    except Exception as e:
        print(f"SEDEF -> Failed to run: {e}")

    # Plotting Sn-Pr Curve (BP and Fragment side-by-side)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    df_sorted = df.sort_values('max_dist')

    # BP plot
    sns.lineplot(data=df_sorted, x='Recall_bp', y='Precision_bp', hue='sub_dist', palette='tab10', legend=False, alpha=0.5, sort=False, ax=ax1)
    sns.scatterplot(data=df, x='Recall_bp', y='Precision_bp', hue='max_dist', style='sub_dist', palette='viridis', s=100, ax=ax1)
    if sedef_f1_bp > 0:
        ax1.scatter([sedef_sn_bp], [sedef_pr_bp], color='red', marker='*', s=300, label='SEDEF', zorder=5)
    ax1.set_xlabel('Recall (BP)')
    ax1.set_ylabel('Precision (BP)')
    ax1.set_ylim(0, 1.1)
    ax1.set_title('BP-level Recall vs Precision')
    ax1.grid(True)

    # Frag plot
    sns.lineplot(data=df_sorted, x='Recall_frag', y='Precision_frag', hue='sub_dist', palette='tab10', legend=False, alpha=0.5, sort=False, ax=ax2)
    sns.scatterplot(data=df, x='Recall_frag', y='Precision_frag', hue='max_dist', style='sub_dist', palette='viridis', s=100, ax=ax2)
    if sedef_f1_frag > 0:
        ax2.scatter([sedef_sn_frag], [sedef_pr_frag], color='red', marker='*', s=300, label='SEDEF', zorder=5)
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.set_ylim(0, 1.1)
    ax2.set_title('Fragment-level Recall vs Precision')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('evaluation_sn_pr_curve.png', dpi=300)
    plt.close()
    print("Saved evaluation_sn_pr_curve.png")
