import random
import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

def generate_simulated_genome(genome_length=50_000_000, num_dups=1000, dup_len=3000):
    bases = ['A', 'C', 'G', 'T']
    genome = [random.choice(bases) for _ in range(genome_length)]
    
    true_pairs = []
    used_intervals = []
    
    def is_overlap(s, e):
        for us, ue in used_intervals:
            if not (e <= us or s >= ue):
                return True
        return False

    for _ in range(num_dups):
        while True:
            s1 = random.randint(0, genome_length - dup_len)
            if not is_overlap(s1, s1 + dup_len):
                used_intervals.append((s1, s1 + dup_len))
                break
        
        while True:
            s2 = random.randint(0, genome_length - dup_len)
            if not is_overlap(s2, s2 + dup_len):
                used_intervals.append((s2, s2 + dup_len))
                break
        
        div = random.uniform(0.0, 0.1)
        for i in range(dup_len):
            if random.random() < div:
                orig = genome[s1 + i]
                mut = random.choice([b for b in bases if b != orig])
                genome[s2 + i] = mut
            else:
                genome[s2 + i] = genome[s1 + i]
                
        true_pairs.append(((s1, s1 + dup_len), (s2, s2 + dup_len)))
        
    with open("sim.fa", "w") as f:
        f.write(">sim\n")
        seq = "".join(genome)
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
            
    return true_pairs

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

    # TP for True Pairs (Recall denominator: total_true_bp)
    tp_true_bp = 0
    for (s1, e1), (s2, e2) in norm_true:
        L_t = max(0, e1 - s1)
        if L_t == 0: continue
        intervals = []
        for (x1, y1), (x2, y2) in norm_pred:
            # Direct orientation
            ks = max(0, x1 - s1, x2 - s2)
            ke = min(L_t, y1 - s1, y2 - s2)
            if ks < ke:
                intervals.append((ks, ke))
            # Cross orientation
            ks2 = max(0, x2 - s1, x1 - s2)
            ke2 = min(L_t, y2 - s1, y1 - s2)
            if ks2 < ke2:
                intervals.append((ks2, ke2))
        tp_true_bp += merge_intervals(intervals)

    # TP for Predicted Pairs (Precision denominator: total_pred_bp)
    tp_pred_bp = 0
    if total_pred_bp > 0:
        for (x1, y1), (x2, y2) in norm_pred:
            L_p = max(0, y1 - x1)
            if L_p == 0: continue
            intervals = []
            for (s1, e1), (s2, e2) in norm_true:
                # Direct orientation
                ks = max(0, s1 - x1, s2 - x2)
                ke = min(L_p, e1 - x1, e2 - x2)
                if ks < ke:
                    intervals.append((ks, ke))
                # Cross orientation
                ks2 = max(0, s2 - x1, s1 - x2)
                ke2 = min(L_p, e2 - x1, e1 - x2)
                if ks2 < ke2:
                    intervals.append((ks2, ke2))
            tp_pred_bp += merge_intervals(intervals)

    Sn = tp_true_bp / total_true_bp if total_true_bp > 0 else 0.0
    Pr = tp_pred_bp / total_pred_bp if total_pred_bp > 0 else 0.0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0.0
    return Sn, Pr, f1

def evaluate(true_pairs, max_dist, sub_dist, flank_ratio):
    start_time = time.time()
    subprocess.run(["./plantsds", "-d", str(max_dist), "-D", str(sub_dist), "-f", str(flank_ratio), "-p", "8", "-w", "1000", "sim.fa", "-o", "sim_out"], 
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    exec_time = time.time() - start_time
    
    clusters = {}
    if os.path.exists("sim_out.dup.bed"):
        with open("sim_out.dup.bed") as f:
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

    Sn, Pr, f1 = evaluate_bp(true_pairs, predicted_pairs)
    return Sn, Pr, f1, exec_time, predicted_pairs

def evaluate_bedpe(true_pairs, filepath):
    if not os.path.exists(filepath): return 0, 0, 0
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
    Sn, Pr, f1 = evaluate_bp(true_pairs, predicted_pairs)
    return Sn, Pr, f1


if __name__ == "__main__":
    true_pairs = generate_simulated_genome()
    print("Genome generated: sim.fa\n")

    save_bedpe(true_pairs, "true.bedpe")
    print("Saved true.bedpe\n")

    d_values = [i / 1000.0 + 0.01 for i in range(201)]
    D_values = [0.1]
    f_values = [0.1]
    
    print(f"{'sub_dist(-D)':>12} | {'max_dist(-d)':>12} | {'flank_ratio(-f)':>15} | {'Recall':>12} | {'Precision':>12} | {'F1-Score':>12} | {'Time(s)':>10}")
    print("-" * 98)
    
    results = []
    best_f1 = -1.0
    best_pred_pairs = []
    
    for sd in D_values:
        for md in d_values:
            for f_val in f_values:
                sn, pr, f1, exec_time, pred_pairs = evaluate(true_pairs, md, sd, f_val)
                if f1 > best_f1:
                    best_f1 = f1
                    best_pred_pairs = pred_pairs
                print(f"{sd:12.2f} | {md:12.2f} | {f_val:15.2f} | {sn:12.4f} | {pr:12.4f} | {f1:12.4f} | {exec_time:10.4f}")
                results.append({
                    'sub_dist': sd,
                    'max_dist': md,
                    'flank_ratio': f_val,
                    'Recall': sn,
                    'Precision': pr,
                    'F1-Score': f1,
                    'Time(s)': exec_time
                })
            print("-" * 98)

    save_bedpe(best_pred_pairs, "predict.bedpe")
    print("Saved predict.bedpe\n")

    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv("evaluation_results.csv", index=False)
    print("Saved evaluation_results.csv")

    if not df.empty:
        best_row = df.loc[df['F1-Score'].idxmax()]
        print(f"\n[Best Combination (Max F1-Score)]")
        print(f"sub_dist(-D)   : {best_row['sub_dist']:.2f}")
        print(f"max_dist(-d)   : {best_row['max_dist']:.2f}")
        print(f"flank_ratio(-f): {best_row['flank_ratio']:.2f}")
        print(f"Recall         : {best_row['Recall']:.4f}")
        print(f"Precision      : {best_row['Precision']:.4f}")
        print(f"F1-Score       : {best_row['F1-Score']:.4f}")
        print(f"Time(s)        : {best_row['Time(s)']:.4f}")

    print("\nRunning SEDEF benchmark...")
    sedef_sn, sedef_pr, sedef_f1 = 0, 0, 0
    try:
        env = os.environ.copy()
        env['PATH'] = "/opt/homebrew/bin:" + os.path.abspath("sedef") + ":" + env.get('PATH', '')
        t0 = time.time()
        subprocess.run(["./sedef/sedef.sh", "-o", "sedef_out", "-f", "-j", "8", "sim.fa"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        sedef_exec_time = time.time() - t0
        if os.path.exists("sedef_out/final.bed"):
            sedef_sn, sedef_pr, sedef_f1 = evaluate_bedpe(true_pairs, "sedef_out/final.bed")
        print(f"SEDEF -> Recall: {sedef_sn:.4f}, Pr: {sedef_pr:.4f}, F1: {sedef_f1:.4f}, Time: {sedef_exec_time:.4f}s")
    except Exception as e:
        print("SEDEF execution failed:", e)

    # Plotting Sn-Pr Curve
    plt.figure(figsize=(8, 6))
    
    df_sorted = df.sort_values('max_dist')
    
    sns.lineplot(
        data=df_sorted,
        x='Recall',
        y='Precision',
        hue='sub_dist',
        palette='tab10',
        legend=False,
        alpha=0.5,
        sort=False
    )
    
    sns.scatterplot(
        data=df, 
        x='Recall', 
        y='Precision', 
        hue='max_dist', 
        style='sub_dist', 
        palette='viridis', 
        s=100
    )
    if sedef_f1 > 0:
        plt.scatter([sedef_sn], [sedef_pr], color='red', marker='*', s=300, label='SEDEF', zorder=5)

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim(0, 1.1)
    plt.title('Recall vs Precision Curve')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('evaluation_sn_pr_curve.png', dpi=300)
    plt.close()
    print("Saved evaluation_sn_pr_curve.png")
