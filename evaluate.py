import random
import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

def generate_simulated_genome(genome_length=5000000, num_dups=200, dup_len=3000):
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

def is_match(ra_s, ra_e, rb_s, rb_e, ta_s, ta_e, tb_s, tb_e, tol=100):
    match1 = abs(ra_s - ta_s) <= tol and abs(ra_e - ta_e) <= tol and abs(rb_s - tb_s) <= tol and abs(rb_e - tb_e) <= tol
    match2 = abs(ra_s - tb_s) <= tol and abs(ra_e - tb_e) <= tol and abs(rb_s - ta_s) <= tol and abs(rb_e - ta_e) <= tol
    return match1 or match2

def merge_predicted_pairs(pairs, dist=100):
    normalized = []
    for p in pairs:
        r1, r2 = p
        if r1[0] > r2[0]:
            r1, r2 = r2, r1
        normalized.append([list(r1), list(r2)])
        
    changed = True
    while changed:
        changed = False
        new_pairs = []
        used = [False] * len(normalized)
        for i in range(len(normalized)):
            if used[i]: continue
            p1 = normalized[i]
            for j in range(i + 1, len(normalized)):
                if used[j]: continue
                p2 = normalized[j]
                
                def can_merge(a, b):
                    return (min(a[1], b[1]) - max(a[0], b[0])) >= -dist
                
                if can_merge(p1[0], p2[0]) and can_merge(p1[1], p2[1]):
                    p1[0][0] = min(p1[0][0], p2[0][0])
                    p1[0][1] = max(p1[0][1], p2[0][1])
                    p1[1][0] = min(p1[1][0], p2[1][0])
                    p1[1][1] = max(p1[1][1], p2[1][1])
                    used[j] = True
                    changed = True
            new_pairs.append(p1)
        normalized = new_pairs
    return [((p[0][0], p[0][1]), (p[1][0], p[1][1])) for p in normalized]


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
                # Store (start, end, subcluster_id)
                clusters[cluster_id].append((start, end, subcluster_id))
                
    predicted_pairs = []
    for cluster_id, regions in clusters.items():
        for i in range(len(regions)):
            for j in range(i + 1, len(regions)):
                ra_s, ra_e, ra_sub = regions[i]
                rb_s, rb_e, rb_sub = regions[j]
                predicted_pairs.append(((ra_s, ra_e), (rb_s, rb_e)))
                
    predicted_pairs = merge_predicted_pairs(predicted_pairs)
    
    T_count = len(true_pairs)
    P_count = len(predicted_pairs)
    
    true_found = 0
    for t_pair in true_pairs:
        (ta_s, ta_e), (tb_s, tb_e) = t_pair
        found = False
        for p_pair in predicted_pairs:
            (ra_s, ra_e), (rb_s, rb_e) = p_pair
            if is_match(ra_s, ra_e, rb_s, rb_e, ta_s, ta_e, tb_s, tb_e, 1000):
                found = True
                break
        if found:
            true_found += 1
            
    pred_correct = 0
    for p_pair in predicted_pairs:
        (ra_s, ra_e), (rb_s, rb_e) = p_pair
        correct = False
        for t_pair in true_pairs:
            (ta_s, ta_e), (tb_s, tb_e) = t_pair
            if is_match(ra_s, ra_e, rb_s, rb_e, ta_s, ta_e, tb_s, tb_e, 1000):
                correct = True
                break
        if correct:
            pred_correct += 1
                    
    Sn = true_found / T_count if T_count > 0 else 0
    Pr = pred_correct / P_count if P_count > 0 else 0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0
    return Sn, Pr, f1, exec_time

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
            
    predicted_pairs = merge_predicted_pairs(predicted_pairs)
            
    T_count = len(true_pairs)
    P_count = len(predicted_pairs)
    
    true_found = 0
    for t_pair in true_pairs:
        (ta_s, ta_e), (tb_s, tb_e) = t_pair
        found = False
        for p_pair in predicted_pairs:
            (ra_s, ra_e), (rb_s, rb_e) = p_pair
            if is_match(ra_s, ra_e, rb_s, rb_e, ta_s, ta_e, tb_s, tb_e, 1000):
                found = True
                break
        if found:
            true_found += 1
            
    pred_correct = 0
    for p_pair in predicted_pairs:
        (ra_s, ra_e), (rb_s, rb_e) = p_pair
        correct = False
        for t_pair in true_pairs:
            (ta_s, ta_e), (tb_s, tb_e) = t_pair
            if is_match(ra_s, ra_e, rb_s, rb_e, ta_s, ta_e, tb_s, tb_e, 1000):
                correct = True
                break
        if correct:
            pred_correct += 1
            
    Sn = true_found / T_count if T_count > 0 else 0
    Pr = pred_correct / P_count if P_count > 0 else 0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0
    return Sn, Pr, f1


if __name__ == "__main__":
    print("Generating simulated genome (5Mb) with 200 duplicates...")
    true_pairs = generate_simulated_genome()
    print("Genome generated: sim.fa\n")

    d_values = [i / 1000.0 + 0.01 for i in range(251)]
    D_values = [0.1]
    f_values = [0.1]
    
    print(f"{'sub_dist(-D)':>12} | {'max_dist(-d)':>12} | {'flank_ratio(-f)':>15} | {'Recall':>12} | {'Precision':>12} | {'F1-Score':>12} | {'Time(s)':>10}")
    print("-" * 98)
    
    results = []
    
    for sd in D_values:
        for md in d_values:
            for f_val in f_values:
                sn, pr, f1, exec_time = evaluate(true_pairs, md, sd, f_val)
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

    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv("evaluation_results.csv", index=False)
    print("\nSaved evaluation_results.csv")

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
        subprocess.run(["./sedef/sedef.sh", "-o", "sedef_out", "-f", "-j", "4", "sim.fa"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        sedef_exec_time = time.time() - t0
        if os.path.exists("sedef_out/final.bed"):
            sedef_sn, sedef_pr, sedef_f1 = evaluate_bedpe(true_pairs, "sedef_out/final.bed")
        print(f"SEDEF -> Recall: {sedef_sn:.4f}, Pr: {sedef_pr:.4f}, F1: {sedef_f1:.4f}, Time: {sedef_exec_time:.4f}s")
    except Exception as e:
        print("SEDEF execution failed:", e)


    # Plotting Sn-Pr Curve
    plt.figure(figsize=(8, 6))
    
    # Sort values by max_dist so that the line plot connects the points in order of max_dist
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
