#!/usr/bin/env python3
"""
cmp_core.py - Core evaluation engine for Segmental Duplication (SD) analysis.

Provides high-speed Base-Pair (BP) footprint calculations, 1D fragment reciprocal overlap (50%),
and binary-search accelerated 2D fragment pair evaluations.
"""

import os
import bisect

def parse_bed_intervals(filepath):
    """Parses any BED/BEDPE file into unique (chrom, start, end) intervals."""
    intervals = []
    if not filepath or not os.path.exists(filepath):
        return intervals
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.strip().split()
            if len(p) >= 3:
                c1 = p[0].split('-')[-1] if '-' in p[0] and not p[0].startswith('chr') else p[0]
                try:
                    intervals.append((c1, int(p[1]), int(p[2])))
                except ValueError:
                    pass
            if len(p) >= 6:
                try:
                    c2 = p[3].split('-')[-1] if '-' in p[3] and not p[3].startswith('chr') else p[3]
                    intervals.append((c2, int(p[4]), int(p[5])))
                except ValueError:
                    pass
    return list(set(intervals))

def merge_intervals_by_chrom(intervals):
    """Merges overlapping intervals per chromosome."""
    by_c = {}
    for c, s, e in intervals:
        by_c.setdefault(c, []).append([s, e])
    merged = {}
    for c, ints in by_c.items():
        ints.sort()
        m = [list(ints[0])]
        for curr in ints[1:]:
            if curr[0] <= m[-1][1]:
                m[-1][1] = max(m[-1][1], curr[1])
            else:
                m.append(list(curr))
        merged[c] = m
    return merged

def calc_bp_metrics(pred_intervals, ref_intervals):
    """Computes Base-Pair (BP) footprint overlap metrics."""
    p_m = merge_intervals_by_chrom(pred_intervals)
    r_m = merge_intervals_by_chrom(ref_intervals)
    p_bp = sum(e - s for ints in p_m.values() for s, e in ints)
    r_bp = sum(e - s for ints in r_m.values() for s, e in ints)
    
    is_bp = 0
    for c in set(p_m.keys()).intersection(r_m.keys()):
        p_list, r_list = p_m[c], r_m[c]
        i = j = 0
        while i < len(p_list) and j < len(r_list):
            overlap = max(0, min(p_list[i][1], r_list[j][1]) - max(p_list[i][0], r_list[j][0]))
            is_bp += overlap
            if p_list[i][1] < r_list[j][1]:
                i += 1
            else:
                j += 1

    rec = is_bp / r_bp if r_bp > 0 else 0.0
    prec = is_bp / p_bp if p_bp > 0 else 0.0
    f1 = 2 * rec * prec / (rec + prec) if (rec + prec) > 0 else 0.0
    union = p_bp + r_bp - is_bp
    return {
        'pred_bp': p_bp, 'ref_bp': r_bp, 'is_bp': is_bp,
        'pred_unique_bp': p_bp - is_bp, 'ref_unique_bp': r_bp - is_bp,
        'recall': rec, 'precision': prec, 'f1': f1,
        'jaccard': is_bp / union if union > 0 else 0.0
    }

def eval_fragment_overlap(pred_intervals, ref_intervals, fraction=0.5):
    """Evaluates 1D fragment intervals using binary search (>= fraction coverage)."""
    p_by_c, r_by_c = {}, {}
    for c, s, e in pred_intervals: p_by_c.setdefault(c, []).append((s, e))
    for c, s, e in ref_intervals: r_by_c.setdefault(c, []).append((s, e))

    p_sorted = {c: sorted(ints, key=lambda x: x[0]) for c, ints in p_by_c.items()}
    p_starts = {c: [x[0] for x in ints] for c, ints in p_sorted.items()}
    r_sorted = {c: sorted(ints, key=lambda x: x[0]) for c, ints in r_by_c.items()}
    r_starts = {c: [x[0] for x in ints] for c, ints in r_sorted.items()}

    tp = 0
    for c, s_r, e_r in ref_intervals:
        l_r = e_r - s_r
        if l_r <= 0 or c not in p_sorted: continue
        idx = bisect.bisect_right(p_starts[c], e_r)
        for k in range(idx - 1, -1, -1):
            s_p, e_p = p_sorted[c][k]
            if e_p <= s_r: continue
            if max(0, min(e_r, e_p) - max(s_r, s_p)) / l_r >= fraction:
                tp += 1
                break

    matched_p = 0
    for c, s_p, e_p in pred_intervals:
        l_p = e_p - s_p
        if l_p <= 0 or c not in r_sorted: continue
        idx = bisect.bisect_right(r_starts[c], e_p)
        for k in range(idx - 1, -1, -1):
            s_r, e_r = r_sorted[c][k]
            if e_r <= s_p: continue
            if max(0, min(e_p, e_r) - max(s_p, s_r)) / l_p >= fraction:
                matched_p += 1
                break

    n_ref, n_pred = len(ref_intervals), len(pred_intervals)
    rec = tp / n_ref if n_ref > 0 else 0.0
    prec = matched_p / n_pred if n_pred > 0 else 0.0
    f1 = 2 * rec * prec / (rec + prec) if (rec + prec) > 0 else 0.0
    return {'recall': rec, 'precision': prec, 'f1': f1, 'tp': tp, 'fp': n_pred - matched_p, 'fn': n_ref - tp, 'total_ref': n_ref, 'total_pred': n_pred}

eval_reciprocal_overlap = eval_fragment_overlap

def load_segtrace_pairs(in_path, use_subclusters=True):
    """Loads paired SDs from Segtrace .dup.bed. Filters same-subcluster orthologs when use_subclusters=True."""
    clusters = {}
    if not in_path or not os.path.exists(in_path): return []
    with open(in_path, 'r') as fin:
        for line in fin:
            if line.startswith('#') or not line.strip(): continue
            p = line.strip().split()
            if len(p) >= 4:
                c = p[0].split('-')[-1] if '-' in p[0] and not p[0].startswith('chr') else p[0]
                clusters.setdefault(p[3], []).append((c, int(p[1]), int(p[2]), p[4] if len(p) >= 5 else "0"))

    pairs = []
    for locs in clusters.values():
        n = len(locs)
        for i in range(n):
            for j in range(i + 1, n):
                c1, s1, e1, sub1 = locs[i]
                c2, s2, e2, sub2 = locs[j]
                if use_subclusters and sub1 == sub2 and sub1 != "0": continue
                if c1 == c2 and max(s1, s2) < min(e1, e2): continue
                pairs.append(((c1, s1, e1), (c2, s2, e2)))
    return pairs

def load_bedpe_pairs(in_path):
    """Loads paired SD regions from BEDPE file (SEDEF, BISER)."""
    pairs = []
    if not in_path or not os.path.exists(in_path): return pairs
    with open(in_path, 'r') as fin:
        for line in fin:
            if line.startswith('#') or not line.strip(): continue
            p = line.strip().split()
            if len(p) >= 6:
                try:
                    c1 = p[0].split('-')[-1] if '-' in p[0] and not p[0].startswith('chr') else p[0]
                    s1, e1 = int(p[1]), int(p[2])
                    s2, e2 = int(p[4]), int(p[5])
                    c2 = p[3].split('-')[-1] if '-' in p[3] and not p[3].startswith('chr') else p[3]
                    pairs.append(((c1, s1, e1), (c2, s2, e2)))
                except (ValueError, IndexError):
                    continue
    return pairs

def evaluate_frag_pairs_fast(ref_pairs, target_pairs, threshold=0.5):
    """Binary-search accelerated 2D fragment pair reciprocal overlap evaluation."""
    n_ref, n_tgt = len(ref_pairs) if ref_pairs else 0, len(target_pairs) if target_pairs else 0
    if n_ref == 0 or n_tgt == 0:
        return {'recall': 0.0, 'precision': 0.0, 'f1': 0.0, 'tp': 0, 'fp': n_tgt, 'fn': n_ref, 'total_ref': n_ref, 'total_pred': n_tgt}

    tgt_by_c = {}
    for t_idx, ((c1, s1, e1), (c2, s2, e2)) in enumerate(target_pairs):
        k = (c1, c2) if c1 <= c2 else (c2, c1)
        tgt_by_c.setdefault(k, []).append(((s1, e1), (s2, e2), t_idx) if c1 <= c2 else ((s2, e2), (s1, e1), t_idx))

    tgt_idx = {k: (sorted(l, key=lambda x: x[0][0]), [x[0][0] for x in sorted(l, key=lambda x: x[0][0])]) for k, l in tgt_by_c.items()}
    m_ref, m_tgt = set(), set()

    for r_idx, ((c1, s1, e1), (c2, s2, e2)) in enumerate(ref_pairs):
        k = (c1, c2) if c1 <= c2 else (c2, c1)
        if k not in tgt_idx: continue
        r1_s, r1_e, r2_s, r2_e = (s1, e1, s2, e2) if c1 <= c2 else (s2, e2, s1, e1)
        l1, l2 = float(r1_e - r1_s), float(r2_e - r2_s)
        if l1 <= 0 or l2 <= 0: continue

        t_list, x_starts = tgt_idx[k]
        for i in range(bisect.bisect_left(x_starts, r1_e)):
            (x1, y1), (x2, y2), orig_idx = t_list[i]
            if y1 <= r1_s: continue
            lp1, lp2 = float(y1 - x1), float(y2 - x2)
            if lp1 <= 0 or lp2 <= 0: continue

            # Direct orientation
            o1, o2 = max(0.0, min(r1_e, y1) - max(r1_s, x1)), max(0.0, min(r2_e, y2) - max(r2_s, x2))
            if o1 / l1 >= threshold and o1 / lp1 >= threshold and o2 / l2 >= threshold and o2 / lp2 >= threshold:
                m_ref.add(r_idx); m_tgt.add(orig_idx); continue

            # Cross orientation
            o12, o21 = max(0.0, min(r1_e, y2) - max(r1_s, x2)), max(0.0, min(r2_e, y1) - max(r2_s, x1))
            if o12 / l1 >= threshold and o12 / lp2 >= threshold and o21 / l2 >= threshold and o21 / lp1 >= threshold:
                m_ref.add(r_idx); m_tgt.add(orig_idx)

    tp, matched_p = len(m_ref), len(m_tgt)
    rec = tp / n_ref if n_ref > 0 else 0.0
    prec = matched_p / n_tgt if n_tgt > 0 else 0.0
    f1 = 2 * rec * prec / (rec + prec) if (rec + prec) > 0 else 0.0
    return {'recall': rec, 'precision': prec, 'f1': f1, 'tp': tp, 'fp': n_tgt - matched_p, 'fn': n_ref - tp, 'total_ref': n_ref, 'total_pred': n_tgt}
