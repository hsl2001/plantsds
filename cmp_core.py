#!/usr/bin/env python3
"""
cmp_core.py - Core evaluation engine for Segmental Duplication (SD) analysis.

and Base-Pair (BP) footprint calculations across datasets.
"""

import os

def parse_bed_intervals(filepath):
    """Parses any BED file into a list of unique (chrom, start, end) intervals."""
    intervals = []
    if not filepath or not os.path.exists(filepath):
        return intervals
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                c1 = parts[0].split('-')[-1] if '-' in parts[0] and not parts[0].startswith('chr') else parts[0]
                intervals.append((c1, int(parts[1]), int(parts[2])))
            if len(parts) >= 6:  # SEDEF / BEDPE pair second arm (chrom2, start2, end2)
                try:
                    c2 = parts[3].split('-')[-1] if '-' in parts[3] and not parts[3].startswith('chr') else parts[3]
                    intervals.append((c2, int(parts[4]), int(parts[5])))
                except ValueError:
                    pass
    return list(set(intervals))

def merge_intervals_by_chrom(intervals):
    """Merges overlapping intervals per chromosome for footprint calculations."""
    by_chrom = {}
    for c, s, e in intervals:
        by_chrom.setdefault(c, []).append([s, e])
    
    merged = {}
    for c, ints in by_chrom.items():
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
    """Computes Base-Pair (BP) level footprint overlap metrics."""
    p_m = merge_intervals_by_chrom(pred_intervals)
    r_m = merge_intervals_by_chrom(ref_intervals)

    p_bp = sum(e - s for ints in p_m.values() for s, e in ints)
    r_bp = sum(e - s for ints in r_m.values() for s, e in ints)
    
    is_bp = 0
    common_chroms = set(p_m.keys()).intersection(r_m.keys())
    for c in common_chroms:
        p_list, r_list = p_m[c], r_m[c]
        i = j = 0
        while i < len(p_list) and j < len(r_list):
            overlap = max(0, min(p_list[i][1], r_list[j][1]) - max(p_list[i][0], r_list[j][0]))
            is_bp += overlap
            if p_list[i][1] < r_list[j][1]:
                i += 1
            else:
                j += 1

    bp_rec = is_bp / r_bp if r_bp > 0 else 0.0
    bp_prec = is_bp / p_bp if p_bp > 0 else 0.0
    bp_f1 = 2 * bp_rec * bp_prec / (bp_rec + bp_prec) if (bp_rec + bp_prec) > 0 else 0.0
    union_bp = p_bp + r_bp - is_bp
    bp_jaccard = is_bp / union_bp if union_bp > 0 else 0.0

    return {
        'pred_bp': p_bp,
        'ref_bp': r_bp,
        'is_bp': is_bp,
        'pred_unique_bp': p_bp - is_bp,
        'ref_unique_bp': r_bp - is_bp,
        'recall': bp_rec,
        'precision': bp_prec,
        'f1': bp_f1,
        'jaccard': bp_jaccard
    }

def eval_fragment_overlap(pred_intervals, ref_intervals, fraction=0.5):
    """
    Evaluates fragment intervals:
    - Recall (Sensitivity): Fraction of reference intervals covered >= fraction (ov / len_r >= fraction).
    - Precision: Fraction of predicted intervals valid SD sequence (ov / len_p >= fraction).
    Calculates TP, FP, FN, Recall, Precision, and F1.
    """
    p_by_c = {}
    for c, s, e in pred_intervals:
        p_by_c.setdefault(c, []).append((s, e))
        
    r_by_c = {}
    for c, s, e in ref_intervals:
        r_by_c.setdefault(c, []).append((s, e))

    # TP & FN: Reference Intervals covered >= fraction by any prediction
    tp = 0
    total_ref = len(ref_intervals)
    for c, s_r, e_r in ref_intervals:
        len_r = e_r - s_r
        if len_r <= 0 or c not in p_by_c:
            continue
        matched = False
        for s_p, e_p in p_by_c[c]:
            ov = max(0, min(e_r, e_p) - max(s_r, s_p))
            if ov / len_r >= fraction:
                matched = True
                break
        if matched:
            tp += 1

    fn = total_ref - tp

    # Matched Predictions & FP: Predictions where >= fraction of prediction is valid reference SD
    matched_pred = 0
    total_pred = len(pred_intervals)
    for c, s_p, e_p in pred_intervals:
        len_p = e_p - s_p
        if len_p <= 0 or c not in r_by_c:
            continue
        matched = False
        for s_r, e_r in r_by_c[c]:
            ov = max(0, min(e_r, e_p) - max(s_r, s_p))
            if ov / len_p >= fraction:
                matched = True
                break
        if matched:
            matched_pred += 1

    fp = total_pred - matched_pred

    recall = tp / total_ref if total_ref > 0 else 0.0
    precision = matched_pred / total_pred if total_pred > 0 else 0.0
    f1 = 2 * recall * precision / (recall + precision) if (recall + precision) > 0 else 0.0

    return {
        'recall': recall,
        'precision': precision,
        'f1': f1,
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'total_ref': total_ref,
        'total_pred': total_pred
    }

# Backwards compatibility alias
eval_reciprocal_overlap = eval_fragment_overlap
