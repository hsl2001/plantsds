#!/usr/bin/env python3
"""
cmp_core.py - Core evaluation engine for Segmental Duplication (SD) analysis.

Provides high-speed base-pair footprint calculations and 1D fragment reciprocal overlap (50%).
"""

import os
import bisect


def _to_int_or_none(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None

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
                s1 = _to_int_or_none(p[1])
                e1 = _to_int_or_none(p[2])
                if s1 is not None and e1 is not None and e1 > s1:
                    intervals.append((c1, s1, e1))

            # Parse second interval only when line is BEDPE-like:
            # fields 4-6 must be (chrom, start, end) with valid coordinates.
            if len(p) >= 6:
                s2 = _to_int_or_none(p[4])
                e2 = _to_int_or_none(p[5])
                if s2 is not None and e2 is not None and e2 > s2:
                    c2 = p[3].split('-')[-1] if '-' in p[3] and not p[3].startswith('chr') else p[3]
                    intervals.append((c2, s2, e2))
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
    """Evaluates 1D fragment intervals using reciprocal overlap (>= fraction on both sides)."""
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
            ov = max(0, min(e_r, e_p) - max(s_r, s_p))
            l_p = e_p - s_p
            if l_p > 0 and (ov / l_r) >= fraction and (ov / l_p) >= fraction:
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
            ov = max(0, min(e_p, e_r) - max(s_p, s_r))
            l_r = e_r - s_r
            if l_r > 0 and (ov / l_p) >= fraction and (ov / l_r) >= fraction:
                matched_p += 1
                break

    n_ref, n_pred = len(ref_intervals), len(pred_intervals)
    rec = tp / n_ref if n_ref > 0 else 0.0
    prec = matched_p / n_pred if n_pred > 0 else 0.0
    f1 = 2 * rec * prec / (rec + prec) if (rec + prec) > 0 else 0.0
    return {'recall': rec, 'precision': prec, 'f1': f1, 'tp': tp, 'fp': n_pred - matched_p, 'fn': n_ref - tp, 'total_ref': n_ref, 'total_pred': n_pred}

eval_reciprocal_overlap = eval_fragment_overlap
