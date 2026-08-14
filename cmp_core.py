#!/usr/bin/env python3
"""
cmp_core.py - Core evaluation engine for Segmental Duplication (SD) analysis.

Provides high-speed binary-search accelerated fragment overlap evaluation (TP, FP, FN, Recall, Precision, F1)
and Base-Pair (BP) footprint calculations across datasets.
"""

import os
import bisect

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
    Evaluates fragment intervals using binary-search accelerated indexing (O(N log N)):
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

    p_sorted = {}
    p_starts = {}
    for c, ints in p_by_c.items():
        sorted_ints = sorted(ints, key=lambda x: x[0])
        p_sorted[c] = sorted_ints
        p_starts[c] = [x[0] for x in sorted_ints]

    r_sorted = {}
    r_starts = {}
    for c, ints in r_by_c.items():
        sorted_ints = sorted(ints, key=lambda x: x[0])
        r_sorted[c] = sorted_ints
        r_starts[c] = [x[0] for x in sorted_ints]

    # TP & FN: Reference Intervals covered >= fraction by any prediction
    tp = 0
    total_ref = len(ref_intervals)
    for c, s_r, e_r in ref_intervals:
        len_r = e_r - s_r
        if len_r <= 0 or c not in p_sorted:
            continue
        matched = False
        idx_end = bisect.bisect_right(p_starts[c], e_r)
        for k in range(idx_end - 1, -1, -1):
            s_p, e_p = p_sorted[c][k]
            if e_p <= s_r:
                continue
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
        if len_p <= 0 or c not in r_sorted:
            continue
        matched = False
        idx_end = bisect.bisect_right(r_starts[c], e_p)
        for k in range(idx_end - 1, -1, -1):
            s_r, e_r = r_sorted[c][k]
            if e_r <= s_p:
                continue
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

def load_segtrace_pairs(in_path, use_subclusters=True):
    """
    Loads paired SD regions from Segtrace .dup.bed file based on cluster_id and subcluster_id.
    If use_subclusters is True, pairs sharing the same non-zero subcluster_id (same insertion locus)
    are filtered out.
    """
    clusters = {}
    if not in_path or not os.path.exists(in_path):
        return []
    with open(in_path, 'r') as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 4:
                chrom = parts[0].split('-')[-1] if '-' in parts[0] and not parts[0].startswith('chr') else parts[0]
                start, end, cid = int(parts[1]), int(parts[2]), parts[3]
                subid = parts[4] if len(parts) >= 5 else "0"
                clusters.setdefault(cid, []).append((chrom, start, end, subid))

    pairs = []
    for cid, locs in clusters.items():
        n = len(locs)
        for i in range(n):
            for j in range(i + 1, n):
                c1, s1, e1, sub1 = locs[i]
                c2, s2, e2, sub2 = locs[j]
                if use_subclusters and sub1 == sub2 and sub1 != "0":
                    continue
                if c1 == c2 and max(s1, s2) < min(e1, e2):
                    continue
                pairs.append(((c1, s1, e1), (c2, s2, e2)))
    return pairs

def load_bedpe_pairs(in_path):
    """Loads paired SD regions from BEDPE file (SEDEF, BISER, etc.)."""
    pairs = []
    if not in_path or not os.path.exists(in_path):
        return pairs
    with open(in_path, 'r') as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    c1 = parts[0].split('-')[-1] if '-' in parts[0] and not parts[0].startswith('chr') else parts[0]
                    s1, e1 = int(parts[1]), int(parts[2])
                    
                    # Standard BEDPE: cols 0-2 (chr1, start1, end1), cols 3-5 (chr2, start2, end2)
                    try:
                        s2, e2 = int(parts[4]), int(parts[5])
                        c2 = parts[3].split('-')[-1] if '-' in parts[3] and not parts[3].startswith('chr') else parts[3]
                    except ValueError:
                        # Fallback for other formats like chr:start-end in parts[3] or 12-col UCSC format
                        if ':' in parts[3] and '-' in parts[3]:
                            c2_raw, span = parts[3].split(':', 1)
                            c2 = c2_raw.split('-')[-1] if '-' in c2_raw and not c2_raw.startswith('chr') else c2_raw
                            s2_str, e2_str = span.split('-', 1)
                            s2, e2 = int(s2_str), int(e2_str)
                        elif len(parts) >= 12:
                            c2 = parts[9].split('-')[-1] if '-' in parts[9] and not parts[9].startswith('chr') else parts[9]
                            s2, e2 = int(parts[10]), int(parts[11])
                        else:
                            continue
                    pairs.append(((c1, s1, e1), (c2, s2, e2)))
                except (ValueError, IndexError):
                    continue
    return pairs

def evaluate_frag_pairs_fast(ref_pairs, target_pairs, threshold=0.5):
    """
    High-speed binary-search accelerated pairwise fragment overlap evaluation (O(N log M)).
    Evaluates reciprocal overlap (>= threshold) on both arms in direct and cross orientations.
    """
    if not ref_pairs or not target_pairs:
        total_ref = len(ref_pairs) if ref_pairs else 0
        total_tgt = len(target_pairs) if target_pairs else 0
        return {
            'recall': 0.0,
            'precision': 0.0,
            'f1': 0.0,
            'tp': 0,
            'fp': total_tgt,
            'fn': total_ref,
            'total_ref': total_ref,
            'total_pred': total_tgt
        }

    target_by_chrom = {}
    for t_idx, ((c1, s1, e1), (c2, s2, e2)) in enumerate(target_pairs):
        key = (c1, c2) if c1 <= c2 else (c2, c1)
        if key not in target_by_chrom:
            target_by_chrom[key] = []
        if c1 <= c2:
            target_by_chrom[key].append(((s1, e1), (s2, e2), t_idx))
        else:
            target_by_chrom[key].append(((s2, e2), (s1, e1), t_idx))

    target_index = {}
    for key, t_list in target_by_chrom.items():
        t_list.sort(key=lambda item: item[0][0])
        x1_starts = [item[0][0] for item in t_list]
        target_index[key] = (t_list, x1_starts)

    matched_ref = set()
    matched_target = set()

    for r_idx, ((c1, s1, e1), (c2, s2, e2)) in enumerate(ref_pairs):
        key = (c1, c2) if c1 <= c2 else (c2, c1)
        if key not in target_index:
            continue

        if c1 <= c2:
            r1_s, r1_e, r2_s, r2_e = s1, e1, s2, e2
        else:
            r1_s, r1_e, r2_s, r2_e = s2, e2, s1, e1

        Lt1 = float(r1_e - r1_s)
        Lt2 = float(r2_e - r2_s)
        if Lt1 <= 0 or Lt2 <= 0:
            continue

        t_list, x1_starts = target_index[key]
        limit_idx = bisect.bisect_left(x1_starts, r1_e)

        for i in range(limit_idx):
            (x1, y1), (x2, y2), t_orig_idx = t_list[i]
            if y1 <= r1_s:
                continue

            Lp1 = float(y1 - x1)
            Lp2 = float(y2 - x2)
            if Lp1 <= 0 or Lp2 <= 0:
                continue

            # Direct orientation
            o1 = max(0.0, min(r1_e, y1) - max(r1_s, x1))
            o2 = max(0.0, min(r2_e, y2) - max(r2_s, x2))

            if (o1 / Lt1 >= threshold and o1 / Lp1 >= threshold and
                o2 / Lt2 >= threshold and o2 / Lp2 >= threshold):
                matched_ref.add(r_idx)
                matched_target.add(t_orig_idx)
                continue

            # Cross orientation (for intra-chromosomal or symmetric)
            o12 = max(0.0, min(r1_e, y2) - max(r1_s, x2))
            o21 = max(0.0, min(r2_e, y1) - max(r2_s, x1))

            if (o12 / Lt1 >= threshold and o12 / Lp2 >= threshold and
                o21 / Lt2 >= threshold and o21 / Lp1 >= threshold):
                matched_ref.add(r_idx)
                matched_target.add(t_orig_idx)

    total_ref = len(ref_pairs)
    total_target = len(target_pairs)

    tp = len(matched_ref)
    fn = total_ref - tp
    matched_pred = len(matched_target)
    fp = total_target - matched_pred

    recall = tp / total_ref if total_ref > 0 else 0.0
    precision = matched_pred / total_target if total_target > 0 else 0.0
    f1 = 2 * recall * precision / (recall + precision) if (recall + precision) > 0 else 0.0

    return {
        'recall': recall,
        'precision': precision,
        'f1': f1,
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'total_ref': total_ref,
        'total_pred': total_target
    }

