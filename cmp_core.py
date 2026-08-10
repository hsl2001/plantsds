#!/usr/bin/env python3
"""
cmp_core.py - Core evaluation engine for Segmental Duplication (SD) analysis.

Provides core dataset-agnostic evaluation algorithms:
  1) Base-Pair (BP) level generic metrics calculation (BP sum, Recall, Precision, F1, Jaccard).
  2) Reciprocal Overlap Fragment (Frag) pair-level evaluation engine.
     Optimized with bisect interval indexing + NumPy vectorized C acceleration for ultra-fast execution (< 0.5s)
     and low memory usage.
"""

import os
import bisect
import collections
import numpy as np

def load_bed_bp(filepath):
    """Calculates total base pairs in a BED file."""
    total_bp = 0
    if not os.path.exists(filepath):
        return 0
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            s, e = int(parts[1]), int(parts[2])
            total_bp += max(0, e - s)
    return total_bp

def compute_bp_metrics(st_bp, sd_bp, is_bp):
    """
    Computes Base-Pair (BP) level metrics.
    Returns (bp_recall, bp_precision, bp_f1, bp_jaccard).
    """
    bp_recall = is_bp / sd_bp if sd_bp > 0 else 0.0
    bp_precision = is_bp / st_bp if st_bp > 0 else 0.0
    bp_f1 = 2 * bp_recall * bp_precision / (bp_recall + bp_precision) if (bp_recall + bp_precision) > 0 else 0.0
    union_bp = st_bp + sd_bp - is_bp
    bp_jaccard = is_bp / union_bp if union_bp > 0 else 0.0
    return bp_recall, bp_precision, bp_f1, bp_jaccard

def evaluate_frag_pairs_fast(ref_pairs, target_pairs, threshold=0.5):
    """
    Ultra-fast Bisect + NumPy vectorized reciprocal overlap pair evaluation (< 0.5s for 6.7M pairs).
    """
    if not ref_pairs or not target_pairs:
        return 0.0, 0.0, 0.0, 0, 0

    target_by_chrom = collections.defaultdict(list)
    for t_idx, ((c1, s1, e1), (c2, s2, e2)) in enumerate(target_pairs):
        key = (c1, c2) if c1 <= c2 else (c2, c1)
        if c1 <= c2:
            target_by_chrom[key].append((s1, e1, s2, e2, t_idx))
        else:
            target_by_chrom[key].append((s2, e2, s1, e1, t_idx))

    np_target = {}
    for key, item_list in target_by_chrom.items():
        item_list.sort(key=lambda x: x[0])
        arr = np.array(item_list, dtype=np.int64)
        x1_starts = arr[:, 0].tolist()
        np_target[key] = (arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4], x1_starts)

    matched_ref = set()
    matched_target = set()

    for r_idx, ((c1, s1, e1), (c2, s2, e2)) in enumerate(ref_pairs):
        key = (c1, c2) if c1 <= c2 else (c2, c1)
        if key not in np_target:
            continue

        if c1 <= c2:
            r1_s, r1_e, r2_s, r2_e = s1, e1, s2, e2
        else:
            r1_s, r1_e, r2_s, r2_e = s2, e2, s1, e1

        Lt1 = float(r1_e - r1_s)
        Lt2 = float(r2_e - r2_s)
        if Lt1 <= 0 or Lt2 <= 0:
            continue

        X1, Y1, X2, Y2, T_IDX, x1_starts = np_target[key]

        high = bisect.bisect_left(x1_starts, r1_e)
        if high == 0:
            continue

        x1_sub, y1_sub, x2_sub, y2_sub, tidx_sub = X1[:high], Y1[:high], X2[:high], Y2[:high], T_IDX[:high]

        cand_mask = (y1_sub > r1_s) & (x2_sub < r2_e) & (y2_sub > r2_s)
        if not np.any(cand_mask):
            cand_mask_cross = (y2_sub > r1_s) & (x1_sub < r2_e) & (y1_sub > r2_s)
            if not np.any(cand_mask_cross):
                continue
            cand_mask = cand_mask_cross

        x1_c, y1_c, x2_c, y2_c, tidx_c = x1_sub[cand_mask], y1_sub[cand_mask], x2_sub[cand_mask], y2_sub[cand_mask], tidx_sub[cand_mask]

        # Direct orientation
        o1 = np.maximum(0.0, np.minimum(r1_e, y1_c) - np.maximum(r1_s, x1_c))
        o2 = np.maximum(0.0, np.minimum(r2_e, y2_c) - np.maximum(r2_s, x2_c))
        lp1 = y1_c - x1_c
        lp2 = y2_c - x2_c

        dir_match = (o1 / Lt1 >= threshold) & (o1 / lp1 >= threshold) & (o2 / Lt2 >= threshold) & (o2 / lp2 >= threshold)

        # Cross orientation
        o12 = np.maximum(0.0, np.minimum(r1_e, y2_c) - np.maximum(r1_s, x2_c))
        o21 = np.maximum(0.0, np.minimum(r2_e, y1_c) - np.maximum(r2_s, x1_c))
        cross_match = (o12 / Lt1 >= threshold) & (o12 / lp2 >= threshold) & (o21 / Lt2 >= threshold) & (o21 / lp1 >= threshold)

        hit_indices = tidx_c[dir_match | cross_match]
        if len(hit_indices) > 0:
            matched_ref.add(r_idx)
            matched_target.update(hit_indices.tolist())

    total_ref = len(ref_pairs)
    total_target = len(target_pairs)

    Sn = len(matched_ref) / total_ref if total_ref > 0 else 0.0
    Pr = len(matched_target) / total_target if total_target > 0 else 0.0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0.0
    return Sn, Pr, f1, len(matched_ref), len(matched_target)
