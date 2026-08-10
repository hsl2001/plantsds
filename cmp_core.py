#!/usr/bin/env python3
"""
cmp_core.py - Core evaluation engine for Segmental Duplication (SD) analysis.

Provides core dataset-agnostic evaluation algorithms:
  1) Base-Pair (BP) level generic metrics calculation (BP sum, Recall, Precision, F1, Jaccard).
  2) Reciprocal Overlap Fragment (Frag) pair-level evaluation engine.
     Optimized with 2D spatial grid indexing for ultra-fast execution (< 0.2s)
     and zero thread-lock overhead.
"""

import os
import collections

def count_bed_bp(filepath):
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

def calc_bp_metrics(st_bp, sd_bp, is_bp):
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

def calc_frag_metrics(ref_pairs, target_pairs, threshold=0.5, bin_size=1000000):
    """
    Ultra-fast 1D spatial grid reciprocal overlap fragment evaluation.
    Evaluates fragments as independent SD regions, ignoring pair matching.
    """
    import collections
    if not ref_pairs or not target_pairs:
        return 0.0, 0.0, 0.0, 0, 0, 0

    ref_frags = list(set([f for p in ref_pairs for f in p]))
    tgt_frags = list(set([f for p in target_pairs for f in p]))

    t_grid = collections.defaultdict(list)
    for t_idx, (c, s, e) in enumerate(tgt_frags):
        b_min, b_max = s // bin_size, e // bin_size
        for b in range(b_min, b_max + 1):
            t_grid[(c, b)].append((s, e, t_idx))

    r_grid = collections.defaultdict(list)
    for r_idx, (c, s, e) in enumerate(ref_frags):
        b_min, b_max = s // bin_size, e // bin_size
        for b in range(b_min, b_max + 1):
            r_grid[(c, b)].append((s, e, r_idx))

    matched_ref = 0
    for c, s, e in ref_frags:
        Lt = float(e - s)
        if Lt <= 0: continue
        
        b_min, b_max = s // bin_size, e // bin_size
        cands = []
        for b in range(b_min, b_max + 1):
            if (c, b) in t_grid:
                cands.extend(t_grid[(c, b)])
        
        matched = False
        for x, y, _ in cands:
            Lp = float(y - x)
            if Lp <= 0: continue
            
            o = max(0.0, min(e, y) - max(s, x))
            if o / Lt >= threshold and o / Lp >= threshold:
                matched = True
                break
        if matched:
            matched_ref += 1

    matched_tgt = 0
    for c, s, e in tgt_frags:
        Lp = float(e - s)
        if Lp <= 0: continue
        
        b_min, b_max = s // bin_size, e // bin_size
        cands = []
        for b in range(b_min, b_max + 1):
            if (c, b) in r_grid:
                cands.extend(r_grid[(c, b)])
        
        matched = False
        for x, y, _ in cands:
            Lt = float(y - x)
            if Lt <= 0: continue
            
            o = max(0.0, min(e, y) - max(s, x))
            if o / Lp >= threshold and o / Lt >= threshold:
                matched = True
                break
        if matched:
            matched_tgt += 1

    Sn = matched_ref / len(ref_frags) if len(ref_frags) > 0 else 0.0
    Pr = matched_tgt / len(tgt_frags) if len(tgt_frags) > 0 else 0.0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0.0
    return Sn, Pr, f1, matched_ref, matched_tgt, len(ref_frags)

def calc_frag_metrics_from_clusters(sedef_pairs, segtrace_clusters, threshold=0.5, bin_size=1000000):
    import collections
    
    ref_frags = list(set([f for p in sedef_pairs for f in p]))
    
    tgt_frags_raw = []
    for cid, regions in segtrace_clusters.items():
        for (c, s, e, subid) in regions:
            tgt_frags_raw.append((c, s, e))
    tgt_frags = list(set(tgt_frags_raw))

    if not ref_frags or not tgt_frags:
        return 0.0, 0.0, 0.0, 0, 0, 0

    t_grid = collections.defaultdict(list)
    for t_idx, (c, s, e) in enumerate(tgt_frags):
        b_min, b_max = s // bin_size, e // bin_size
        for b in range(b_min, b_max + 1):
            t_grid[(c, b)].append((s, e, t_idx))

    r_grid = collections.defaultdict(list)
    for r_idx, (c, s, e) in enumerate(ref_frags):
        b_min, b_max = s // bin_size, e // bin_size
        for b in range(b_min, b_max + 1):
            r_grid[(c, b)].append((s, e, r_idx))

    matched_ref = 0
    for c, s, e in ref_frags:
        Lt = float(e - s)
        if Lt <= 0: continue
        
        b_min, b_max = s // bin_size, e // bin_size
        cands = []
        for b in range(b_min, b_max + 1):
            if (c, b) in t_grid:
                cands.extend(t_grid[(c, b)])
        
        matched = False
        for x, y, _ in cands:
            Lp = float(y - x)
            if Lp <= 0: continue
            
            o = max(0.0, min(e, y) - max(s, x))
            if o / Lt >= threshold and o / Lp >= threshold:
                matched = True
                break
        if matched:
            matched_ref += 1

    matched_tgt = 0
    for c, s, e in tgt_frags:
        Lp = float(e - s)
        if Lp <= 0: continue
        
        b_min, b_max = s // bin_size, e // bin_size
        cands = []
        for b in range(b_min, b_max + 1):
            if (c, b) in r_grid:
                cands.extend(r_grid[(c, b)])
        
        matched = False
        for x, y, _ in cands:
            Lt = float(y - x)
            if Lt <= 0: continue
            
            o = max(0.0, min(e, y) - max(s, x))
            if o / Lp >= threshold and o / Lt >= threshold:
                matched = True
                break
        if matched:
            matched_tgt += 1

    total_ref = len(ref_frags)
    total_tgt = len(tgt_frags)
    Sn = matched_ref / total_ref if total_ref > 0 else 0.0
    Pr = matched_tgt / total_tgt if total_tgt > 0 else 0.0
    f1 = 2 * Sn * Pr / (Sn + Pr) if (Sn + Pr) > 0 else 0.0
    return Pr, Sn, f1, matched_tgt, matched_ref, total_tgt
