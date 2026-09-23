#!/usr/bin/awk -f
# annotate-failed-validation.awk
#
# Profile SegTrace validation clusters that did NOT reach 100% member recovery.
# It joins the per-cluster summary (cluster_match_ratio.csv) with the per-member
# detail table (cluster_member_recovery.tsv) and, for every cluster whose
# match_ratio < 1, explains WHY each unmatched member was dropped.
#
# Failure categories (as emitted by validate_segtrace.sh classify_member):
#   no_blast_hit    - no BLAST hit >= MIN_IDENT% identity overlaps the member.
#                     Sub-split here into:
#                       * no_hit_anywhere  - no BLAST hit at all for the member.
#                       * hit_off_target   - hits exist but none cover the member
#                                            coordinates (positional mismatch).
#   fail_min_hit_bp - best overlapping hit spans < MIN_HIT_BP bases (short core).
#   fail_min_overlap- best hit covers < MIN_OVERLAP of the shorter of member/hit.
#
# Usage:
#   awk -f annotate-failed-validation.awk \
#       cluster_match_ratio.csv cluster_member_recovery.tsv
#
# Optional thresholds (only affect the human-readable text, must match the run):
#   awk -v MIN_OVERLAP=0.5 -v MIN_HIT_BP=100 -v MIN_IDENT=60 \
#       -f annotate-failed-validation.awk <csv> <tsv>

BEGIN {
    FS = "\t"
    OFS = "\t"
    if (MIN_OVERLAP == "") MIN_OVERLAP = 0.5
    if (MIN_HIT_BP  == "") MIN_HIT_BP  = 100
    if (MIN_IDENT   == "") MIN_IDENT   = 60
}

# --- pass 1: cluster_match_ratio.csv (comma-separated) --------------------
FILENAME ~ /\.csv$/ {
    n = split($0, c, ",")
    if (c[1] == "cluster_id") next          # header
    cid = c[1]
    csv_seen[cid]   = 1
    n_members[cid]  = c[2] + 0
    n_matched[cid]  = c[3] + 0
    ratio[cid]      = c[4] + 0
    qchrom[cid]     = c[5]
    f_bp[cid]       = c[6] + 0
    f_ov[cid]       = c[7] + 0
    f_no[cid]       = c[8] + 0
    next
}

# --- pass 2: cluster_member_recovery.tsv (tab-separated) ------------------
{
    if ($1 == "cluster_id") next            # header
    cid    = $1
    status = $7
    if (status == "matched") next           # only profile the misses

    chrom     = $2
    start     = $3
    end       = $4
    isq       = $5
    mlen      = $6
    pident    = $8
    hitlen    = $9
    ov_bp     = $10
    ov_frac   = $11

    # Refine no_blast_hit: a present pident means hits existed but off target.
    tag = status
    if (status == "no_blast_hit") {
        if (pident == "") tag = "no_hit_anywhere"
        else              tag = "hit_off_target"
    }

    # Track the "closest miss" per category to gauge how marginal the failure is.
    if (status == "fail_min_overlap" && (ov_frac + 0) > best_frac[cid]) {
        best_frac[cid] = ov_frac + 0
    }
    if (status == "fail_min_hit_bp" && (ov_bp + 0) > best_bp[cid]) {
        best_bp[cid] = ov_bp + 0
    }

    desc = chrom ":" start "-" end
    if (isq == 1) desc = desc "*"           # '*' marks the representative member
    desc = desc "(" tag
    if (pident != "") desc = desc " pid=" pident
    if (hitlen != "") desc = desc " hitbp=" hitlen
    if (ov_bp  != "") desc = desc " ovbp=" ov_bp
    if (ov_frac != "") desc = desc " ovfrac=" ov_frac
    desc = desc ")"

    if (cid in miss_detail) miss_detail[cid] = miss_detail[cid] "; " desc
    else                    miss_detail[cid] = desc
}

END {
    print "cluster_id", "query_chrom", "n_members", "n_matched", "match_ratio", \
          "primary_reason", "n_fail_min_hit_bp", "n_fail_min_overlap", \
          "n_no_blast_hit", "diagnosis", "unmatched_members"

    n_failed = 0
    for (cid in csv_seen) {
        if (n_matched[cid] >= n_members[cid]) continue   # fully recovered
        n_failed++

        # Primary reason = dominant failure category (tie-break: no > ov > bp).
        primary = "no_blast_hit"; pmax = f_no[cid]
        if (f_ov[cid] > pmax) { primary = "fail_min_overlap"; pmax = f_ov[cid] }
        if (f_bp[cid] > pmax) { primary = "fail_min_hit_bp";  pmax = f_bp[cid] }
        reason_count[primary]++

        diag = diagnose(primary, cid)
        details = (cid in miss_detail) ? miss_detail[cid] : "(no per-member row)"

        print cid, qchrom[cid], n_members[cid], n_matched[cid], \
              sprintf("%.4f", ratio[cid]), primary, \
              f_bp[cid], f_ov[cid], f_no[cid], diag, details
    }

    # Summary to stderr so it never pollutes the parseable TSV on stdout.
    printf("\n[summary] failed clusters: %d\n", n_failed) > "/dev/stderr"
    printf("  primary no_blast_hit    : %d\n", reason_count["no_blast_hit"])    > "/dev/stderr"
    printf("  primary fail_min_overlap: %d\n", reason_count["fail_min_overlap"]) > "/dev/stderr"
    printf("  primary fail_min_hit_bp : %d\n", reason_count["fail_min_hit_bp"])  > "/dev/stderr"
    printf("  thresholds: MIN_OVERLAP=%s MIN_HIT_BP=%s MIN_IDENT=%s\n", \
           MIN_OVERLAP, MIN_HIT_BP, MIN_IDENT) > "/dev/stderr"
}

function diagnose(primary, cid,   msg) {
    if (primary == "no_blast_hit") {
        msg = "no homology recovered: member has no BLAST hit >= " MIN_IDENT \
              "% identity overlapping its coordinates (novel/diverged segment or off-target hits)"
    } else if (primary == "fail_min_hit_bp") {
        msg = "shared core too short: best overlap < " MIN_HIT_BP " bp"
        if (cid in best_bp) msg = msg " (closest = " best_bp[cid] " bp)"
    } else if (primary == "fail_min_overlap") {
        msg = "partial overlap: best hit covers < " MIN_OVERLAP \
              " of the member/hit (window-boundary or fragmentary homology)"
        if (cid in best_frac) msg = msg " (closest = " sprintf("%.3f", best_frac[cid]) ")"
    } else {
        msg = "unclassified"
    }
    return msg
}
