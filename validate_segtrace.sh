#!/usr/bin/env bash
# SegTrace cluster validation with BLAST.
set -euo pipefail

if [ -z "${FASTAS+x}" ]; then
  FASTAS=(tair12.fasta t2t_nip.fasta)
fi
OUTDIR="${OUTDIR:-segtrace_validation}"
PREFIX="${PREFIX:-$OUTDIR/segtrace}"
THREADS="${THREADS:-64}"
MIN_OVERLAP="${MIN_OVERLAP:-0.5}"   # member covered fraction to count as a match
MIN_IDENT="${MIN_IDENT:-90}"        # BLAST identity floor; 
MIN_MEMBER_LEN=1024                  # minimum SegTrace member length in the recovery population
SEGTRACE_EXTRA="${SEGTRACE_EXTRA:-}"
N_CLUSTERS="${N_CLUSTERS:-1000}"    # analyze only the N clusters with the shortest longest-member

BED="$PREFIX.seg.bed"
COMBINED="$OUTDIR/combined.fa"
QUERY="$OUTDIR/cluster_reps.fa"
MEMBERS="$OUTDIR/cluster_members.tsv"
DB="$OUTDIR/blastdb/combined"
BLAST="$OUTDIR/blast_hits.tsv"
CSV="$OUTDIR/cluster_match_ratio.csv"
MEMBER_RESULTS="$OUTDIR/cluster_member_recovery.tsv"
PLOT="$OUTDIR/cluster_match_ratio.png"

mkdir -p "$OUTDIR" "$OUTDIR/blastdb"

# ------------------------------------------------------------------ 1. SegTrace
echo "[1/5] Building and running segtrace..."
make -s segtrace
./segtrace -p "$THREADS" $SEGTRACE_EXTRA -o "$PREFIX" "${FASTAS[@]}"
[ -s "$BED" ] || { echo "[ERROR] no segtrace output: $BED" >&2; exit 1; }
echo "       clusters: $(awk 'NR>1{print $4}' "$BED" | sort -u | wc -l | tr -d ' '), regions: $(($(wc -l < "$BED") - 1))"

# ---------------------------------- 2. combined FASTA (genome-seq ids) + queries
echo "[2/5] Building combined FASTA and per-cluster representative queries..."
python3 - "$COMBINED" "$QUERY" "$MEMBERS" "$BED" "$N_CLUSTERS" "$MIN_MEMBER_LEN" "${FASTAS[@]}" <<'PY'
import os
import sys

combined_path, query_path, members_path, bed_path, top_n_s, min_member_len_s = sys.argv[1:7]
fastas = sys.argv[7:]
top_n = int(top_n_s)
min_member_len = int(min_member_len_s)
W = 60  # FASTA line width; fixed so we can seek into records by base offset


def genome_label(path):
    """Reproduce segtrace get_basename(): strip dir and known FASTA extensions."""
    name = os.path.basename(path)
    for ext in (".gz", ".bgz"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    for ext in (".fa", ".fna", ".fasta", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return name


index = {}  # label -> (byte offset of sequence, base length)


def write_record(out, label, seq):
    out.write(b">" + label.encode() + b"\n")
    offset = out.tell()
    for i in range(0, len(seq), W):
        out.write(seq[i : i + W])
        out.write(b"\n")
    index[label] = (offset, len(seq))


with open(combined_path, "wb") as out:
    for fa in fastas:
        g = genome_label(fa)
        cur = None
        chunks = []
        with open(fa, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    if cur is not None:
                        write_record(out, cur, b"".join(chunks).upper())
                    cur = f"{g}-{line[1:].split()[0].decode()}"
                    chunks = []
                else:
                    chunks.append(line.strip())
            if cur is not None:
                write_record(out, cur, b"".join(chunks).upper())


def extract(fh, offset, length, a, b):
    a = max(0, a)
    b = min(length, b)
    if b <= a:
        return b""
    start_byte = offset + (a // W) * (W + 1) + (a % W)
    last = b - 1
    end_byte = offset + (last // W) * (W + 1) + (last % W) + 1
    fh.seek(start_byte)
    return fh.read(end_byte - start_byte).replace(b"\n", b"")


clusters = {}
total_members = 0
excluded_members = 0
with open(bed_path) as bh:
    for line in bh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        total_members += 1
        member = (f[0], int(f[1]), int(f[2]))
        if member[2] - member[1] < min_member_len:
            excluded_members += 1
            continue
        clusters.setdefault(int(f[3]), []).append(member)

# Keep only clusters with at least two eligible members; a single-member
# cluster would be trivially recovered by its self-hit. Then keep the N
# clusters whose longest eligible member is the shortest.
eligible_cluster_count = len(clusters)
single_eligible_clusters = sum(len(members) < 2 for members in clusters.values())
length_eligible_members = sum(len(members) for members in clusters.values())
clusters = {cid: members for cid, members in clusters.items() if len(members) >= 2}

rep = {}  # cid -> (index of longest member, its length)
for cid, members in clusters.items():
    qi = max(range(len(members)), key=lambda i: members[i][2] - members[i][1])
    rep[cid] = (qi, members[qi][2] - members[qi][1])
selected = sorted(clusters, key=lambda cid: (rep[cid][1], cid))[:top_n]

with open(combined_path, "rb") as cf, open(query_path, "wb") as qf, \
        open(members_path, "w") as mf:
    mf.write("cluster_id\tchrom\tstart\tend\tis_query\n")
    for cid in sorted(selected):
        members = clusters[cid]
        qi = rep[cid][0]
        qchrom, qs, qe = members[qi]
        off, length = index[qchrom]
        seq = extract(cf, off, length, qs, qe)
        if seq:
            qf.write(f">c{cid}\n".encode())
            for i in range(0, len(seq), W):
                qf.write(seq[i : i + W])
                qf.write(b"\n")
        for i, (c, s, e) in enumerate(members):
            mf.write(f"{cid}\t{c}\t{s}\t{e}\t{1 if i == qi else 0}\n")

selected_members = sum(len(clusters[cid]) for cid in selected)
print(f"       members >= {min_member_len} bp: {length_eligible_members}/{total_members} (excluded {excluded_members} shorter members)")
print(f"       clusters with >=2 eligible members: {len(clusters)}/{eligible_cluster_count} (excluded {single_eligible_clusters} single-member clusters)")
print(f"       recovery population: {selected_members} members across {len(selected)} selected clusters")
PY

# --------------------------------------------------------- 3. BLAST DB + search
if [ -s "$QUERY" ]; then
    echo "[3/5] Building BLAST database..."
    makeblastdb -dbtype nucl -in "$COMBINED" -out "$DB" >/dev/null

    echo "[4/5] BLASTing eligible cluster representatives against the genomes..."
    # dc-megablast (discontiguous seeds) finds diverged homology that plain megablast
    # (word size 28) misses entirely for <~90% identity duplications.
    blastn -task dc-megablast -query "$QUERY" -db "$DB" -num_threads "$THREADS" \
        -perc_identity "$MIN_IDENT" -evalue 1e-5 -max_target_seqs 100000 \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
        -out "$BLAST"
else
    echo "[3/5] No eligible cluster members; skipping BLAST."
    : > "$BLAST"
fi

# ------------------------------------------------ 5. match ratio + CSV + graph
echo "[5/5] Scoring cluster recovery and plotting..."
python3 - "$MEMBERS" "$BLAST" "$CSV" "$MEMBER_RESULTS" "$PLOT" "$MIN_OVERLAP" "$MIN_IDENT" <<'PY'
import csv
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

members_path, blast_path, csv_path, member_results_path, plot_path, min_overlap_s, min_ident_s = sys.argv[1:8]
min_overlap = float(min_overlap_s)
min_ident = float(min_ident_s)

members = defaultdict(list)  # cid -> [(chrom, start, end, is_query)]
with open(members_path) as mh:
    next(mh)
    for line in mh:
        cid, chrom, s, e, isq = line.rstrip("\n").split("\t")
        members[int(cid)].append((chrom, int(s), int(e), int(isq)))

hits = defaultdict(list)  # (cid, chrom) -> [(start, end, pident, hsp_length)] on subject
with open(blast_path) as bh:
    for line in bh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or float(f[2]) < min_ident:
            continue
        sstart, send = int(f[6]), int(f[7])
        hits[(int(f[0][1:]), f[1])].append((
            min(sstart, send) - 1,
            max(sstart, send),
            float(f[2]),
            int(f[3]),
        ))


def classify_member(start, end, intervals):
    mlen = end - start
    if mlen <= 0:
        return "no_blast_hit", None

    overlapping = [
        (overlap, overlap / min(mlen, b - a), pident, length, a, b)
        for a, b, pident, length in intervals
        if (overlap := max(0, min(end, b) - max(start, a))) > 0
    ]
    if not overlapping:
        return "no_blast_hit", max(intervals, key=lambda h: (h[2], h[3]), default=None)

    key = lambda h: (h[0], h[1], h[2], h[3])
    valid = [h for h in overlapping if h[1] >= min_overlap]
    if valid:
        return "matched", max(valid, key=key)
    best = max(overlapping, key=key)
    return "fail_min_overlap", best


rows = []
unmatched = Counter()
with open(member_results_path, "w", newline="") as mr:
    writer = csv.writer(mr, delimiter="\t")
    writer.writerow(["cluster_id", "chrom", "start", "end", "is_query", "member_length", "status",
                     "best_pident", "best_hit_length", "overlap_bp", "overlap_fraction",
                     "best_subject_start", "best_subject_end"])
    for cid, mem in sorted(members.items()):
        counts = Counter()
        for chrom, start, end, is_query in mem:
            status, best = classify_member(start, end, hits.get((cid, chrom), []))
            counts[status] += 1
            if status != "matched":
                unmatched[status] += 1
            if best is None:
                details = ("", "", "", "", "", "")
            elif len(best) == 4:
                a, b, pident, length = best
                details = (pident, length, 0, 0.0, a, b)
            else:
                overlap, fraction, pident, length, a, b = best
                details = (pident, length, overlap, fraction, a, b)
            writer.writerow((cid, chrom, start, end, is_query, end - start, status, *details))

        matched = counts["matched"]
        rows.append((cid, len(mem), matched, matched / len(mem),
                     next(c for c, _, _, q in mem if q),
                     counts["fail_min_overlap"], counts["no_blast_hit"]))

with open(csv_path, "w", newline="") as ch:
    writer = csv.writer(ch)
    writer.writerow(["cluster_id", "n_members", "n_matched", "match_ratio", "query_chrom",
                     "n_fail_min_overlap", "n_no_blast_hit"])
    writer.writerows((cid, n, m, f"{ratio:.6f}", chrom, fraction, no_hit)
                     for cid, n, m, ratio, chrom, fraction, no_hit in rows)

ratios = np.array([r[3] for r in rows], dtype=float)
if ratios.size:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hist(ratios, bins=20, range=(0, 1), color="#4C72B0", edgecolor="white")
    axes[0].axvline(ratios.mean(), color="crimson", ls="--",
                    label=f"mean = {ratios.mean():.3f}")
    axes[0].set_xlabel("per-cluster match ratio")
    axes[0].set_ylabel("number of clusters")
    axes[0].set_title(f"SegTrace cluster BLAST recovery (n={ratios.size})")
    axes[0].legend()

    xs = np.sort(ratios)
    ys = np.arange(1, xs.size + 1) / xs.size
    axes[1].plot(xs, ys, color="#55A868")
    axes[1].set_xlabel("match ratio")
    axes[1].set_ylabel("cumulative fraction of clusters")
    axes[1].set_title("CDF")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)

    print(f"       clusters={ratios.size}  mean_ratio={ratios.mean():.4f}  "
          f"median={np.median(ratios):.4f}  fully_recovered="
          f"{int((ratios >= 1.0).sum())} ({(ratios >= 1.0).mean() * 100:.1f}%)")
    print("       unmatched members: " + "  ".join(
        f"{k}={unmatched[k]}" for k in ("fail_min_overlap", "no_blast_hit")))
else:
    print("       no clusters to score")
print(f"       CSV : {csv_path}")
print(f"       MEMBER RESULTS: {member_results_path}")
print(f"       PLOT: {plot_path}")
PY

echo "Done."
