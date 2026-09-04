#!/usr/bin/env bash
# SegTrace cluster validation with BLAST.
#
# For every SegTrace cluster it takes the longest member as the representative,
# BLASTs it against the whole genome set, and measures how many of the cluster
# members are recovered by a BLAST hit (per-cluster match ratio). The ratios are
# written to a CSV and plotted.
set -euo pipefail

if [ -z "${FASTAS+x}" ]; then
  FASTAS=(tair12.fasta t2t_nip.fasta)
fi
OUTDIR="${OUTDIR:-segtrace_validation}"
PREFIX="${PREFIX:-$OUTDIR/segtrace}"
THREADS="${THREADS:-8}"
MIN_OVERLAP="${MIN_OVERLAP:-0.5}"   # member covered fraction to count as a match
MIN_IDENT="${MIN_IDENT:-80}"        # BLAST percent-identity cutoff (segtrace ~0.8)
SEGTRACE_EXTRA="${SEGTRACE_EXTRA:--c 1}"

BED="$PREFIX.seg.bed"
COMBINED="$OUTDIR/combined.fa"
QUERY="$OUTDIR/cluster_reps.fa"
MEMBERS="$OUTDIR/cluster_members.tsv"
DB="$OUTDIR/blastdb/combined"
BLAST="$OUTDIR/blast_hits.tsv"
CSV="$OUTDIR/cluster_match_ratio.csv"
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
python3 - "$COMBINED" "$QUERY" "$MEMBERS" "$BED" "${FASTAS[@]}" <<'PY'
import os
import sys

combined_path, query_path, members_path, bed_path = sys.argv[1:5]
fastas = sys.argv[5:]
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
with open(bed_path) as bh:
    for line in bh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        clusters.setdefault(int(f[3]), []).append((f[0], int(f[1]), int(f[2])))

with open(combined_path, "rb") as cf, open(query_path, "wb") as qf, \
        open(members_path, "w") as mf:
    mf.write("cluster_id\tchrom\tstart\tend\tis_query\n")
    for cid in sorted(clusters):
        members = clusters[cid]
        qi = max(range(len(members)), key=lambda i: members[i][2] - members[i][1])
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

print(f"       representatives: {len(clusters)}")
PY

# --------------------------------------------------------- 3. BLAST DB + search
echo "[3/5] Building BLAST database..."
makeblastdb -dbtype nucl -in "$COMBINED" -out "$DB" >/dev/null

echo "[4/5] BLASTing cluster representatives against the genomes..."
blastn -task megablast -query "$QUERY" -db "$DB" -num_threads "$THREADS" \
  -perc_identity "$MIN_IDENT" -evalue 1e-5 -max_target_seqs 100000 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
  -out "$BLAST"

# ------------------------------------------------ 5. match ratio + CSV + graph
echo "[5/5] Scoring cluster recovery and plotting..."
python3 - "$MEMBERS" "$BLAST" "$CSV" "$PLOT" "$MIN_OVERLAP" "$MIN_IDENT" <<'PY'
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

members_path, blast_path, csv_path, plot_path, min_overlap_s, min_ident_s = sys.argv[1:7]
min_overlap = float(min_overlap_s)
min_ident = float(min_ident_s)

members = defaultdict(list)  # cid -> [(chrom, start, end, is_query)]
with open(members_path) as mh:
    next(mh)
    for line in mh:
        cid, chrom, s, e, isq = line.rstrip("\n").split("\t")
        members[int(cid)].append((chrom, int(s), int(e), int(isq)))

hits = defaultdict(list)  # (cid, chrom) -> [(start, end)] on subject
with open(blast_path) as bh:
    for line in bh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or float(f[2]) < min_ident:
            continue
        sstart, send = int(f[6]), int(f[7])
        hits[(int(f[0][1:]), f[1])].append((min(sstart, send) - 1, max(sstart, send)))


def covered(chrom, start, end, intervals):
    mlen = end - start
    if mlen <= 0:
        return False
    for a, b in intervals:
        overlap = min(end, b) - max(start, a)
        if overlap > 0 and overlap / mlen >= min_overlap:
            return True
    return False


rows = []
for cid in sorted(members):
    mem = members[cid]
    matched = sum(covered(c, s, e, hits.get((cid, c), [])) for c, s, e, _ in mem)
    qchrom = next((c for c, _, _, q in mem if q), mem[0][0])
    rows.append((cid, len(mem), matched, matched / len(mem), qchrom))

with open(csv_path, "w", newline="") as ch:
    writer = csv.writer(ch)
    writer.writerow(["cluster_id", "n_members", "n_matched", "match_ratio", "query_chrom"])
    for cid, n, m, ratio, qchrom in rows:
        writer.writerow([cid, n, m, f"{ratio:.6f}", qchrom])

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
else:
    print("       no clusters to score")
print(f"       CSV : {csv_path}")
print(f"       PLOT: {plot_path}")
PY

echo "Done."
