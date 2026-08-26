#!/usr/bin/env python3
"""NUMT/NUPT analysis of segtrace hits via bedtools getfasta + blastn.

1) Extract all segtrace hit sequences with bedtools getfasta (with automatic contig header resolution).
2) blastn each query set against reference cp and mt genomes with a single evalue threshold.
3) Calculate match proportions of segtrace hits (none / cp / mt / both).
4) 100-bin profile of hit counts across cp and mt reference genomes + hotspot test.

Requires: bedtools, blastn, makeblastdb, samtools in PATH.
"""
import argparse
import functools
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chisquare

REFS = {
    "cp": ("NC_000932.1", "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_000932.1&db=nuccore&report=fasta&retmode=text"),
    "mt": ("NC_037304.1", "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_037304.1&db=nuccore&report=fasta&retmode=text"),
}
COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", default="ANGIOSPERM_FAMILY.dup.bed")
    p.add_argument("--summary-tsv", default="angiosperm_family_min.tsv", help="TSV with accession,fasta columns")
    p.add_argument("--output-dir", default="results/numt_nupt")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--evalue", type=float, default=1e-5, help="single blastn evalue threshold")
    p.add_argument("--nbins", type=int, default=100, help="number of bins for organelle reference genomes (default: 100)")
    p.add_argument("--mc-permutations", type=int, default=100000)
    p.add_argument("--max-queries", type=int, default=0, help="subsample BED rows (0 = all)")
    return p.parse_args()


def run(cmd):
    subprocess.check_call(cmd)


def run_capture(cmd):
    return subprocess.check_output(cmd, text=True)


def parse_fasta(text):
    header, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def ensure_refs(out_dir):
    ref_dir = out_dir / "ref"
    ref_dir.mkdir(parents=True, exist_ok=True)
    dbs = {}
    for name, (acc, url) in REFS.items():
        fa = ref_dir / f"{acc}.fa"
        if not fa.exists():
            with urllib.request.urlopen(url, timeout=60) as r:
                txt = r.read().decode()
            assert txt.startswith(">"), f"download failed for {acc}"
            fa.write_text(txt)
        db = ref_dir / acc
        if not list(ref_dir.glob(f"{acc}.n*")):
            run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(db)])
        dbs[name] = db
    return dbs


def resolve_fasta_path(raw_path: str, summary_tsv: Path) -> Path:
    p = Path(raw_path)
    if p.is_absolute() and p.exists():
        return p
    candidates = [
        Path.cwd() / p,
        p,
        summary_tsv.resolve().parent / p,
        summary_tsv.resolve().parent.parent / p,
        Path.cwd() / "selected" / p,
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    return (Path.cwd() / p).resolve()


@functools.lru_cache(maxsize=None)
def get_fai_contigs(fasta_path: str):
    fpath = Path(fasta_path)
    fai = fpath.with_name(fpath.name + ".fai")
    if not fai.exists():
        run(["samtools", "faidx", str(fpath)])
    with open(fai) as fh:
        return frozenset(line.split("\t", 1)[0] for line in fh if line.strip())


def resolve_chrom_name(raw_chrom: str, seq_id: str, fai_names):
    if raw_chrom in fai_names:
        return raw_chrom
    if seq_id in fai_names:
        return seq_id
    suffix = raw_chrom.split("-")[-1]
    if suffix in fai_names:
        return suffix
    for name in fai_names:
        if name.endswith(suffix) or suffix.endswith(name):
            return name
    return None


def load_bed(bed_path, max_queries):
    bed = pd.read_csv(bed_path, sep="\t")
    bed["accession"] = bed["#chrom"].str.extract(r"^(GC[AF]_\d+\.\d+)")
    bed["seq_id"] = bed["#chrom"].str.split("-").str[-1]
    assert not bed[["accession", "seq_id"]].isna().any().any(), "cannot parse #chrom as ACC-seqid"
    bed["n_species"] = bed.groupby("cluster_id")["accession"].transform("nunique")
    if max_queries > 0 and len(bed) > max_queries:
        bed = bed.sample(n=max_queries, random_state=42).sort_index()
    bed["qid"] = ("cl" + bed["cluster_id"].astype(str)
                  + "|nsp" + bed["n_species"].astype(str)
                  + "|" + bed["accession"].astype(str)
                  + "|" + bed["seq_id"].astype(str)
                  + "|" + bed["start"].astype(str) + "-" + bed["end"].astype(str))
    return bed


def extract_fasta(bed, acc_to_fasta, summary_tsv, out_fa):
    written = 0
    missing_fastas = set()
    no_contig_fastas = set()
    with open(out_fa, "w") as out:
        for acc, sub in bed.groupby("accession"):
            raw_f = acc_to_fasta.get(acc)
            if not raw_f:
                print(f"[warn] no fasta mapping for accession {acc}", file=sys.stderr)
                continue
            fasta = resolve_fasta_path(str(raw_f), summary_tsv)
            if not fasta.exists():
                missing_fastas.add(f"{acc} -> {fasta}")
                continue

            fai_names = get_fai_contigs(str(fasta))
            tmp = out_fa.with_suffix(f".{acc}.bed")
            valid_rows = 0
            with open(tmp, "w") as fh:
                for r in sub.itertuples():
                    cname = resolve_chrom_name(r._1, r.seq_id, fai_names)  # r._1 is '#chrom'
                    if cname is not None:
                        fh.write(f"{cname}\t{r.start}\t{r.end}\t{r.qid}\n")
                        valid_rows += 1

            if valid_rows == 0:
                no_contig_fastas.add(f"{acc} ({fasta.name})")
                tmp.unlink(missing_ok=True)
                continue

            txt = run_capture(["bedtools", "getfasta", "-fi", str(fasta), "-bed", str(tmp), "-name", "-fo", "-"])
            tmp.unlink(missing_ok=True)
            for header, seq in parse_fasta(txt):
                out.write(f">{header.split('::')[0]}\n{seq}\n")
                written += 1

    if missing_fastas:
        print(f"[warn] {len(missing_fastas)} FASTA files were not found on disk:", file=sys.stderr)
        for m in sorted(missing_fastas)[:5]:
            print(f"       {m}", file=sys.stderr)
        if len(missing_fastas) > 5:
            print(f"       ... and {len(missing_fastas)-5} more", file=sys.stderr)
    if no_contig_fastas:
        print(f"[warn] {len(no_contig_fastas)} FASTA files had 0 matching contigs with BED:", file=sys.stderr)
        for m in sorted(no_contig_fastas)[:5]:
            print(f"       {m}", file=sys.stderr)
    return written


def run_blastn(query_fa, db, out_tsv, threads, evalue):
    if not out_tsv.exists():
        run(["blastn", "-task", "blastn", "-query", str(query_fa), "-db", str(db), "-out", str(out_tsv),
             "-outfmt", "6 " + " ".join(COLS), "-evalue", str(evalue), "-dust", "no",
             "-max_target_seqs", "5", "-num_threads", str(threads)])
    if out_tsv.stat().st_size == 0:
        return pd.DataFrame(columns=COLS)
    hits = pd.read_csv(out_tsv, sep="\t", names=COLS)
    best = hits.sort_values(["qseqid", "bitscore"], ascending=[True, False]).groupby("qseqid", as_index=False).first()
    return best


def classify_queries(bed_df, best_cp, best_mt):
    """Classify each query into none, cp_only, mt_only, both (and primary assignment by bitscore)."""
    df = bed_df[["qid", "n_species", "cluster_id"]].drop_duplicates("qid").copy()
    cp_map = best_cp.set_index("qseqid")["bitscore"].to_dict() if not best_cp.empty else {}
    mt_map = best_mt.set_index("qseqid")["bitscore"].to_dict() if not best_mt.empty else {}

    df["cp_bitscore"] = df["qid"].map(cp_map).fillna(0.0)
    df["mt_bitscore"] = df["qid"].map(mt_map).fillna(0.0)
    df["has_cp"] = df["cp_bitscore"] > 0
    df["has_mt"] = df["mt_bitscore"] > 0

    def category(row):
        if row["has_cp"] and row["has_mt"]:
            return "both"
        if row["has_cp"]:
            return "cp"
        if row["has_mt"]:
            return "mt"
        return "none"

    def primary(row):
        if not row["has_cp"] and not row["has_mt"]:
            return "none"
        return "cp" if row["cp_bitscore"] >= row["mt_bitscore"] else "mt"

    df["class_detailed"] = df.apply(category, axis=1)
    df["class_primary"] = df.apply(primary, axis=1)
    return df


def calc_proportions(classified_df):
    results = {}
    groups = {
        "all": classified_df,
        "shared": classified_df[classified_df["n_species"] >= 2],
        "singleton": classified_df[classified_df["n_species"] == 1],
    }
    for gname, sub in groups.items():
        tot = len(sub)
        counts_det = sub["class_detailed"].value_counts().to_dict()
        counts_pri = sub["class_primary"].value_counts().to_dict()
        results[gname] = {
            "total_queries": tot,
            "none": {"count": int(counts_pri.get("none", 0)), "rate": float(counts_pri.get("none", 0) / tot) if tot else 0.0},
            "cp": {"count": int(counts_pri.get("cp", 0)), "rate": float(counts_pri.get("cp", 0) / tot) if tot else 0.0},
            "mt": {"count": int(counts_pri.get("mt", 0)), "rate": float(counts_pri.get("mt", 0) / tot) if tot else 0.0},
            "detailed": {
                "none": int(counts_det.get("none", 0)),
                "cp_only": int(counts_det.get("cp", 0)),
                "mt_only": int(counts_det.get("mt", 0)),
                "both": int(counts_det.get("both", 0)),
            },
        }
    return results


def binned_hotspot(best, sseqid, nbins, nperm):
    """Compute hit count in each of the nbins along the reference genome."""
    if best.empty:
        return None
    sub = best[best["sseqid"] == sseqid].copy()
    if len(sub) == 0:
        return None
    slen = int(sub["slen"].iloc[0])
    bins = np.linspace(0, slen, nbins + 1)
    mid = (sub["sstart"] + sub["send"]) / 2.0
    counts, _ = np.histogram(mid, bins=bins)

    all_bins = [
        {"bin": int(i + 1), "start": int(bins[i]), "end": int(bins[i + 1]), "count": int(counts[i])}
        for i in range(nbins)
    ]
    exp = np.full(nbins, len(sub) / nbins)
    chi2 = float(((counts - exp) ** 2 / exp).sum())
    p_chi2 = float(chisquare(counts, exp).pvalue)
    rng = np.random.default_rng(42)
    sims = rng.multinomial(len(sub), np.full(nbins, 1 / nbins), size=nperm)
    p_mc = float(((((sims - exp) ** 2 / exp).sum(axis=1) >= chi2).sum() + 1) / (nperm + 1))
    top = np.argsort(counts)[::-1][:5]
    top_bins = [all_bins[i] for i in top]

    return {
        "n_hits": int(len(sub)),
        "ref_len": slen,
        "nbins": nbins,
        "chi2": chi2,
        "p_chi2": p_chi2,
        "p_mc": p_mc,
        "top_bins": top_bins,
        "bins": all_bins,
    }


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dbs = ensure_refs(out_dir)

    bed = load_bed(Path(args.bed), args.max_queries)
    print(f"[info] bed rows={len(bed)} clusters={bed['cluster_id'].nunique()}")

    summary_tsv = Path(args.summary_tsv)
    meta = pd.read_csv(summary_tsv, sep="\t")
    acc_to_fasta = {str(r.accession): r.fasta for r in meta.itertuples() if pd.notna(r.fasta)}

    query_fa = out_dir / "queries.fa"
    if query_fa.exists() and query_fa.stat().st_size > 0:
        n_written = sum(1 for line in open(query_fa) if line.startswith(">"))
        print(f"[info] reused existing {n_written} queries -> {query_fa}")
    else:
        query_fa.unlink(missing_ok=True)
        n_written = extract_fasta(bed, acc_to_fasta, summary_tsv, query_fa)
        print(f"[info] extracted {n_written}/{len(bed)} sequences -> {query_fa}")

    if n_written == 0:
        query_fa.unlink(missing_ok=True)
        raise SystemExit(f"[error] no sequences extracted from {len(bed)} BED rows. Check FASTA file existence and contig names.")

    # Step 2: BLAST vs cp and mt with single evalue
    best_cp = run_blastn(query_fa, dbs["cp"], out_dir / "hits_cp.tsv", args.threads, args.evalue)
    best_mt = run_blastn(query_fa, dbs["mt"], out_dir / "hits_mt.tsv", args.threads, args.evalue)
    if not best_cp.empty:
        best_cp.to_csv(out_dir / "best_cp.tsv", sep="\t", index=False)
    if not best_mt.empty:
        best_mt.to_csv(out_dir / "best_mt.tsv", sep="\t", index=False)
    print(f"[info] cp hits: {len(best_cp)} queries (evalue<={args.evalue})")
    print(f"[info] mt hits: {len(best_mt)} queries (evalue<={args.evalue})")

    # Step 3: Match proportions (none, cp, mt, both)
    classified = classify_queries(bed, best_cp, best_mt)
    classified.to_csv(out_dir / "queries_classified.tsv", sep="\t", index=False)
    proportions = calc_proportions(classified)

    # Step 4: 100-bin profile for cp and mt
    hotspots = {
        "cp": binned_hotspot(best_cp, REFS["cp"][0], args.nbins, args.mc_permutations),
        "mt": binned_hotspot(best_mt, REFS["mt"][0], args.nbins, args.mc_permutations),
    }
    for name, h in hotspots.items():
        if h:
            pd.DataFrame(h["bins"]).to_csv(out_dir / f"bins_{name}_{args.nbins}.tsv", sep="\t", index=False)

    out = {
        "input": {"bed": str(args.bed), "n_queries": n_written},
        "evalue": args.evalue,
        "nbins": args.nbins,
        "proportions": proportions,
        "hotspots": hotspots,
    }
    with open(out_dir / "summary.json", "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"\n[done] output_dir={out_dir}")
    print("=" * 70)
    print("1. Match proportions (none / cp / mt):")
    for gname, r in proportions.items():
        print(f"  [{gname}] total={r['total_queries']} | "
              f"none={r['none']['count']} ({r['none']['rate']:.2%}) | "
              f"cp={r['cp']['count']} ({r['cp']['rate']:.2%}) | "
              f"mt={r['mt']['count']} ({r['mt']['rate']:.2%}) | "
              f"(both_overlap={r['detailed']['both']})")

    print("\n2. Organelle genome 100-bin hotspot summary:")
    for name, h in hotspots.items():
        if h:
            top_desc = ", ".join(f"bin{b['bin']}({b['start']}-{b['end']}bp)={b['count']}" for b in h["top_bins"][:3])
            print(f"  [{name}] n={h['n_hits']}, chi2={h['chi2']:.2f}, p_chi2={h['p_chi2']:.3e}, p_mc={h['p_mc']:.3e}")
            print(f"       Top bins: {top_desc}")
            print(f"       Saved all {args.nbins} bins to: {out_dir}/bins_{name}_{args.nbins}.tsv")
        else:
            print(f"  [{name}] no hits found")
    print("=" * 70)


if __name__ == "__main__":
    main()
