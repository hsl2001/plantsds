#!/usr/bin/env python3
"""Analysis of SegTrace plant shared genomic segments:
1) Genome & SD proportions: Total genome, 1-family, >=2-family shared segments (bp & %).
2) Comprehensive hit classification of >=2-family shared segments:
   - Independent hit metrics (Protein-coding 648, cp 53, mt 31, Unassigned 220).
   - Mutually exclusive Venn partition explaining organellar-nuclear coding overlaps.
3) 100-bin Chi-square hotspot test on cp (53 hits) and mt (31 hits) genomes.
4) Extraction of cross-family shared clusters containing Arabidopsis thaliana.
"""
import argparse
import gzip
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chisquare

ORGANELLES = {
    "cp": ("NC_000932.1", "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_000932.1&db=nuccore&report=fasta&retmode=text"),
    "mt": ("NC_037304.1", "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_037304.1&db=nuccore&report=fasta&retmode=text"),
}

CORE_PLANT_URLS = {
    "ath_rna": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_rna.fna.gz",
    "ath_pep": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz",
    "osa_rna": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_rna.fna.gz",
    "osa_pep": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_protein.faa.gz",
}

REF_FAMILIES = {
    "GCF_000471905.2": "Amborellaceae",
    "GCF_000001735.4": "Brassicaceae",
    "GCF_000005505.3": "Poaceae",
    "GCF_000309985.2": "Brassicaceae",
    "GCA_037997075.1": "Poaceae",
    "GCF_000002775.5": "Salicaceae",
    "GCF_036512215.1": "Solanaceae",
    "GCF_030704535.1": "Vitaceae",
}

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", default="ANGIOSPERM_FAMILY.dup.bed", help="SegTrace dup.bed")
    p.add_argument("--meta", default="angiosperm_family_complete.tsv", help="Complete summary TSV")
    p.add_argument("--output-dir", default="results/numt_nupt_cross_family", help="Output directory")
    p.add_argument("--threads", type=int, default=8, help="BLAST threads")
    p.add_argument("--evalue", type=float, default=1e-5, help="BLAST evalue threshold")
    p.add_argument("--nbins", type=int, default=100, help="Number of bins for hotspot test")
    return p.parse_args()


def download_file(url: str, dest_path: Path):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(dest_path, "wb") as out_f:
        shutil.copyfileobj(resp, out_f)


def prepare_databases(ref_dir: Path):
    ref_dir.mkdir(parents=True, exist_ok=True)
    dbs = {}

    for name, (acc, url) in ORGANELLES.items():
        fa = ref_dir / f"{acc}.fa"
        if not fa.exists():
            download_file(url, fa)
        db = ref_dir / acc
        if not list(ref_dir.glob(f"{acc}.n*")):
            subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(db)], check=True)
        dbs[name] = db

    core_rna = ref_dir / "plant_core_rna.fa"
    if not core_rna.exists():
        for tag in ["ath_rna", "osa_rna"]:
            gz = ref_dir / f"{tag}.fna.gz"
            fna = ref_dir / f"{tag}.fna"
            if not fna.exists():
                download_file(CORE_PLANT_URLS[tag], gz)
                with gzip.open(gz, "rb") as f_in, open(fna, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                gz.unlink(missing_ok=True)
        with open(core_rna, "w") as out:
            for tag in ["ath_rna", "osa_rna"]:
                with open(ref_dir / f"{tag}.fna") as f_in:
                    shutil.copyfileobj(f_in, out)
    if not list(ref_dir.glob("plant_core_rna.n*")):
        subprocess.run(["makeblastdb", "-in", str(core_rna), "-dbtype", "nucl", "-out", str(ref_dir / "plant_core_rna")], check=True)
    dbs["rna"] = ref_dir / "plant_core_rna"

    core_pep = ref_dir / "plant_core_pep.faa"
    if not core_pep.exists():
        for tag in ["ath_pep", "osa_pep"]:
            gz = ref_dir / f"{tag}.faa.gz"
            faa = ref_dir / f"{tag}.faa"
            if not faa.exists():
                download_file(CORE_PLANT_URLS[tag], gz)
                with gzip.open(gz, "rb") as f_in, open(faa, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                gz.unlink(missing_ok=True)
        with open(core_pep, "w") as out:
            for tag in ["ath_pep", "osa_pep"]:
                with open(ref_dir / f"{tag}.faa") as f_in:
                    shutil.copyfileobj(f_in, out)
    if not list(ref_dir.glob("plant_core_pep.p*")):
        subprocess.run(["makeblastdb", "-in", str(core_pep), "-dbtype", "prot", "-out", str(ref_dir / "plant_core_pep")], check=True)
    dbs["pep"] = ref_dir / "plant_core_pep"

    return dbs


def run_blast(program: str, query_fa: Path, db: Path, out_tsv: Path, threads: int, evalue: float):
    if not out_tsv.exists():
        cmd = [program, "-query", str(query_fa), "-db", str(db), "-out", str(out_tsv),
               "-outfmt", "6 " + " ".join(BLAST_COLS), "-evalue", str(evalue),
               "-max_target_seqs", "3", "-num_threads", str(threads)]
        if program == "blastn":
            cmd += ["-task", "blastn", "-dust", "no"]
        subprocess.run(cmd, check=True)
    if out_tsv.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLS)
    hits = pd.read_csv(out_tsv, sep="\t", names=BLAST_COLS)
    best = hits.sort_values(["qseqid", "bitscore"], ascending=[True, False]).groupby("qseqid", as_index=False).first()
    return best


def calc_hotspot(best_df: pd.DataFrame, sseqid: str, nbins: int):
    sub = best_df[best_df["sseqid"] == sseqid]
    if sub.empty:
        return None
    slen = int(sub["slen"].iloc[0])
    bins = np.linspace(0, slen, nbins + 1)
    mid = (sub["sstart"] + sub["send"]) / 2.0
    counts, _ = np.histogram(mid, bins=bins)

    exp = np.full(nbins, len(sub) / nbins)
    chi2 = float(((counts - exp) ** 2 / exp).sum())
    p_chi2 = float(chisquare(counts, exp).pvalue)

    bin_records = [
        {"bin": i + 1, "start": int(bins[i]), "end": int(bins[i + 1]), "count": int(counts[i])}
        for i in range(nbins)
    ]
    top_indices = np.argsort(counts)[::-1][:5]
    top_bins = [bin_records[i] for i in top_indices]

    return {"n_hits": len(sub), "chi2": chi2, "p_chi2": p_chi2, "top_bins": top_bins, "bins": bin_records}


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_dir = out_dir / "ref"

    # 1. Load Metadata & BED
    meta = pd.read_csv(args.meta, sep="\t")
    for acc, fam in REF_FAMILIES.items():
        meta.loc[meta["accession"] == acc, "family"] = fam

    bed = pd.read_csv(args.bed, sep="\t")
    bed["accession"] = bed["#chrom"].str.extract(r"^(GC[AF]_\d+\.\d+)")
    bed["seq_id"] = bed["#chrom"].str.split("-").str[-1]
    bed["length"] = bed["end"] - bed["start"]
    bed = bed.merge(meta[["accession", "family", "organism", "size_bp"]], on="accession", how="left")

    total_genome_bp = int(meta["size_bp"].sum())
    total_sd_bp = int(bed["length"].sum())

    fam_per_clus = bed.groupby("cluster_id")["family"].nunique()
    bed["n_fam"] = bed["cluster_id"].map(fam_per_clus)

    intra_mask = bed["n_fam"] == 1
    cross_mask = bed["n_fam"] >= 2
    intra_bp = int(bed.loc[intra_mask, "length"].sum())
    cross_bp = int(bed.loc[cross_mask, "length"].sum())

    print("=" * 80)
    print(" 1. GENOME & SD BP PROPORTIONS (TASK 1)")
    print("=" * 80)
    print(f"Total Genome Size     : {total_genome_bp:15,d} bp ({total_genome_bp/1e9:8.2f} Gb) | 100.0000%")
    print(f"Total SD              : {total_sd_bp:15,d} bp ({total_sd_bp/1e6:8.2f} Mb) | {total_sd_bp/total_genome_bp*100:8.4f}% of genome")
    print(f"  - 1-Family (Intra)  : {intra_bp:15,d} bp ({intra_bp/1e6:8.2f} Mb) | {intra_bp/total_genome_bp*100:8.4f}% of genome | {intra_bp/total_sd_bp*100:6.2f}% of SD")
    print(f"  - >=2-Family (Cross): {cross_bp:15,d} bp ({cross_bp/1e6:8.2f} Mb) | {cross_bp/total_genome_bp*100:8.4f}% of genome | {cross_bp/total_sd_bp*100:6.2f}% of SD")

    # 2. Filter Cross-Family (>=2 families) Segments & Extract Queries
    cross_bed = bed[cross_mask].copy()
    cross_bed["n_species"] = cross_bed.groupby("cluster_id")["accession"].transform("nunique")
    cross_bed["qid"] = ("cl" + cross_bed["cluster_id"].astype(str)
                        + "|nsp" + cross_bed["n_species"].astype(str)
                        + "|" + cross_bed["accession"].astype(str)
                        + "|" + cross_bed["seq_id"].astype(str)
                        + "|" + cross_bed["start"].astype(str) + "-" + cross_bed["end"].astype(str))

    query_fa = out_dir / "queries.fa"
    if not query_fa.exists() or query_fa.stat().st_size == 0:
        print(f"\n[info] Downloading {len(cross_bed)} query sequences via NCBI...")
        with open(query_fa, "w") as out:
            for r in cross_bed.itertuples():
                url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&val={r.seq_id}&from={r.start+1}&to={r.end}&fmt_mask=0&report=fasta&retmode=text"
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=15) as resp:
                    txt = resp.read().decode("utf-8").strip()
                seq = "".join(txt.splitlines()[1:])
                out.write(f">{r.qid}\n{seq}\n")

    # 3. Multi-Tier BLAST Searching
    dbs = prepare_databases(ref_dir)
    best_cp = run_blast("blastn", query_fa, dbs["cp"], out_dir / "hits_cp.tsv", args.threads, args.evalue)
    best_mt = run_blast("blastn", query_fa, dbs["mt"], out_dir / "hits_mt.tsv", args.threads, args.evalue)
    best_rna = run_blast("blastn", query_fa, dbs["rna"], out_dir / "hits_plant_core_rna.tsv", args.threads, args.evalue)
    best_pep = run_blast("blastx", query_fa, dbs["pep"], out_dir / "hits_plant_core_pep.tsv", args.threads, args.evalue)

    # 4. Independent Hit Metrics & Complete Disjoint Partition
    cp_set = set(best_cp["qseqid"])
    mt_set = set(best_mt["qseqid"])
    prot_set = set(best_rna["qseqid"]).union(set(best_pep["qseqid"]))

    cross_bed["has_cp"] = cross_bed["qid"].isin(cp_set)
    cross_bed["has_mt"] = cross_bed["qid"].isin(mt_set)
    cross_bed["has_prot"] = cross_bed["qid"].isin(prot_set)
    cross_bed["is_unassigned"] = (~cross_bed["has_cp"]) & (~cross_bed["has_mt"]) & (~cross_bed["has_prot"])

    print("\n" + "=" * 80)
    print(" 2. >=2-FAMILY SHARED SEGMENTS: INDEPENDENT HIT BREAKDOWN (TASK 2)")
    print("=" * 80)
    for name, col in [("Protein-coding gene hit", "has_prot"),
                      ("Chloroplast (cp) hit", "has_cp"),
                      ("Mitochondria (mt) hit", "has_mt"),
                      ("Unassigned (No hit to any)", "is_unassigned")]:
        sub = cross_bed[cross_bed[col]]
        n = len(sub)
        bp = int(sub["length"].sum())
        print(f"  - {name:30s}: {bp:10,d} bp ({bp/1e3:6.1f} kb) | {n:4d} segs ({n/len(cross_bed)*100:5.2f}%) | {bp/cross_bp*100:6.2f}% of cross-SD | {bp/total_genome_bp*100:8.6f}% of genome")

    def get_venn_cat(r):
        c, m, p = r["has_cp"], r["has_mt"], r["has_prot"]
        if p and not c and not m:
            return "Nuclear Protein-coding only"
        elif c and not m and not p:
            return "Chloroplast (cp) only"
        elif m and not c and not p:
            return "Mitochondria (mt) only"
        elif c and p and not m:
            return "cp + Protein-coding shared"
        elif m and p and not c:
            return "mt + Protein-coding shared"
        elif c and m and p:
            return "cp + mt + Protein-coding shared"
        else:
            return "Unassigned"

    cross_bed["venn_category"] = cross_bed.apply(get_venn_cat, axis=1)
    cross_bed.to_csv(out_dir / "cross_family_classified.tsv", sep="\t", index=False)

    print("\n" + "-" * 80)
    print(" [Detailed Mutually Exclusive Partition]")
    print("-" * 80)
    for cat, sub in cross_bed.groupby("venn_category"):
        n = len(sub)
        bp = int(sub["length"].sum())
        print(f"  * {cat:32s}: {bp:10,d} bp ({bp/1e3:6.1f} kb) | {n:4d} segs ({n/len(cross_bed)*100:5.2f}%) | {bp/cross_bp*100:5.2f}% of cross-SD")

    # 3. 100-bin Chi-square Hotspot Test
    hs_cp = calc_hotspot(best_cp, ORGANELLES["cp"][0], args.nbins)
    hs_mt = calc_hotspot(best_mt, ORGANELLES["mt"][0], args.nbins)

    print("\n" + "=" * 80)
    print(" 3. ORGANELLE HOTSPOT CHI-SQUARE TEST (100 BINS) (TASK 3)")
    print("=" * 80)
    if hs_cp:
        top_str = ", ".join(f"bin{b['bin']}({b['start']}-{b['end']}bp)={b['count']}" for b in hs_cp["top_bins"][:3])
        print(f"  - Chloroplast   (cp): hits={hs_cp['n_hits']:3d} | chi2={hs_cp['chi2']:8.2f} | p={hs_cp['p_chi2']:.3e}")
        print(f"    * Top hotspots: {top_str}")
        pd.DataFrame(hs_cp["bins"]).to_csv(out_dir / f"bins_cp_{args.nbins}.tsv", sep="\t", index=False)
    if hs_mt:
        top_str = ", ".join(f"bin{b['bin']}({b['start']}-{b['end']}bp)={b['count']}" for b in hs_mt["top_bins"][:3])
        print(f"  - Mitochondrion (mt): hits={hs_mt['n_hits']:3d} | chi2={hs_mt['chi2']:8.2f} | p={hs_mt['p_chi2']:.3e}")
        print(f"    * Top hotspots: {top_str}")
        pd.DataFrame(hs_mt["bins"]).to_csv(out_dir / f"bins_mt_{args.nbins}.tsv", sep="\t", index=False)

    # 4. Extract Arabidopsis-Containing Cross-Family Clusters
    ath_clusters = cross_bed[cross_bed["organism"].str.contains("Arabidopsis thaliana", case=False, na=False)]["cluster_id"].unique()
    ath_sub = cross_bed[cross_bed["cluster_id"].isin(ath_clusters)]

    ath_fa_path = out_dir / "arabidopsis_cross_family_clusters.fa"
    with open(ath_fa_path, "w") as out:
        for r in ath_sub.itertuples():
            url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&val={r.seq_id}&from={r.start+1}&to={r.end}&fmt_mask=0&report=fasta&retmode=text"
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                txt = resp.read().decode("utf-8").strip()
            seq = "".join(txt.splitlines()[1:])
            header = f"{r.organism.split('(')[0].strip()}_{r.family}_cl{r.cluster_id}_{r.seq_id}_{r.start}_{r.end}"
            out.write(f">{header}\n{seq}\n")

    print("\n" + "=" * 80)
    print(" 4. ARABIDOPSIS-CONTAINING CROSS-FAMILY CLUSTERS (TASK 4)")
    print("=" * 80)
    print(f"  - Clusters found: {len(ath_clusters)} (Cluster IDs: {list(ath_clusters)})")
    for r in ath_sub.itertuples():
        print(f"    * [{r.family}] {r.organism}: {r.seq_id}:{r.start}-{r.end} ({r.length} bp)")
    print(f"  - Saved all sequence FASTA to: {ath_fa_path}")
    print("=" * 80)


if __name__ == "__main__":
    main()
