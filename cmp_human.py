#!/usr/bin/env python3
"""
cmp_human.py - Comparison CLI for real human / pangenome datasets (Segtrace vs SEDEF / BISER).

Usage:
  python3 cmp_human.py --segtrace results/t2t-chm13_sd.dup.bed --sedef sedef-human-t2tchm13.bed
"""

import sys
import os
import argparse
from cmp_core import run_sd_core_comparison

def main():
    parser = argparse.ArgumentParser(description="Compare Segtrace and SEDEF BED files on Human/Real genome datasets.")
    parser.add_argument("--segtrace", default="results/t2t-chm13_sd.dup.bed", help="Path to Segtrace dup.bed file")
    parser.add_argument("--sedef", default="sedef-human-t2tchm13.bed", help="Path to SEDEF bed file")
    parser.add_argument("--out-renamed", default=None, help="Optional output path to save renamed Segtrace BED file")
    parser.add_argument("--work-dir", default="_cmp_tmp", help="Temporary working directory")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary intermediate files")
    args = parser.parse_args()

    if not os.path.exists(args.segtrace):
        print(f"[ERROR] Segtrace file '{args.segtrace}' not found.")
        sys.exit(1)
    if not os.path.exists(args.sedef):
        print(f"[ERROR] SEDEF file '{args.sedef}' not found.")
        sys.exit(1)

    run_sd_core_comparison(args.segtrace, args.sedef, out_renamed=args.out_renamed, work_dir=args.work_dir, keep_temp=args.keep_temp)

if __name__ == "__main__":
    main()
