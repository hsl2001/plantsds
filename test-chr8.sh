#!/usr/bin/env bash
set -e

make clean && make
echo "=== Running Segtrace with Memory & Time Profiling ==="
/usr/bin/time -l ./segtrace -o data/t2t_chm13v2.0_chr8_sd data/t2t_chm13v2.0_chr8.fna.gz
echo "=== Evaluation vs SEDEF Ground Truth ==="
python3 cmp_human.py --segtrace data/t2t_chm13v2.0_chr8_sd.dup.bed --sedef data/t2t_chm13v2.0_chr8_SD.bed