#!/usr/bin/env bash
set -e

make clean && make
echo "=== Running Segtrace with Memory & Time Profiling ==="
/usr/bin/time -l ./segtrace -o t2t-nip data/t2t_nip.fasta
echo "=== Evaluation vs SEDEF Ground Truth ==="
python cmp_human.py --segtrace t2t-nip.dup.bed --sedef t2t_nip_sedef.sorted.bed