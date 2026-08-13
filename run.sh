#!/usr/bin/env bash
# ==============================================================================
# run.sh - Segtrace qsub Job Submission Script for Plant Pangenomes (11 Projects)
# Distributed across node02 and node03
# ==============================================================================

WORKDIR="$(pwd)"
mkdir -p ~/log

# PBS Node & PPN configuration
NODE02="nodes=node02:ppn=128"
NODE03="nodes=node03:ppn=128"

# ------------------------------------------------------------------------------
# 1. Maize (Zea mays NAM/T2T Pangenome) -> node02
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o MAIZE ${WORKDIR}/pangenome_data/maize/assemblies/*.fna.gz" | qsub -N segtrace-maize -l ${NODE02} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-maize.log

# ------------------------------------------------------------------------------
# 2. Wheat (Triticum aestivum 10+ Pangenome) -> node03
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o WHEAT ${WORKDIR}/pangenome_data/wheat/assemblies/*.fa.gz" | qsub -N segtrace-wheat -l ${NODE03} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-wheat.log

# ------------------------------------------------------------------------------
# 3. Wild Grape (Vitis spp. Super-Pangenome) -> node02
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o WILD_GRAPE ${WORKDIR}/pangenome_data/wild_grape/assemblies/*.fa*" | qsub -N segtrace-grape -l ${NODE02} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-grape.log

# ------------------------------------------------------------------------------
# 4. Rapeseed (Brassica napus SV Pangenome) -> node03
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o RAPESEED ${WORKDIR}/pangenome_data/rapeseed/assemblies/*.fa.gz" | qsub -N segtrace-rapeseed -l ${NODE03} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-rapeseed.log

# ------------------------------------------------------------------------------
# 5. Potato (Solanum tuberosum Tetraploid Pangenome) -> node02
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o POTATO ${WORKDIR}/pangenome_data/potato/assemblies/*.fna.gz" | qsub -N segtrace-potato -l ${NODE02} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-potato.log

# ------------------------------------------------------------------------------
# 6. Eggplant (Solanum melongena Pangenome) -> node03
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o EGGPLANT ${WORKDIR}/pangenome_data/eggplant/assemblies/*.fna.gz" | qsub -N segtrace-eggplant -l ${NODE03} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-eggplant.log

# ------------------------------------------------------------------------------
# 7. Tea Plant (Camellia sinensis Haplotype Pangenome) -> node02
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o TEA ${WORKDIR}/pangenome_data/tea_plant/assemblies/*.fa*" | qsub -N segtrace-tea -l ${NODE02} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-tea.log

# ------------------------------------------------------------------------------
# 8. Watermelon (Citrullus lanatus Super-Pangenome) -> node03
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o WATERMELON ${WORKDIR}/pangenome_data/watermelon/assemblies/*.fa*" | qsub -N segtrace-watermelon -l ${NODE03} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-watermelon.log

# ------------------------------------------------------------------------------
# 9. Tomato (Solanum lycopersicum T2T Super-Pangenome) -> node02
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o TOMATO-149 ${WORKDIR}/pangenome_data/tomato/assemblies/*.fa.gz" | qsub -N Sl -l ${NODE02} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-Sl.log

# ------------------------------------------------------------------------------
# 10. Marchantia (Marchantia polymorpha Pangenome) -> node03
# ------------------------------------------------------------------------------
# echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o MARCHANTIA ${WORKDIR}/pangenome_data/marchantia/assemblies/*.fa*" | qsub -N segtrace-marchantia -l ${NODE03} -v WORKDIR=${WORKDIR} -j oe -o ~/log/segtrace-marchantia.log

# ------------------------------------------------------------------------------
# 11. Citrus (Citrus spp. Pangenome) -> node02
# ------------------------------------------------------------------------------
echo "cd ${WORKDIR}; ./time ./segtrace -p 128 -o CITRUS ${WORKDIR}/pangenome_data/citrus/assemblies/*.genome.fa ${WORKDIR}/pangenome_data/citrus/assemblies/*.chromosome.fa" | qsub -N Citrus -l ${NODE02} -v WORKDIR=${WORKDIR} -j oe -o ~/log/CITRUS-Vv.log