echo "cd `pwd`; ./time ./run_sd.sh -t 128 -s t2t-nip" | qsub -N segtrace-species -l nodes=node02:ppn=128 -v WORKDIR=`pwd` -j oe -o ~/log/segtrace-species.log
