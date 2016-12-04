#!/bin/bash 
#PBS -N myjob
#PBS -o /home/kong/github/interface_HW/info/${a}_${b}_${k}_${t}_${x}_${d}
#PBS -e /home/kong/github/interface_HW/err/${a}_${b}_${k}_${t}_${x}_${d}
#PBS -l nodes=1:ppn=4,walltime=2:00:00
cd /home/kong/github/interface_HW
./main $a $b $k $t $x $d ./data/output_${a}_${b}_${k}_${t}_${x}_${d}
