#!/bin/bash
for patient in  `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do 
   sbatch -JPREAST_${patient} -N1 -t 06:00 --mail-user=tholzman --mail-type=FAIL ./runList.bash ${1} $patient 
done
