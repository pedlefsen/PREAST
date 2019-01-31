#!/bin/bash
for patient in  `ls -c1 ${1}/*.fasta | egrep --only "_[0-9]{5,}_" | tr -d "_" | sort -u`
do 
   sbatch -JPREAST_${patient} -N1 -t 06:00 --mail-user=tholzman --mail-type=FAIL ./runList.bash ${1} $patient 
done
