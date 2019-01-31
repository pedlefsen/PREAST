#!/bin/bash
##
# First argument is the tabfile with the bounds and PIDs in it.
# Second argument is the directory with the lists of fasta files in it
#
# --TAH 3/16
##
export tabFile=$1
export tabBase=`basename $tabFile .tab`
export dataDir=$2
while read pinfo
do
   if [[ "$pinfo" == *"lower"* ]]; then
      continue
   fi
   if [[ "$pinfo" == *"NA"* ]]; then
      continue
   fi
   export pn=`echo "$pinfo" | cut  -f1`
   export lb=`echo "$pinfo" | cut  -f2`
   export ub=`echo "$pinfo" | cut  -f3`

   if [[ "$lb" == *"-"* ]]; then
      export lb=1
   fi   

   sbatch -JPREAST_${pn} -N1 -t 06:00 -A edlefsen_p --mail-user=tholzman --mail-type=FAIL ./runListBounds.bash $dataDir $tabBase $pn $lb $ub 

done < $tabFile
