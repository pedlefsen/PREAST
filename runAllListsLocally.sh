#!/bin/bash
for patient in  `ls -c1 ${1}/*.list  | egrep --only "[0-9]+\.list" | egrep --only "[0-9]+"  | sort -u`
do 
   ./runList.bash ${1} $patient 
done
