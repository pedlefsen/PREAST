#!/bin/bash
. <(./setup.sh)
export R_LIBS_USER=/home/tholzman/R/Library
export mainDir=$1
export patient=$2
export outputDir=${mainDir}/test-founder-inference-bakeoff_${patient}
mkdir $outputDir
rm -rf ${outputDir}/*
export listFile=${mainDir}/${patient}.list
export homeDir=`pwd`
while read fname
do
   echo $fname
   bn=`basename $fname .fasta`
   export outputFile=${outputDir}/${bn}.out
   export errFile=${outputDir}/${bn}.err
   echo $fname > $outputFile
   pushd .
   cd $outputDir
   /usr/bin/time -a -o $errFile $homeDir/bin/infer.sh -v -k  $fname >>$outputFile 2>>$errFile
   mv infer_multiple.fa ${outputFile}multiple.fa
   mv infer_single.fa   ${outputFile}single.fa
   mv infer_toi.csv     ${outputFile}toi.csv
   popd
done < $listFile


