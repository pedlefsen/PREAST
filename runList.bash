#!/bin/bash
. <(./setup.sh)
export R_LIBS_USER=/home/tholzman/R/Library
export mainDir=$1
export patient=$2
export outputDir=${mainDir}/founder-inference-bakeoff_${patient}
#rm -rf ${outputDir}
mkdir $outputDir
export listFile=${mainDir}/processed_${patient}.list
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
   $homeDir/bin/infer.sh $fname >>$outputFile 2>>$errFile
   mv infer_multiple.fa  ${outputFile}multiple.fa
   mv infer_single.fa    ${outputFile}single.fa
   mv infer_toi.csv      ${outputFile}toi.csv
   mv infer_beast_in.xml ${outputFile}_beast_in.xml
   popd
done < $listFile


