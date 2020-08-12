#!/bin/bash

while getopts "r" option; do
  case ${option} in
  h) echo 'Usage: processRef.sh -r <reference>'
     echo 'Creates an index, a sequence dictionary, and a bed file for the reference.'
     exit 0 ;;
  r) REF=${OPTARG};;
esac
done
shift "$((OPTIND -1))"

refDir=$((dirname "$REF"))

if [[ -f ${REF} ]]
then
  echo "Found ${REF}"
else 
  echo "Reference not found"
  exit 0
fi

python createBed.py ${REF} ${refDir}/ref.bed >> out 2>&1

refPrefix=$(echo "${REF%.*}")

picard CreateSequenceDictionary R=${REF} O=${refPrefix}.dict #create ref dictionary

samtools faidx ${REF} #index the reference

bwa index ${REF} #index the reference for purposes of aligning 
