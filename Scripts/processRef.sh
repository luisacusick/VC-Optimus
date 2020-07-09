#!/bin/bash

while getopts "r" option; do
  case ${option} in
  h) echo 'Usage: processRef.sh -r <reference>'
     echo 'Creates an index, a sequence dictionary, and a bed file for the reference.'
  r) REF=${OPTARG};;
esac
done
shift "$((OPTIND -1))"

refDir=$((dirname "$REF"))


python createBed.py ${REF} ${refDir}/ref.bed

refPrefix=$(echo "${REF}" | cut -f 1 -d '.')

picard CreateSequenceDictionary R=${REF} O=${refPrefix}.dict #create ref dictionary

samtools faidx ${REF} #index the reference

bwa index ${REF} #index the reference for purposes of aligning 
