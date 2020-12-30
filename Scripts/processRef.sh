#!/bin/bash

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

while getopts ":hr:p:" option; do
  case ${option} in
  h) echo ""
     echo 'Usage: processRef.sh -r reference.fa'
   	 echo ""
     echo 'Creates a sequence dictionary and indexes for the reference'
   	 echo ""
     echo "---Required---"
     echo "-r  [str] Path to reference genome fasta file"
   	 echo ""
     exit 0;;
  r) REF=${OPTARG};;
  p) PARAM=${OPTARG};;
  esac
done
shift "$((OPTIND -1))" 

if [[ -f ${REF} ]]
then
  echo "Found ${REF}"
else
  echo "Reference not found"
  exit 0
fi

#set prefix as all chars up to the terminal '.suffix'
refDir="$(dirname $REF)"
refFilename=$(basename -- "${REF}")
refPrefix="${filename%.*}"

if [[ ! -f ${refPrefix}.dict ]]
then
  picard CreateSequenceDictionary R=${REF} O=${refPrefix}.dict #create ref dictionary
fi

if [[ ! -f ${refPrefix}.fai ]]
then
  samtools faidx ${REF} #index the reference
fi

bwa index ${REF} #index the reference for purposes of aligning 
