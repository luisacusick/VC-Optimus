#!/bin/bash

OUTDIR='vc-optimus-output/' # set default output directory, user input can override 

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

while getopts ":hr:p:o:" option; do
  case ${option} in
  h) echo ""
     echo 'Usage: processRef.sh -r reference.fa'
   	 echo ""
     echo 'Creates a sequence dictionary and indexes for the reference'
   	 echo ""
     echo "---Required---"
     echo "-r  [str] Path to reference genome fasta file"
   	 echo ""
     echo "---Optional---"
     exit 0;;
  r) REF=${OPTARG};;
  p) PARAM=${OPTARG};;
  o) OUTDIR=${OPTARG};;
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

if [[ ! -d ${OUTDIR} ]]
then
  mkdir ${OUTDIR} # if output directory does not exist create it
fi

mkdir ${OUTDIR}/logs # create log directry

if [[ ! -f ${refPrefix}.dict ]]
then
  picard CreateSequenceDictionary R=${REF} O=${refPrefix}.dict 1>> ${OUTDIR}/logs/processRef.log 2>> ${OUTDIR}/logs/processRef.err #create ref dictionary
fi

if [[ ! -f ${refPrefix}.fai ]]
then
  samtools faidx ${REF} 1>> ${OUTDIR}/logs/processRef.log 2>> ${OUTDIR}/logs/processRef.err #index the reference
fi

bwa index ${REF} 1>> ${OUTDIR}/logs/processRef.log 2>> ${OUTDIR}/logs/processRef.err #index the reference for purposes of aligning 
