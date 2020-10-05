#!/bin/bash

printUsage(){
  echo "USAGE"
  exit 0
}

REF=''
DIV=0.1
INDEL=0
PAIRED=false
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

while getopts "hpr:d:s:t:g:v:i:o:" option; do
  case ${option} in
  h) printUsage;;
  p) PAIRED=true;;
  r) REF=${OPTARG};; #reference
  d) DIV=${OPTARG};; #divergence, optional
  i) INDEL=${OPTARG};; #indel count, optional
  s) SAMPLE1=${OPTARG};; #fastq file 
  t) SAMPLE2=${OPTARG};; #optional second fastq file
  v) VCF=${OPTARG};; #vcf, optional
  o) output=${OPTARG};; #output directory, optional
esac
done
shift "$((OPTIND -1))"

#check user input to ensure files exist
if [[ -f ${REF} ]]
then
  echo "Found ${REF}"
else
  echo "Reference sequence not found"
  exit 0
fi

if [[ -f ${SAMPLE1} ]]
then
  echo "Found ${SAMPLE1}"
else
  echo "Couldn't find sample"
  exit 0
fi

if [ -n ${SAMPLE2} ]
then
  if [[ -f ${SAMPLE2} ]]
  then
    echo "Found ${SAMPLE2}"
  else
    echo "Couldn't find second sample"
    exit 0
  fi
fi

chars=($(wc -m ${REF})) #counts number of characters in reference genome
snps=$(echo "(($chars*$DIV)+0.5)/1" | bc) #calculate number of snps, rounded to nearest integer

simDir=${output}/sim.$(date "+%Y.%m.%d-%H.%M.%S")
mkdir ${simDir}

refFile=$(basename -- "$REF")
refPrefix="${refFile%.*}"

simuG.pl -refseq ${REF} -snp_count ${snps} -prefix ${simDir}/${refPrefix} 1> ${simDir}/simuG.out 2> ${simDir}/simuG.err 

START=$(head -1 ${simDir}/${refPrefix}.simseq.genome.fa)

java -jar ArtificialFastqGenerator.jar -O ${simDir}/${refPrefix}.simseq.reads -R ${simDir}/${refPrefix}.simseq.genome.fa -F1 ${SAMPLE1} -F2 ${SAMPLE2} -URQS true -SE true -S ${START}

${scriptDir}/processSample.sh -r ${REF} -s ${simDir}/${refPrefix}.simseq.reads.1.fastq -s ${simDir}/${refPrefix}.simseq.reads.2.fastq -o ${simDir}

mkdir ${simDir}/vcfs
${scriptDir}/runVCs.sh -r ${REF} -b ${simDir}/alnFinal.bam -o ${simDir}/vcfs

DICT=$(echo "${REF%.*}").dict 

${scriptDir}/normAndCombineVCF.sh -d ${simDir} -r ${REF} -v ${simDir}/${refPrefix}.refseq2simseq.SNP.vcf -s ${DICT} -o ${simDir}/result_summary.txt

