#!/bin/bash

printUsage(){
  echo "USAGE"
  exit 0
}

declare -a SAMPLES
REF=''
DIV=0.1
INDEL=0

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

while getopts "hr:d:s:g:v:i:" option; do
  case ${option} in
  h) printUsage ;;
  r) REF=${OPTARG};; #reference
  d) DIV=${OPTARG};; #divergence, optional
  i) INDEL=${OPTARG};; #indel count, optional
  s) SAMPLES+=("$OPTARG");; #fastq files
  v) VCF=${OPTARG};; #vcf, optional
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

if [[ -f ${SAMPLES[0]} ]]
then
  echo "Found ${SAMPLES[0]}"
else
  echo "Couldn't find sample"
  exit 0
fi

if [ ${#SAMPLES[@]} -gt 1 ]
then
  if [[ -f ${SAMPLES[1]} ]]
  then
    echo "Found ${SAMPLES[1]}"
  else
    echo "Couldn't find second sample"
    exit 0
  fi
fi

chars=($(wc -m ${REF})) #counts number of characters in reference genome
snps=$(echo "(($chars*$DIV)+0.5)/1" | bc) #calculate number of snps, rounded to nearest integer

simDir = sim.$(date "+%Y.%m.%d-%H.%M.%S") #a unique name with a timestamp to enable multiple simultation runs
mkdir ${simDir}  

refPrefix=$(echo "${REF}" | cut -f 1 -d '.')

perl ${scriptDir}/simuG.pl -refseq ${REF} -snp_count ${snps} -prefix ${simDir}/${refPrefix} 

START=$(head -1 ${simDir}/${refPrefix}.simseq.genome.fa)

java -jar ${scriptDir}/ArtificialFastqGenerator.jar -O ${simDir}/${refPrefix}.simseq.reads -R ${simDir}/${refPrefix}.simseq.genome.fa -F1 ${SAMPLES[0]} -F2 ${SAMPLES[1]} -URQS true -SE true -S ${START} 1> ${simDir}/fastq.out 2> ${simDir}/fastq.err

${scriptDir}/processSample.sh -r ${REF} -s ${simDir}/${refPrefix}.simseq.reads.1.fastq -s ${simDir}/${refPrefix}.simseq.reads.fastq -o ${simDir}

mkdir sim/results
${scriptDir}/runVCs.sh -r ${REF} -b ${simDir}/alnFinal.bam -o ${simDir}/results

DICT=$(echo "${REF%.*}").dict 

${scriptDir}/normAndCombine.sh -d ${simDir} -r ${REF} -v ${simDir}/${refPrefix}.refseq2simseq.SNP.vcf -s ${DICT}

