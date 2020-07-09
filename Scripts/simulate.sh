#!/bin/bash

printUsage(){
  echo "USAGE"
}

declare -a FASTQ
REF=''
DIV=0.1
INDEL=0

while getopts "hr:d:s:g:v:i:" option; do
  case ${option} in
  h) printUsage ;;
  r) REF=${OPTARG};; #reference
  d) DIV=${OPTARG};; #divergence, optional
  i) INDEL=${OPTARG};; #indel count, optional
  s) FASTQ+=("$OPTARG");; #fastq files
  v) VCF=${OPTARG};; #vcf, optional
esac
done
shift "$((OPTIND -1))"

chars=($(wc -m ${REF})) #counts number of characters in reference genome
snps=$(echo "$chars*$DIV" | bc) #calculate number of snps 

mkdir sim  

perl simuG/simuG.pl -refseq ${REF} -snp_count ${snps} -prefix sim/mutated 

START=$(head -1 sim/mutated.simseq.genome.fa)

java -jar FastqGenerator/ArtificialFastqGenerator.jar -O sim/read -R sim/mutated.simseq.genome.fa -F1 ${SAMPLES[0]} -F2 ${SAMPLES[1]} -URQS true -SE true -S ${START} > sim/fastq.out

./processSample.sh -r ${REF} -s sim/read.1.fastq -s sim/read.2.fastq -o sim

./runVCs.sh -r ${REF} 
