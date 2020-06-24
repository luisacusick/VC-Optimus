#!/bin/bash

printUsage(){
  echo "USAGE"
}

declare -a FASTQ
REF=''
DIV=0
INDEL=0

while getopts "hr:d:f:g:v:i:" option; do
  case ${option} in
  h) printUsage ;;
  r) REF=${OPTARG};; #reference
  d) DIV=${OPTARG};; #divergence, optional
  i) INDEL=${OPTARG};; #indel count, optional
  f) FASTQ1=${OPTARG};; #fastq file1
  g) FASTQ2=${OPTARG};; #fastq file2
  v) VCF=${OPTARG};; #vcf, optional
esac
done
shift "$((OPTIND -1))"

mkdir sim  

perl simuG/simuG.pl -refseq ${REF} -snp_count ${DIV} -prefix sim/mutated 

START=$(head -1 sim/mutated.simseq.genome.fa)

echo $START
echo $FASTQ1
echo $FASTQ2

time java -jar FastqGenerator/ArtificialFastqGenerator.jar -O sim/read -R sim/mutated.simseq.genome.fa -F1 ${FASTQ1} -F2 ${FASTQ2} -URQS true -SE true -S ${START} > sim/fastq.out 
