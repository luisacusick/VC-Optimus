#!/bin/bash

while getopts "hr:b:o:" option; do
  case ${option} in
  h) echo ""
     echo "Usage: runVCs.sh -r Ref.fa -b sample.bam -o outfolder" ;;
  r) REF=${OPTARG};;
  b) BAM=${OPTARG};;
  o) OUT=${OPTARG};;
  esac
done
shift "$(($OPTIND -1))"

mkdir ${OUT}

gatk/gatk HaplotypeCaller -R ${REF} -I ${BAM} -O ${OUT}/gatk.vcf

python createBed.py ${REF} ${OUT}/ref.bed

BED=${OUT}/ref.bed

vardict-java -U -G ${REF} -f .01 -N ${sampleName} -b ${BAM} -c 1 -S 2 -E 3 -g 4 ${BED} | teststrandbias.R | var2vcf_valid.pl -N ${sampleName} -E > ${OUT}/vardict.vcf

freebayes -f ${REF} ${BAM} > ${OUT}/freeBayes.vcf


