#!/bin/bash

GATK=true
VARDICT=true
FREEBAYES=true

while getopts "hr:b:o:afgv" option; do
  case ${option} in
  h) echo ""
     echo "Usage: runVCs.sh -r Ref.fa -b sample.bam -o outfolder -f (excludes free bayes) -g (excludes gatk) -v (excludes vardict)"
     echo "The default is to run gatk, free bayes, and vardict. Use flags to exclude any of them from the analysis."
  r) REF=${OPTARG};;
  b) BAM=${OPTARG};;
  o) OUT=${OPTARG};;
  f) FREEBAYES=false;;
  g) GATK=false;;
  v) VARDICT=false;;
  esac
done
shift "$(($OPTIND -1))"

mkdir ${OUT}


if [ "$GATK" = true ] ; then 
  gatk/gatk HaplotypeCaller -R ${REF} -I ${BAM} -O ${OUT}/gatk.vcf
fi

if [ "$VARDICT" = true ] ; then
  python createBed.py ${REF} ${OUT}/ref.bed

  BED=${OUT}/ref.bed

  vardict-java -U -G ${REF} -f .01 -N ${sampleName} -b ${BAM} -c 1 -S 2 -E 3 -g 4 ${BED} | teststrandbias.R | var2vcf_valid.pl -N ${sampleName} -E > ${OUT}/vardict.vcf
fi

if [ "$FREEBAYES" = true ] ; then
  freebayes -f ${REF} ${BAM} > ${OUT}/freeBayes.vcf
fi

