#!/bin/bash

GATK=false
VARDICT=false
FREEBAYES=false

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

while getopts "hr:b:o:afgv" option; do
  case ${option} in
  h) echo ""
     echo "Usage: runVCs.sh -r Ref.fa -b sample.bam -o outfolder -f (excludes free bayes) -g (excludes gatk) -v (excludes vardict)"
     echo "The default is to run gatk, freebayes, and vardict. Use flags to exclude any of them from the analysis."
     exit 0;;
  r) REF=${OPTARG};;
  b) BAM=${OPTARG};;
  o) OUT=${OPTARG};;
  f) FREEBAYES=true;;
  g) GATK=true;;
  v) VARDICT=true;;
  esac
done
shift "$(($OPTIND -1))"

mkdir ${OUT}

if [ "$GATK" = true ] ; then 
  gatk HaplotypeCaller -R ${REF} -I ${BAM} -O ${OUT}/gatk.vcf
fi

if [ "$VARDICT" = true ] ; then
  python ${scriptDir}/createBed.py ${REF} ${OUT}/ref.bed
  bamPrefix=$(echo "${BAM}" | cut -f 1 -d '.')
  BED=${OUT}/ref.bed
  vardict-java -U -G ${REF} -f .01 -N ${bamPrefix} -b ${BAM} -c 1 -S 2 -E 3 -g 4 ${BED} | teststrandbias.R | var2vcf_valid.pl -N ${bamPrefix} -E > ${OUT}/vardict.vcf
fi

if [ "$FREEBAYES" = true ] ; then
  freebayes -f ${REF} ${BAM} > ${OUT}/freeBayes.vcf
fi

