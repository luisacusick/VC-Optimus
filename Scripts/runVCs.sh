#!/bin/bash

GATK=false
VARDICT=false
FREEBAYES=false
OUT=''
BAM=''
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

while getopts "r:b:o:f:g:v:h" option; do
  case ${option} in
  h) echo ""
     echo "Usage: runVCs.sh -r Ref.fa -b sample.bam -o outfolder -f (includes free bayes) -g (includes gatk) -v (includes vardict)"
     exit 0;;
  r) REF=${OPTARG};;
  b) BAM=${OPTARG};;
  o) OUT=${OPTARG};;
  f) FREEBAYES=${OPTARG};;
  g) GATK=${OPTARG};;
  v) VARDICT=${OPTARG};;
  esac
done
shift "$(( OPTIND -1))"

mkdir ${OUT}/vcfs #create a vcf directory 

if [ "$GATK" = true ] ; then 
  gatk HaplotypeCaller -R ${REF} -I ${BAM} -O ${OUT}/vcfs/gatk.vcf
fi

if [ "$VARDICT" = true ] ; then
  python ${scriptDir}/createBed.py ${REF} ${OUT}/tmp/ref.bed
  sampleName="$(samtools view -H ${BAM} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)"
  BED=${OUT}/tmp/ref.bed
  vardict-java -U -G ${REF} -f .01 -b ${BAM} -c 1 -S 2 -E 3 -g 4 ${BED} -N ${sampleName} | teststrandbias.R | var2vcf_valid.pl -E -N ${sampleName} > ${OUT}/vcfs/vardict.vcf
fi

if [ "$FREEBAYES" = true ] ; then
  freebayes -f ${REF} ${BAM} > ${OUT}/vcfs/freeBayes.vcf
fi

