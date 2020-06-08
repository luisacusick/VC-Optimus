#!/bin/bash

declare -a SAMPLES
PAIRED=false
ZIPPED=false
CLEAN_NAMES=false
SAMPLE_NAME='SAMPLE'

while getopts ":r:s:n:cpzh" option; do
  case ${option} in
  h | ?) 
     echo "Usage: program [-c|-p|-z] -r Ref.fa -s read1.fq read2.fq"
     echo "-c Paired end reads downloaded from fastq-dump have .1 or .2 appended to name. These must be removed to properly align reads. Default: False"
     echo "-p Paired end reads? Default: False"
     echo "-z Zipped reads? Default: False"
     echo "-n Sample name. Default: Sample"
     echo "-r Reference to align reads to (Required) "
     echo "-s Reads in fastq format. (Required)"
     exit 0;;
  c) CLEAN_NAMES=true;;
  p) PAIRED=true;;
  n) SAMPLE_NAME=${OPTARG};;
  z) ZIPPED=true;;
  r) REF=${OPTARG};;
  s) SAMPLES+=("$OPTARG");;
esac
done
shift "$(($OPTIND -1))"

mkdir tmp

if ${CLEAN_NAMES} #If the reads need ".1" and ".2" removed from the ends of their names
then
  for sample in "${SAMPLES}"
  do
    gzip -c -d $i | sed -E 's/(^[@+]SRR[0-9]+\.[0-9]+)\.[12]/\1/' | gzip -c > "tmp/${sample}"
  done
fi

#Align reads to reference#

#Index the reference:
bwa index ${REF}
bwa mem ${REF} ${SAMPLES[1]} ${SAMPLES[2]} -o tmp/aln.sam -R "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA"

#Process SAM file
samtools view -S -b tmp/aln.sam > tmp/aln.bam #convert sam to bam
samtools sort tmp/aln.bam -o tmp/alnSort.bam
samtools index tmp/alnSort.bam #creates .bam.bai file in input directory

echo "Created alnSort.bam in tmp directory"


