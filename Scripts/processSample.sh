#!/bin/bash

declare -a SAMPLES
PAIRED=false
ZIPPED=false
CLEAN_NAMES=false
SAMPLE_NAME='SAMPLE'
ID='id'
RGLB='lib1'
RGPL='blank'
RGPU='000'
OUTDIR='vc-optimus-output/'
THREADS=1

while getopts "hcpi:l:pl:pu:n:zr:s:o:" option; do
  case ${option} in
   h) 
     echo ""
     echo "Usage: processSample.sh [-c|-p|-z] -r Ref.fa -s read1.fq -s read2.fq"
     echo "-c True if paired end reads are followed by ".1" and ".2. "These must be removed align reads. Default: False"
     echo "-p Paired end reads? Default: False"
     echo "-z Zipped reads? Default: False"
     echo "-n Sample name. Default: Sample"
     echo "-r Reference to align reads to (Required) "
     echo "-s Reads in fastq format. (Required)"
     echo ""
     echo "---Read Group Naming Options (GATK requires named read groups)"
     echo "-i Read group id. Default: id"
     echo "-l Read group library name. Default: lib1"
     echo "-pl Read group platform. Default: blank"
     echo "-pu Read group platform unit (ex. SRRXXXXXXX) Default: 000"
     echo "-o Path to output directory. If it doesn't already exist program creates it."
     echo "-t Number of threads for aligner"
     exit 0;;
  c) CLEAN_NAMES=true;;
  p) PAIRED=true;;
  i) ID=${OPTARG};;
  l) RGLB=${OPTARG};;
  pl) RGPL=${OPTARG};;
  pu) RGPU=${OPTARG};;
  n) SAMPLE_NAME=${OPTARG};;
  z) ZIPPED=true;;
  r) REF=${OPTARG};;
  s) SAMPLES+=("$OPTARG");;
  o) OUTDIR=${OPTARG};;
  t) THREADS=${OPTARG};;
esac
done
shift "$(($OPTIND -1))"

echo "REF ${REF}"
if [[ -f ${REF} ]]
then
  echo "Found ${REF}"
else 
  echo "Reference not found"
  exit 0
fi

if [[ -f ${SAMPLES[0]} ]]
then
  echo "Found ${SAMPLES[0]}"
else
  echo "Couldn't find sample ${SAMPLES[0]}"
  exit 0
fi

if [ ${#SAMPLES[@]} -gt 1 ]
then
  if [[ -f ${SAMPLES[1]} ]]
  then
    echo "Found ${SAMPLES[1]}"
  else
    echo "Couldn't find second sample ${SAMPLES[1]}"
    exit 0
  fi
fi

if [ ${#SAMPLES[@]} -gt 2 ]
then 
  echo "Use a maximum of 2 samples"
  exit 0
fi

#Verify thread number is > 1
if [ ${THREADS} -lt 1 ]
then
  echo "Thread must be an int greater than or equal to 1"
  exit 0
fi

OUTPATH=${OUTDIR}/tmp

mkdir ${OUTDIR}/tmp #store non-result intermmediate files here

if ${CLEAN_NAMES} #If the reads need ".1" and ".2" removed from the ends of their names
then
  for sample in "${SAMPLES}"
  do
    gzip -c -d $i | sed -E 's/(^[@+]SRR[0-9]+\.[0-9]+)\.[12]/\1/' | gzip -c > "${OUTPATH}/${sample}"
  done
fi

#Align reads to reference#

#Align sample to ref
bwa mem ${REF} ${SAMPLES[0]} ${SAMPLES[1]} -t ${THREADS} -M -o ${OUTPATH}/aln.sam
samtools sort ${OUTPATH}/aln.sam -o ${OUTPATH}/aln.bam -@ 8

rm ${OUTPATH}/aln.sam #remove sam file 

#Insert read group naming
picard AddOrReplaceReadGroups INPUT=${OUTPATH}/aln.bam OUTPUT=${OUTPATH}/alnRG.bam RGID=${ID} RGLB=${RGLB} RGPL=${RGPL} RGPU=${RGPU} RGSM=${SAMPLE_NAME}

#Mark duplicates
picard MarkDuplicates I=${OUTPATH}/alnRG.bam O=${OUTPATH}/alnFinal.bam M=${OUTPATH}/marked_dup_metrics.txt

echo "Created alnFinal.bam in output tmp directory"

samtools index ${OUTPATH}/alnFinal.bam #creates .bam.bai (indexed bam) in input directory

picard ValidateSamFile I=${OUTPATH}/alnFinal.bam R=${REF} MODE=SUMMARY
