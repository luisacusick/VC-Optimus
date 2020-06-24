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
OUTPATH='tmp'

while getopts "hcpi:l:pl:pu:n:zr:s:o:" option; do
  case ${option} in
   h) 
     echo ""
     echo "Usage: process.sh [-c|-p|-z] -r Ref.fa -s read1.fq -s read2.fq"
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
  o) OUTPATH=${OPTARG};;
esac
done
shift "$(($OPTIND -1))"

mkdir $OUTPATH #Temporary directory to keep intermediate files in 

if ${CLEAN_NAMES} #If the reads need ".1" and ".2" removed from the ends of their names
then
  for sample in "${SAMPLES}"
  do
    gzip -c -d $i | sed -E 's/(^[@+]SRR[0-9]+\.[0-9]+)\.[12]/\1/' | gzip -c > "${OUTPATH}/${sample}"
  done
fi

#Align reads to reference#

#Index the reference:
bwa index ${REF}
bwa mem ${REF} ${SAMPLES[1]} ${SAMPLES[2]} -t 8 -M | samtools sort -@8 -o ${OUTPATH}/aln.bam - 
#rm ${OUTPATH}/aln.sam #remove sam file 

#Insert read group naming
picard AddOrReplaceReadGroups INPUT=${OUTPATH}/aln.bam OUTPUT=${OUTPATH}/alnRG.bam RGID=${ID} RGLB=${RGLB} RGPL=${RGPL} RGPU=${RGPU} RGSM=${SAMPLE_NAME}

#Sort bam file
#samtools sort ${OUTPATH}/alnRG.bam -o ${OUTPATH}/alnSort.bam

#Mark duplicates
picard MarkDuplicates I=${OUTPATH}/alnRG.bam O=${OUTPATH}/alnFinal.bam M=${OUTPATH}/marked_dup_metrics.txt
samtools index ${OUTPATH}/alnFinal.bam #creates .bam.bai file in input directory

echo "Created alnFinal.bam in output directory"

#Recalibrate base scores
gatk/gatk BaseRecalibrator -I ${OUTPATH}/alnFinal.bam -R ${REF} -O ${OUTPATH}/recal_data.table

picard ValidateSamFile I=${OUTPATH}/alnFinal.bam MODE=SUMMARY

#These steps are necessary for each referenece for variant calling, but don't want to perform them more than once per reference
#PREFIX=$(echo "${REF}" | cut -f 1 -d '.')
#picard CreateSequenceDictionary R=${REF} O=${PREFIX}.dict #create reference dictionary
#samtools faidx ${REF} #index reference  
