#!/bin/bash

declare -a SAMPLES
CLEAN_NAMES=false
ID='id'
RGLB='lib1'
RGPL='blank'
RGPU='SRR000'
OUTDIR='vc-optimus-output/'
THREADS=1

while getopts "h:n:r:s:c:t:o:i:l:pl:pu:" option; do
  case ${option} in
   h) 
     echo ""
     echo "Usage: processSample.sh -n name -r ref.fa -s read1.fq -s read2.fq [OPTIONS]"
   	 echo ""
     echo "---Required---"
     echo "-n  [str] Unique name to prefix sampleProcess.sh output files"
     echo "-r  [str] Path to reference genome fasta file"
     echo "-s  [str] Path(s) to reads in fastq format for mapping to reference"
   	 echo ""
     echo "---Optional---"
     echo "-c        Indicates paired end reads are followed by ".1" and ".2. "These must be removed in order to properly align reads (default: False)"
     echo "-t  [int] Number of threads/processors to pass to multi-threaded programs (default: 1)"
     echo "-o  [str] Path to output directory, will be created if it doesn't exist (default: vc-optimus-output/ in working directory)"
   	 echo ""
     echo "---Read Group Naming Options (GATK requires named read groups)---"
     echo "-i  [str] Read group id (default: id)"
     echo "-l  [str] Read group library name (default: lib1)"
     echo "-pl [str] Read group platform (default: blank)"
     echo "-pu [str] Read group platform unit (ex. SRRXXXXXX; default: SRR000)"
   	 echo ""
     exit 0;;
  n) SAMPLE_NAME=${OPTARG};;
  r) REF=${OPTARG};;
  s) SAMPLES+=("$OPTARG");;
  c) CLEAN_NAMES=true;;
  t) THREADS=${OPTARG};;
  o) OUTDIR=${OPTARG};;
  i) ID=${OPTARG};;
  l) RGLB=${OPTARG};;
  pl) RGPL=${OPTARG};;
  pu) RGPU=${OPTARG};;
esac
done
shift "$(($OPTIND -1))"

#check user input to ensure files exist
if [[ -f ${REF} ]]
then
  echo "Found reference ${REF}"
else 
  echo "Reference not found, exiting.."
  exit 0
fi

if [[ -f ${SAMPLES[0]} ]]
then
  echo "Found sample reads ${SAMPLES[0]}"
else
  echo "Couldn't find sample ${SAMPLES[0]}, exiting.."
  exit 0
fi

if [ ${#SAMPLES[@]} -gt 1 ]
then
  if [[ -f ${SAMPLES[1]} ]]
  then
    echo "Found sample reads ${SAMPLES[1]}"
  else
    echo "Couldn't find second sample ${SAMPLES[1]}, exiting.."
    exit 0
  fi
fi

if [ ${#SAMPLES[@]} -gt 2 ]
then 
  echo "You can only have a maximum of 2 read samples, exiting.."
  exit 0
fi

#Verify thread number is > 1
if [ ${THREADS} -lt 1 ]
then
  echo "Thread must be an integer greater than or equal to 1, exiting.."
  exit 0
fi

#Verify sample name is provided
if [ -z ${SAMPLE_NAME} ]
then
  echo "You must provide a unique sample name that will be used to prefix output files, exiting.."
  exit 0
fi

#create processedSamples in provided output directory
mkdir ${OUTDIR}/
OUTPATH=${OUTDIR}/processedSamples
mkdir ${OUTDIR}/processedSamples/ #store non-result intermmediate files here

if ${CLEAN_NAMES} #If the reads need ".1" and ".2" removed from the ends of their names
then
  for sample in "${SAMPLES}"
  do
    gzip -c -d $i | sed -E 's/(^[@+]SRR[0-9]+\.[0-9]+)\.[12]/\1/' | gzip -c > "${OUTPATH}/${sample}"
  done
fi

#Align reads to reference with bwa mem and sort with samtools#
bwa mem ${REF} ${SAMPLES[0]} ${SAMPLES[1]} -t ${THREADS} -M -o ${OUTPATH}/${SAMPLE_NAME}.aln.sam 1> ${OUTPATH}/${SAMPLE_NAME}.processSample.log 2> ${OUTPATH}/${SAMPLE_NAME}.processSample.err
samtools sort ${OUTPATH}/${SAMPLE_NAME}.aln.sam -o ${OUTPATH}/${SAMPLE_NAME}.aln.bam -@ ${THREADS} 1>> ${OUTPATH}/${SAMPLE_NAME}.processSample.log 2>> ${OUTPATH}/${SAMPLE_NAME}.processSample.err
rm ${OUTPATH}/${SAMPLE_NAME}.aln.sam #remove sam file 

#Insert read group naming with picard
picard AddOrReplaceReadGroups INPUT=${OUTPATH}/${SAMPLE_NAME}.aln.bam OUTPUT=${OUTPATH}/${SAMPLE_NAME}.alnRG.bam RGID=${ID} RGLB=${RGLB} RGPL=${RGPL} RGPU=${RGPU} RGSM=${SAMPLE_NAME} 1>> ${OUTPATH}/${SAMPLE_NAME}.processSample.log 2>> ${OUTPATH}/${SAMPLE_NAME}.processSample.err

#Mark duplicates with picard
picard MarkDuplicates I=${OUTPATH}/${SAMPLE_NAME}.alnRG.bam O=${OUTPATH}/${SAMPLE_NAME}.alnFinal.bam M=${OUTPATH}/${SAMPLE_NAME}.marked_dup_metrics.txt 1>> ${OUTPATH}/${SAMPLE_NAME}.processSample.log 2>> ${OUTPATH}/${SAMPLE_NAME}.processSample.err

echo "Created ${SAMPLE_NAME}.alnFinal.bam in processedSamples results directory"

#index bam and validate with picard
samtools index ${OUTPATH}/${SAMPLE_NAME}.alnFinal.bam 1>> ${OUTPATH}/${SAMPLE_NAME}.processSample.log 2>> ${OUTPATH}/${SAMPLE_NAME}.processSample.err #creates .bam.bai (indexed bam) in input directory
picard ValidateSamFile I=${OUTPATH}/${SAMPLE_NAME}.alnFinal.bam R=${REF} MODE=SUMMARY 1>> ${OUTPATH}/${SAMPLE_NAME}.processSample.log 2>> ${OUTPATH}/${SAMPLE_NAME}.processSample.err
