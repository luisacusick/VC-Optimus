#!/bin/bash

GATK=false
VARDICT=false
FREEBAYES=false
OUT=''
BAM=''
THREADS=1
GATK_EXE='gatk' #indicate path to gatk here, does not come packaged with VC-optimus env
VARDICT_EXE='vardict-java' #indicate path to vardict here, should come pre-loaded in VC-optimus env
FREEBAYES_EXE='freebayes' #indicate path to freebayes exe here, should come pre-loaded in VC-optimus env
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

while getopts "r:b:o:f:g:v:t:h" option; do
  case ${option} in
  h) echo ""
     echo "Usage: runVCs.sh -r ref.fa -b sample.bam -o output.directory [OPTIONS]"
   	 echo ""
     echo "---Required---"
     echo "-r  [str] Path to reference genome fasta file"
     echo "-b  [str] Path to mapped sample.bam file"
     echo "-o  [str] Path to output directory"
   	 echo ""
     echo "---Optional---"
     echo "-f  [bln] includes freebayes (default: false)"
     echo "-g  [bln] includes gatk (default: false)"
     echo "-v  [bln] includes vardict (default: false)"
     echo "-t  [int] Number of threads/processors to pass to multi-threaded programs (default: 1)"
   	 echo ""
     exit 0;;
  r) REF=${OPTARG};;
  b) BAM=${OPTARG};;
  o) OUT=${OPTARG};;
  t) THREADS=${OPTARG};;
  f) FREEBAYES=${OPTARG};;
  g) GATK=${OPTARG};;
  v) VARDICT=${OPTARG};;
  esac
done
shift "$(( OPTIND -1))"

#check user input to ensure files exist
if [[ ! -f ${REF} ]]
then
  echo "Reference not found, exiting.."
  exit 0
fi

if [[ ! -f ${BAM} ]]
then
  echo "mapped sample.bam not found, exiting.."
  exit 0
fi

if [[ ! -d ${OUT} ]]
then
  echo "mapped sample.bam not found, exiting.."
  exit 0
fi


mkdir ${OUT}/vcfs #create a vcf directory 

if [ "$GATK" = true ] ; then 
  if [ ! -f ${OUT}/vcfs/gatk.vcf ]
  then
     ${GATK_EXE} HaplotypeCaller -R ${REF} -I ${BAM} -O ${OUT}/vcfs/gatk.vcf --native-pair-hmm-threads ${THREADS}
  else
    echo "Skipping gatk variant calling step, ${OUT}/vcfs/gatk.vcf already exists.."
  fi
fi

if [ "$VARDICT" = true ] ; then
  if [ ! -f ${OUT}/vcfs/vardict.vcf ]
  then
     python ${scriptDir}/createBed.py ${REF} ${OUT}/vcfs/ref.bed
     sampleName="$(samtools view -H ${BAM} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)"
     BED=${OUT}/vcfs/ref.bed
     ${VARDICT_EXE} -U -G ${REF} -f .01 -b ${BAM} -c 1 -S 2 -E 3 -g 4 ${BED} -N ${sampleName} -th ${THREADS} | teststrandbias.R | var2vcf_valid.pl -E -N ${sampleName} > ${OUT}/vcfs/vardict.temp.vcf
     gatk UpdateVCFSequenceDictionary -V ${OUT}/vcfs/vardict.temp.vcf -R ${REF} --output ${OUT}/vcfs/vardict.vcf #add sequence dictionary to vardict vcf header
     rm ${OUT}/vcfs/vardict.temp.vcf
 else
    echo "Skipping vardict variant calling step, ${OUT}/vcfs/vardict.vcf already exists.."
  fi
fi

if [ "$FREEBAYES" = true ] ; then
  if [ ! -f ${OUT}/vcfs/freeBayes.vcf ]
  then
     ${FREEBAYES_EXE} -f ${REF} ${BAM} > ${OUT}/vcfs/freeBayes.vcf
  else
    echo "Skipping freebayes variant calling step, ${OUT}/vcfs/freeBayes.vcf already exists.."
  fi
fi

