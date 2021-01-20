#!/bin/bash

GATK=false
VARDICT=false
FREEBAYES=false
OUT=''
BAM=''
THREADS=1
TRUE_VCF=false # false by default, set to path to gold-standard VCF for comparison runs

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

configFile=$(dirname $(dirname $(readlink -f "$0")))/config/paths.config # file with paths to executables
source ${configFile}

while getopts "c:r:b:o:t:-:n:hgvf" option; do
  if [ "$option" = "-" ]; then
    option="${OPTARG%%=*}"
    OPTARG="${OPTARG#$option}"
    OPTARG="${OPTARG#=}" # if long option argument, remove assigning `=`
  fi
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
     echo "-f  [no arg] includes freebayes (default: false)"
     echo "-g  [no arg] includes gatk (default: false)"
     echo "-v  [no arg] includes vardict (default: false)"
     echo "-t  [int] Number of threads/processors to pass to multi-threaded programs (default: 1)"
   	 echo ""
     exit 0;;
  r) REF=${OPTARG};;
  b) BAM=${OPTARG};;
  o) OUT=${OPTARG};;
  t) THREADS=${OPTARG};;
  f) FREEBAYES=true;;
  g) GATK=true;;
  v) VARDICT=true;;
  n) SAMPLE_NAME=${OPTARG};;
  c) TRUE_VCF=${OPTARG};; 
  gatk) GATK_PARAMS=${OPTARG};;
  freebayes) FB_PARAMS=${OPTARG};;
  vardict) VARDICT_PARAMS=${OPTARG};;
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
     ${GATK_EXE} HaplotypeCaller -R ${REF} -I ${BAM} -O ${OUT}/vcfs/gatk.vcf --native-pair-hmm-threads ${THREADS} ${GATK_PARAMS}
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
     ${FREEBAYES_EXE} -f ${REF} ${BAM} ${FB_PARAMS} > ${OUT}/vcfs/freeBayes.vcf
  else
    echo "Skipping freebayes variant calling step, ${OUT}/vcfs/freeBayes.vcf already exists.."
  fi
fi

#combine VCF files from all simulated reads
echo "Combining VCF files from simulated read variant calling using normAndCombineVCF.sh.."
DICT=$(echo "${REF%.*}").dict 

#TODO: modify call to normAndCombine to work for simulateSample calls to runVCs and cmd line calls to runVCs

${scriptDir}/normAndCombineVCF.sh -d ${OUT} -r ${REF} -c ${TRUE_VCF} -s ${DICT} -o ${OUT}/${SAMPLE_NAME}.result_summary.txt -g true -v true -f true 1>>${OUT}/${SAMPLE_NAME}.simulateSample.log 2>>${OUT}/${SAMPLE_NAME}.simulateSample.err

