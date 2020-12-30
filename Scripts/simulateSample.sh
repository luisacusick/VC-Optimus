#!/bin/bash

DIV=0.01
INDEL=0
OUTDIR='vc-optimus-output/'
THREADS=1
SIMUG_EXE='/users/PAS1046/osu9029/analyses/optimus/reads-to-variants/Scripts/simuG.pl' #indicate path to simuG.pl here
AFQG_EXE='/users/PAS1046/osu9029/analyses/optimus/reads-to-variants/Scripts/ArtificialFastqGenerator.jar' #indicate path to ArtificialFastqGenerator.jar here
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" #all VC-optimus scripts live in the same directory; save the path to call the others

if [[ ! -f ${SIMUG_EXE} ]]
then
  echo "Error: cannot find the simuG.pl executable file; please specify its location in the header of this script, exiting.. "
  exit 0
fi

if [[ ! -f ${AFQG_EXE} ]]
then
  echo "Error: cannot find the ArtificialFastqGenerator.jar executable file; please specify its location in the header of this script, exiting.. "
  exit 0
fi


while getopts "n:r:s:t:o:d:i:l:h" option; do
  case ${option} in
   h)echo ""
     echo "Usage: simulateSample.sh -n name -r ref.fa -s read1.fq -s read2.fq [OPTIONS]"
   	 echo ""
     echo "---Required---"
     echo "-n  [str] Unique name to prefix simulateSample.sh output files and directories"
     echo "-r  [str] Path to reference genome fasta file to be used as starting point for simulations"
     echo "-s  [str] Path(s) to reads in fastq format for estimating sequencing error profiles. MUST be unzipped."
   	 echo ""
     echo "---Optional---"
     echo "-l  [int] Length of simulated sequence reads (default: the mode read length from first -s fastq file)"     
     echo "-t  [int] Number of threads/processors to pass to multi-threaded programs (default: 1)"
     echo "-o  [str] Path to output directory, will be created if it doesn't exist (default: vc-optimus-output/ in working directory)"
	 echo "-d  [flt] Amount of SNPs to simulate, as a fraction of total genome size (default: 0.01, equivalent to 1% sequence divergence)"
	 echo "-i  [int] Number of indels to simulate (NOT IMPLEMENTED; default: 0)"
   	 echo ""
     exit 0;;
  n) SAMPLE_NAME=${OPTARG};;
  r) REF=${OPTARG};;
  s) SAMPLES+=("$OPTARG");;
  l) SIMREAD_LENGTH=${OPTARG};;
  p) PAIRED=true;;
  t) THREADS=${OPTARG};;
  o) OUTDIR=${OPTARG};;
  d) DIV=${OPTARG};;
  i) INDEL=${OPTARG};;
esac
done
shift "$((OPTIND -1))"

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
  if [ ${SAMPLES[0]: -3} == ".gz" ]
  then
    echo "${SAMPLES[0]} cannot be a zipped file, exiting.."
    exit 0
  fi
else
  echo "Couldn't find sample ${SAMPLES[0]}, exiting.."
  exit 0
fi

if [ ${#SAMPLES[@]} -gt 1 ]
then
  if [[ -f ${SAMPLES[1]} ]]
  then
    echo "Found sample reads ${SAMPLES[1]}"
	  if [ ${SAMPLES[1]: -3} == ".gz" ]
	  then
		echo "${SAMPLES[1]} cannot be a zipped file, exiting.."
		exit 0
	  fi
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

#Verify sample name is provided
if [ -z ${SAMPLE_NAME} ]
then
  echo "You must provide a unique sample name that will be used to prefix output files, exiting.."
  exit 0
fi

if [ ! -d ${OUTDIR}/simulatedSamples ]
then
  mkdir ${OUTDIR}/simulatedSamples/
  echo "Created directory ${OUTDIR}/simulatedSamples"
fi
simDir=${OUTDIR}/simulatedSamples/${SAMPLE_NAME}
if [ ! -d ${simDir} ]
then
  mkdir ${simDir}/
  echo "Created output directory directory ${simDir}"
fi

#run simuG.pl
if [ ! -f ${simDir}/${SAMPLE_NAME}.simseq.genome.fa ]
then
  #calculate the number of snps to simulate
  chars=($(wc -m ${REF})) #counts number of characters in reference genome
  snps=$(echo "(($chars*$DIV)+0.5)/1" | bc) #calculate number of snps, rounded to nearest integer
  echo "Reference genome size: ${chars}bp, generating ${snps} SNPs in simulated genome using simuG.pl.."
  ${SIMUG_EXE} -refseq ${REF} -snp_count ${snps} -prefix ${simDir}/${SAMPLE_NAME} 1>${simDir}/${SAMPLE_NAME}.simulateSample.log 2>${simDir}/${SAMPLE_NAME}.simulateSample.err 
else
  echo "Skipping simulate genome step with simuG.pl, ${simDir}/${SAMPLE_NAME}.simseq.genome.fa already exists.."
fi

#run ArtificialFastqSimulator
if [ ! -f ${simDir}/${SAMPLE_NAME}.simseq.reads.1.fastq ]
then
  echo "Running ArtificialFastqSimulator.jar.."
  if [ -z $SIMREAD_LENGTH ]
  then
	SIMREAD_LENGTH=$(cat ${SAMPLES[0]} | awk '{if(NR%4==2) print length($1)}' | uniq -c | perl -p -e 's/^ +//' | sort -k1,1 -n | cut -f2 -d ' ' | tail -1) #calculate the most frequent read length
    echo "Setting simulated read length to $SIMREAD_LENGTH, the most frequent read length in ${SAMPLES[0]}"
  else
	echo "Simulated read length specified as $SIMREAD_LENGTH"
  fi
  START=$(head -1 ${simDir}/${SAMPLE_NAME}.simseq.genome.fa)
  java -jar ${AFQG_EXE} -O ${simDir}/${SAMPLE_NAME}.simseq.reads -R ${simDir}/${SAMPLE_NAME}.simseq.genome.fa -F1 ${SAMPLES[0]} -F2 ${SAMPLES[1]} -URQS true -SE true -S ${START} -RL ${SIMREAD_LENGTH} 1>>${simDir}/${SAMPLE_NAME}.simulateSample.log 2>>${simDir}/${SAMPLE_NAME}.simulateSample.err
else
  echo "Skipping ArtificialFastqSimulator.jar step, ${simDir}/${SAMPLE_NAME}.simseq.reads.1.fastq already exists.."
fi

#process reads using processSample.sh
if [ ! -f ${simDir}/processedSamples/${SAMPLE_NAME}.alnFinal.bam.bai ]
then
  echo "Processing simulated read samples using processSample.sh.."
  ${scriptDir}/processSample.sh -n ${SAMPLE_NAME} -r ${REF} -s ${simDir}/${SAMPLE_NAME}.simseq.reads.1.fastq -s ${simDir}/${SAMPLE_NAME}.simseq.reads.2.fastq -o ${simDir} -t ${THREADS}
else
  echo "Skipping processSample.sh step on simulated samples, ${simDir}/processedSamples/${SAMPLE_NAME}.alnFinal.bam.bai already exists.."
fi

#run all variant callers on simulated reads; runVC.sh will check if .vcf files already exist for each variant caller, and will skip that step if they do
echo "Running GATK, Vardict, and FreeBayes variant callers on simulated read samples using runVCs.sh.."
${scriptDir}/runVCs.sh -t ${THREADS} -r ${REF} -b ${simDir}/processedSamples/${SAMPLE_NAME}.alnFinal.bam -o ${simDir} -g true -v true -f true 1>>${simDir}/${SAMPLE_NAME}.simulateSample.log 2>>${simDir}/${SAMPLE_NAME}.simulateSample.err

#combine VCF files from all simulated reads
echo "Combining VCF files from simulated read variant calling using normAndCombineVCF.sh.."
DICT=$(echo "${REF%.*}").dict 
${scriptDir}/normAndCombineVCF.sh -d ${simDir} -r ${REF} -c ${simDir}/${SAMPLE_NAME}.refseq2simseq.SNP.vcf -s ${DICT} -o ${simDir}/${SAMPLE_NAME}.result_summary.txt -g true -v true -f true 1>>${simDir}/${SAMPLE_NAME}.simulateSample.log 2>>${simDir}/${SAMPLE_NAME}.simulateSample.err

