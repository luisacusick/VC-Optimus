#!/bin/bash

g=false
v=false
f=false
o=false

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" 

configFile=$(dirname $(dirname $(readlink -f "$0")))/config/paths.config # source config file with path to executables here
source ${configFile} 

if [[ ! -f ${VT_EXE} ]]
then
  echo "Error: cannot find the vt executable file; please specify its location in the header of this script, exiting.. "
  exit 0
fi

while getopts "c:r:s:d:o:g:v:f:h" option; do
  case ${option} in
  h) echo ""
     echo "Usage: normAndCombineVCF.sh [-d analysis_directory -r ref.fa -c ground_truth.vcf -s reference_sequence.dictionary -o output -f true -g true -v true]"
   	 echo ""
     echo "---Required---"
     echo "-r  [str] Path to reference genome fasta file"
     echo "-d  [str] Path to analysis directory; must contain directory 'vcfs'"
     echo "-s  [str] Path to reference sequence dictionary"
     echo "-o  [str] Path to output vcf file"
   	 echo ""
     echo "---Optional---"
     echo "-f  [bln] includes freebayes (default: false)"
     echo "-g  [bln] includes gatk (default: false)"
     echo "-v  [bln] includes vardict (default: false)"
     echo "-c  [str] Path to ground truth simulated vcf file"
   	 echo ""
     exit 0;;
  d) DIR=${OPTARG};; #VC-optimus output directory. must contain folder "vcfs" 
  r) REF=${OPTARG};;
  c) trueVCF=${OPTARG};;
  s) seqDict=${OPTARG};;
  g) g=${OPTARG};;
  v) v=${OPTARG};;
  f) f=${OPTARG};;
  o) one_vcf=${OPTARG};; #this creates a single normalized and merged vcf (i.e. doesn't create all combinations). should be false when performing simulations and likely true when a single results file is desired 
esac
done
shift "$((OPTIND -1))"


#echo "SEQ DICT ${seqDict}"
pathBase=${DIR}/vcfs

if [[ ! -f $trueVCF ]] && [[ ! -f ${trueVCF}.gz ]] #if there's no true vcf, or its gz compressed version, to compare results to, exit
then
  echo "$trueVCF, the ground truth vcf, does not exist, exiting.."
  exit 0
fi

#Don't have to normalize gatk vcfs because they are written in normalized form

if [[ $f -eq true ]]
then
  picard FixVcfHeader I=${pathBase}/freeBayes.vcf O=${pathBase}/freeBayes.fixed.vcf
  ${VT_EXE} normalize ${pathBase}/freeBayes.fixed.vcf -o ${pathBase}/freeBayesNorm.vcf -r ${REF}
  vcffilter -f "QUAL > 25" ${pathBase}/freeBayesNorm.vcf > ${pathBase}/freeBayesNormFilt.unsorted.vcf
  bcftools sort ${pathBase}/freeBayesNormFilt.unsorted.vcf -O v -o ${pathBase}/freeBayesNormFilt.vcf
  rm ${pathBase}/freeBayes.fixed.vcf*
  rm ${pathBase}/freeBayesNorm.vcf* #remove normalized unfiltered vcf  
  rm ${pathBase}/freeBayesNormFilt.unsorted.vcf*
fi

if [[ $v -eq true ]]
then
  picard FixVcfHeader I=${pathBase}/vardict.vcf O=${pathBase}/vardict.fixed.vcf
  ${VT_EXE} normalize ${pathBase}/vardict.fixed.vcf -o ${pathBase}/vardictNorm.vcf -r ${REF}
  vcffilter -f "QUAL > 25" ${pathBase}/vardictNorm.vcf > ${pathBase}/vardictNormFilt.unsorted.vcf
  bcftools sort ${pathBase}/vardictNormFilt.unsorted.vcf -O v -o ${pathBase}/vardictNormFilt.vcf
  rm ${pathBase}/vardict.fixed.vcf*
  rm ${pathBase}/vardictNorm.vcf*
  rm ${pathBase}/vardictNormFilt.unsorted.vcf*
fi

if [[ $g -eq true ]]
then 
  picard FixVcfHeader I=${pathBase}/gatk.vcf O=${pathBase}/gatk.fixed.vcf
  vcffilter -f "QUAL > 25" ${pathBase}/gatk.fixed.vcf > ${pathBase}/gatkFilt.unsorted.vcf
  bcftools sort ${pathBase}/gatkFilt.unsorted.vcf -O v -o ${pathBase}/gatkFilt.vcf
  rm ${pathBase}/gatk.fixed.vcf*
  rm ${pathBase}/gatkFilt.unsorted.vcf*
fi

skip=false
if [[ $f -eq true ]] && [[ $v -eq true ]] && [[ $g -eq true ]]
then
  #VCF with all 3 vc results
  picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/vardictNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/all.vcf COMMENT='freebayes, gatk and vardict'
  if [[ $o -eq true ]]
  then
    skip=true #switch that prevents us from creating every combination of vcf if we're just interested in one results file
  fi
fi

if [[ $f -eq true ]] && [[ $v -eq true ]] && [[ $skip -eq false ]]
then
  #VCF with free bayes and vardict
  picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/vardictNormFilt.vcf O=${pathBase}/fbVard.vcf COMMENT='freebayes and vardict'
fi

if [[ $f -eq true ]] && [[ $g -eq true ]] && [[ $skip -eq false ]]
then
  #VCF with free bayes and gatk
  picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/fbGatk.vcf D=${seqDict} COMMENT='freebayes and gatk'
fi

if [[ $v -eq true ]] && [[ $g -eq true ]] && [[ $skip -eq false ]]
then
  #VCF with vardict and gatk
  picard MergeVcfs I=${pathBase}/vardictNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/vardGatk.vcf COMMENT='vardict and gatk' D=${seqDict} #must specify the sequence dictionary when combining vardict and gatk vcfs (neither contain that info.)
fi

#process merged VCFs and compare to ground truth VCF

if [[ ! -f ${trueVCF}.gz ]] #check if a .gz file already exists or not
then
  bgzip ${trueVCF} #will remove original vcf file
fi

bcftools index ${trueVCF}.gz # index the true vcf

if [[ $f -eq true ]] && [[ $v -eq true ]] && [[ $g -eq true ]]
then
  bgzip ${pathBase}/all.vcf
  bcftools index ${pathBase}/all.vcf.gz
  bcftools stats ${pathBase}/all.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/all.log
fi

if [[ $f -eq true ]] && [[ $v -eq true ]]
then
  bgzip ${pathBase}/fbVard.vcf
  bcftools index ${pathBase}/fbVard.vcf.gz
  bcftools stats ${pathBase}/fbVard.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/fbVard.log
fi

if [[ $f -eq true ]] && [[ $g -eq true ]]
then
  bgzip ${pathBase}/fbGatk.vcf
  bcftools index ${pathBase}/fbGatk.vcf.gz
  bcftools stats ${pathBase}/fbGatk.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/fbGatk.log
fi

if [[ $v -eq true ]] && [[ $g -eq true ]]
then
  bgzip ${pathBase}/vardGatk.vcf
  bcftools index ${pathBase}/vardGatk.vcf.gz
  bcftools stats ${pathBase}/vardGatk.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/vardGatk.log
fi

if [[ $g -eq true ]]
then
  bgzip ${pathBase}/gatkFilt.vcf 
  bcftools index ${pathBase}/gatkFilt.vcf.gz
  bcftools stats ${pathBase}/gatkFilt.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/gatk.log
fi

if [[ $v -eq true ]]
then
  bgzip ${pathBase}/vardictNormFilt.vcf
  bcftools index -f ${pathBase}/vardictNormFilt.vcf.gz
  bcftools stats ${pathBase}/vardictNormFilt.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/vardict.log
fi

if [[ $f -eq true ]]
then
  bgzip ${pathBase}/freeBayesNormFilt.vcf
  bcftools index ${pathBase}/freeBayesNormFilt.vcf.gz
  bcftools stats ${pathBase}/freeBayesNormFilt.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/fb.log
fi

python ${scriptDir}/parseDiff.py ${pathBase} ${DIR}/run_summaries.csv
