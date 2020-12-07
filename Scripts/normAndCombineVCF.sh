#!/bin/bash

g=false
v=false
f=false
o=false
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" 

while getopts "c:r:s:d:o:g:v:f:h" option; do
  case ${option} in
  h) "Usage: " ;;
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


echo "SEQ DICT ${seqDict}"
pathBase=${DIR}/vcfs

#Don't have to normalize gatk vcfs because they are written in normalized form

if [[ $f -eq true ]]
then
  vt normalize ${pathBase}/freeBayes.vcf -o ${pathBase}/freeBayesNorm.vcf -r ${REF}
  vcffilter -f "QUAL > 25" ${pathBase}/freeBayesNorm.vcf > ${pathBase}/freeBayesNormFilt.vcf
  rm ${pathBase}/freeBayesNorm.vcf #remove normalized unfiltered vcf  
fi

if [[ $v -eq true ]]
then
  vt normalize ${pathBase}/vardict.vcf -o ${pathBase}/vardictNorm.vcf -r ${REF}
  vcffilter -f "QUAL > 25" ${pathBase}/vardictNorm.vcf > ${pathBase}/vardictNormFilt.vcf
  rm ${pathBase}/vardictNorm.vcf
fi

if [[ $g -eq true ]]
then 
  vcffilter -f "QUAL > 25" ${pathBase}/gatk.vcf > ${pathBase}/gatkFilt.vcf
fi

skip=false
if [[ $f -eq true ]] && [[ $v -eq true ]] && [[ $g -eq true ]]
then
  #VCF with all 3 vc results
  picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/vardictNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/all.vcf
  if [[ $o -eq true ]]
  then
    skip=true #switch that prevents us from creating every combination of vcf if we're just interested in one results file
  fi
fi

if [[ $f -eq true ]] && [[ $v -eq true ]] && [[ $skip -eq false ]]
then
  #VCF with free bayes and vardict
  picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/vardictNormFilt.vcf O=${pathBase}/fbVard.vcf COMMENT='free bayes and vardict'
fi

if [[ $f -eq true ]] && [[ $g -eq true ]] && [[ $skip -eq false ]]
then
  #VCF with free bayes and gatk
  picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/fbGatk.vcf D=${seqDict} COMMENT='free bayes and gatk'
fi

if [[ $v -eq true ]] && [[ $g -eq true ]] && [[ $skip -eq false ]]
then
  #VCF with vardict and gatk
  picard MergeVcfs I=${pathBase}/vardictNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/vardGatk.vcf COMMENT='vardict and gatk' D=${seqDict} #must specify the sequence dictionary when combining vardict and gatk vcfs (neither contain that info.)
fi

if [[ ! -f $trueVCF ]] #if there's no true vcf to compare results to, exit
then
  echo $trueVCF 
  echo "exiting as there is no true VCF to compare to"
  exit 0
fi

bgzip ${trueVCF}
bcftools index ${trueVCF}.gz

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
