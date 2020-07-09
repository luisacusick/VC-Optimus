#!/bin/bash

#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -t 00:30:00
#SBATCH --output=/pylon5/eb5phrp/luc32/MedTrunSimulations_Full50/%a/combine.out
#SBATCH --array=0

module load anaconda2/5.2.0
source activate snvCalling

index=${SLURM_ARRAY_TASK_ID}
trueVCF=/pylon5/eb5phrp/luc32/MedTrunSimulations_Full50/${index}/${index}.refseq2simseq.SNP.vcf
pathBase=/pylon5/eb5phrp/luc32/MedTrunSimulations_Full50/${index}
seqDict=/pylon5/eb5phrp/luc32/Data/Reference/MedTrunRefShort.dict
ref=/pylon5/eb5phrp/luc32/Data/Reference/MedTrunRefShort.fasta

#Don't have to normalize gatk vcfs because they are written in normalized form
#vt/vt normalize ${pathBase}/freeBayes.vcf -o ${pathBase}/freeBayesNorm.vcf -r ${ref}

#vt/vt normalize ${pathBase}/vardict.vcf -o ${pathBase}/vardictNorm.vcf -r ${ref}

##FILTER RESULTS##
#vcffilter -f "QUAL > 10" ${pathBase}/freeBayesNorm.vcf > ${pathBase}/freeBayesNormFilt.vcf
#vcffilter -f "QUAL > 10" ${pathBase}/vardictNorm.vcf > ${pathBase}/vardictNormFilt.vcf
#vcffilter -f "QUAL > 10" ${pathBase}/gatk2.vcf > ${pathBase}/gatkFilt.vcf

#rm ${pathBase}/freeBayesNorm.vcf #removed normalized but unfiltered vcfs
#rm ${pathBase}/vardictNorm.vcf

#VCF with all 3 vc results
#picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/vardictNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/all.vcf

#VCF with free bayes and vardict
#picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/vardictNormFilt.vcf O=${pathBase}/fbVard.vcf COMMENT='free bayes and vardict'

#VCF with free bayes and gatk
#picard MergeVcfs I=${pathBase}/freeBayesNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/fbGatk.vcf COMMENT='free bayes and gatk'

#VCF with vardict and gatk
#picard MergeVcfs I=${pathBase}/vardictNormFilt.vcf I=${pathBase}/gatkFilt.vcf O=${pathBase}/vardGatk.vcf COMMENT='vardict and gatk' D=${seqDict} #must specify the sequence dictionary when combining vardict and gatk vcfs (neither contain that info.)

#bgzip ${pathBase}/all.vcf
#bcftools index ${pathBase}/all.vcf.gz
#bgzip ${pathBase}/fbVard.vcf
#bcftools index ${pathBase}/fbVard.vcf.gz
#bgzip ${pathBase}/fbGatk.vcf
#bcftools index ${pathBase}/fbGatk.vcf.gz
#bgzip ${pathBase}/vardGatk.vcf
#bcftools index ${pathBase}/vardGatk.vcf.gz
bgzip ${pathBase}/gatkFilt.vcf 
#bcftools index ${pathBase}/gatkFilt.vcf.gz
bgzip ${pathBase}/vardictNormFilt.vcf
#bcftools index -f ${pathBase}/vardictNormFilt.vcf.gz
bgzip ${pathBase}/freeBayesNormFilt.vcf
#bcftools index ${pathBase}/freeBayesNormFilt.vcf.gz

#bgzip ${trueVCF}
#bcftools index ${trueVCF}.gz

#mkdir ${pathBase}/resultLogs

#bcftools stats ${pathBase}/all.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/resultLogs/all.log
#bcftools stats ${pathBase}/fbVard.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/resultLogs/fbVard.log
#bcftools stats ${pathBase}/fbGatk.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/resultLogs/fbGatk.log
#bcftools stats ${pathBase}/vardGatk.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/resultLogs/vardGatk.log
bcftools stats ${pathBase}/gatkFilt.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/resultLogs/gatk.log
bcftools stats ${pathBase}/vardictNormFilt.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/resultLogs/vardict.log
bcftools stats ${pathBase}/freeBayesNormFilt.vcf.gz ${trueVCF}.gz -e 'QUAL<20' > ${pathBase}/resultLogs/fb.log
