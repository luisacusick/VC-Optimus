# reads-to-variants
Bioinformatics pipeline that uses the results from SNP calls on simulated data to select optimal combination of variant callers (from GATK, Free Bayes, and Vardict) for a given data set. Performs ensemble calling on data set using those callers. 

Required Dependencies:

1. anaconda 
Installation:

git clone https://github.com/luisacusick/reads-to-variants.git
conda env create -f environment.yml 

Workflow:

1. Process Reference (index, create sequence dictionary, and create bed file)
      
      ./processRef.sh -r <reference>
        --this writes the dictionary/indices to the reference's file
  
1. Process Sample (normalize read names and align reads to reference)

      For paired end reads:
      ./processSample.sh -r <reference> -s <read1.fq> -s <read2.fq> -p 

      For unpaired reads:
      ./processSample.sh -r <reference> -s <read.fq>
 
2. Simulate a genome with desired divergence and reads with the same error profile as the real data. Process/align and call variants on simulated reads. Compares results of variant callers to true variants and outputs table with each methods' sensitivity, pvv, and F1.

      ./simulate.sh -r <reference> -d <divergence proportion> -s <read1.fq> -s <read2.fq>
        --default divergence proportion is .01 (means .01 of the sites in the simulated genome will be different from the reference)
        --read error profiles are used to simulate reads with identical error profiles
        --creates a directory called 'sim' in the current working directory 
  
 3. Run the optimal combination of variant callers on your real data set. 
 
      ./callVariants.sh -r <reference> -b <bam file> -o <output directory> -g (include gatk) -f (include free bayes) -v (include vardict) -a (include all)
      --default output directory is 'snpResults.' program will create the directory if it does not already exist.
      --the flags g, f, and v include gatk, free bayes, and vardict in the analysis. -a includes all variant callers
      --the flag -n normalizes (as per GATK vcf rules) and combines variants from the selected methods into one vcf
      

 
 
 
 
  
  
