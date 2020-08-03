# VC-optimus: Variant Calling Optimization Using Simulations
Population genetic analyses typically begin by identifying genetic variants in a sequenced population using any one of several variant calling tools. However, the performance of any single tool is sensitive to the particularities of individual datasets, potentially resulting in the non-optimal recovery of true variants. VC-optimus is  pipeline that enables users to make an *a priori* informed choice about what tool, or combinations of tools, are best suited for analyzing their data of interest. Briefly, VC-optimus simulates a dataset based on user-supplied criteria and compares the precision and recall of three widely-used tools (GATK, Free Bayes, and Vardict) both individually and in combination, on that dataset. The main sources of error simulated by the VC-optimus pipeline are sequencing coverage, platform error, and average nucleotide diversity of the population. Users can then go on to analyze their dataset of interest using VC-optimus. 

Required Dependencies:

1. anaconda 
2. gatk (tested with v4)

Installation:

```bash
git clone https://github.com/luisacusick/reads-to-variants.git
cd reads-to-variants/
conda env create -f environment.yml
```

This will create a conda environment called snvCalling

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
      

 
 
 
 
  
  
