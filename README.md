# reads-to-variants
Bioinformatics pipeline that uses the results from SNP calls on simulated data to select optimal combination of variant callers (from GATK, Free Bayes, and Vardict) for a given data set. Performs ensemble calling on data set using those callers. 

Installation:

Workflow:

1. Process Reference (index, create sequence dictionary, and create bed file)
      
      ./processRef.sh -r <reference>
        --this writes the dictionary/indices to the reference's file
  
1. Process Sample (normalize read names and align reads to reference)

      For paired end reads:
      ./processSample.sh -r <reference> -s <read1.fq> -s <read2.fq> -p 

      For unpaired reads:
      ./processSample.sh -r <reference> -s <read.fq>
 
2. Simulate a genome with desired divergence and reads with the same error profile as the real data. Process/align and call variants on simulated reads.

      ./simulate.sh -r <reference> -d <divergence proportion> -s <read1.fq> -s <read2.fq>
        --default divergence proportion is .01 (means .01 of the sites in the simulated genome will be different from the reference)
        --read error profiles are used to simulate reads with identical error profiles
        --creates a directory called 'sim' in the current working directory 
  
 
  
  
