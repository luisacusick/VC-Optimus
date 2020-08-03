# VC-optimus: Variant Calling Optimization Using Simulations
Population genetic analyses typically begin by identifying genetic variants in a sequenced population using any one of several variant calling tools. However, the performance of any single tool is sensitive to the particularities of individual datasets, potentially resulting in the non-optimal recovery of true variants. VC-optimus is  pipeline that enables users to make an *a priori* informed choice about what tool, or combinations of tools, are best suited for analyzing their data of interest. Briefly, VC-optimus simulates a dataset based on user-supplied criteria and compares the precision and recall of three widely-used tools (GATK, Free Bayes, and Vardict) both individually and in combination, on that dataset. The main sources of error simulated by the VC-optimus pipeline are sequencing coverage, sequencing error, and expected sequence divergence within population. Users can then go on to analyze their dataset of interest using the optimal variant caller or combination of variant callers identified based on the results of the simulation.

VC-optimus makes use of the following tools: <INSERT_LIST_WITH_REFS>

**Dependencies requiring user installation:**

1. anaconda 
2. gatk (tested with v4)

**Conda installation for all other dependencies:**

```bash
git clone https://github.com/luisacusick/reads-to-variants.git
cd reads-to-variants/
conda env create -f environment.yml
source activate snvCalling
```

**Workflow:**

1. Process reference
```bash
./processRef.sh -r <reference.fa>
``` 
   --this writes the dictionary/indices and bed file to the reference file directory
  
2. Process sample

   For paired end reads:
```bash
./processSample.sh -r <reference.fa> -s <read1.fq> -s <read2.fq> -p
```
   For unpaired reads:
```bash
./processSample.sh -r <reference.fa> -s <read.fq>
```
   --this normalizes read names and align reads to reference
   
3. Simulate reads, call variants, and compare called variants to known variants

```bash
./simulate.sh -r <reference> -d <divergence proportion> -s <read1.fq> -s <read2.fq>
```
   --This simulates a genome with desired sequence divergence from reference (known variants) and simulates sequence reads from that genome using error profiles sampled from the real read dataset. Simulated reads are processed and aligned back to the reference, and subsequently analyzed using the variant calling software GATK, Free Bayes, and VarDict (called variants). Outputs a table and figure describing each methods' sensitivity, positive predictive value (PPV), and F1 score <DESCRIBE_F1> either singly or in combination
   --default divergence proportion is .01 (means .01 of the sites in the simulated genome will be different from the reference)
   --creates a directory called 'sim' in the current working directory 
  
4. Run the optimal combination of variant callers on your real data set

```bash
./callVariants.sh -r <reference> -b <bam file> -o <output directory> -g (include gatk) -f (include free bayes) -v (include vardict) -a (include all)
```
   --default output directory is 'snpResults.' program will create the directory if it does not already exist.
   --the flags g, f, and v include gatk, free bayes, and vardict in the analysis. -a includes all variant callers
   --the flag -n normalizes (as per GATK vcf rules) and combines variants from the selected methods into one vcf
      

 
 
 
 
  
  
