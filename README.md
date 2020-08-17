# VC-optimus: Variant Calling Optimization Using Simulations
Population genetic analyses typically begin by identifying genetic variants in a sequenced population using one of several variant calling tools. However, the performance of any single tool is sensitive to the particularities of individual datasets (e.g., sequencing coverage and error), potentially biasing the recovery of true variants. VC-optimus is  pipeline that enables users to make an *a priori* informed choice about what tool, or combinations of tools, are best suited for analyzing their data. Briefly, VC-optimus first simulates a read dataset based on the characteristics of user-supplied sequencing read data and a reference genome. It then compares the precision and recall of three widely-used programs (gatk, Free Bayes, and Vardict) both individually and in combination on the simulated dataset, such that users are able to identify the optimal variant caller or combination of variant callers to analyze their data with. Users can then continue using VC-optimus to analyze their real data with their variant callers of choice.

VC-optimus makes use of the following tools: <INSERT_LIST_WITH_REFS>

**Dependencies requiring user installation:**

1. anaconda 
2. gatk (tested with v4)

**Conda installation for all other dependencies:**

```bash
git clone https://github.com/luisacusick/reads-to-variants.git
cd reads-to-variants/
chmod 755 Scripts/*
conda env create -f environment.yml
source activate snvCalling
```

**Workflow:**

1. Process reference
```bash
./processRef.sh -r <reference.fa>
``` 
   * this writes the dictionary/indices and bed file to the reference file directory
  
2. Process sample

   For paired end reads:
```bash
./processSample.sh -r <reference.fa> -s <read1.fq> -s <read2.fq> -p
```
   For unpaired reads:
```bash
./processSample.sh -r <reference.fa> -s <read.fq>
```
   * this normalizes read names and align reads to reference
   
3. Simulate reads, call variants, and compare called variants to known variants

```bash
./simulate.sh -r <reference> -d <divergence proportion> -s <read1.fq> -s <read2.fq>
```
   * This simulates a genome with desired sequence divergence from reference (known variants) and simulates sequence reads from that genome using error profiles sampled from the real read dataset. Simulated reads are processed and aligned back to the reference, and subsequently analyzed using the variant calling software gatk, Free Bayes, and VarDict (called variants). Outputs a table and figure describing each methods' sensitivity, positive predictive value (PPV), and F1 score <DESCRIBE_F1> either singly or in combination
   * default divergence proportion is .01 (means .01 of the sites in the simulated genome will be different from the reference)
   * creates a directory called 'sim' in the current working directory 
  
4. Run the optimal combination of variant callers on your real data set

```bash
./callVariants.sh -r <reference> -b <bam file> -o <output directory> -g (include gatk) -f (include free bayes) -v (include vardict) -a (include all)
```
   * default output directory is 'snpResults.' program will create the directory if it does not already exist.
   * the flags g, f, and v include gatk, free bayes, and vardict in the analysis. -a includes all variant callers
   * the flag -n normalizes (as per GATK vcf rules) and combines variants from the selected methods into one vcf
      

 
 
 
 
  
  
