# VC-optimus: Variant Calling Optimization Using Simulations
Population genetic analyses typically begin by identifying genetic variants in a sequenced population using one of several variant calling tools. However, the performance of any single tool is sensitive to the particularities of individual datasets (e.g., sequencing coverage and error), potentially biasing the recovery of true variants. VC-optimus is pipeline that enables users to make an *a priori* informed choice about what tool, or combinations of tools, are best suited for analyzing their data. Briefly, VC-optimus first simulates a read dataset based on the characteristics of user-supplied sequencing read data and a reference genome. It then compares the precision and recall of three widely-used programs (gatk, Free Bayes, and Vardict) both individually and in combination on the simulated dataset, such that users are able to identify the optimal variant caller or combination of variant callers to analyze their data with. Users can then continue using VC-optimus to analyze their real data with their variant callers of choice.

**Dependencies requiring user installation:**

1. [anaconda](https://docs.anaconda.com/anaconda/install/)<sup>1</sup>
2. [gatk](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)<sup>2</sup> (tested with v4)
3. [vt](https://genome.sph.umich.edu/wiki/Vt)<sup>3</sup> (tested with 0.5772)

**Dependencies prepackaged in the /Scripts directory:**

4. simuG<sup>4</sup> (slightly modified from original for execution purposes)
5. Artificial Fastq Generator<sup>5</sup>

**Conda installation for all remaining dependencies<sup>6</sup><sup>,7</sup><sup>,8</sup><sup>,9</sup><sup>,10</sup><sup>,11</sup>:**

```bash
git clone https://github.com/luisacusick/reads-to-variants.git
cd reads-to-variants/
chmod 755 Scripts/*
conda env create -f environment.yml
```

**Output:**

This is a sample directory structure for a vc-optimus run with 3 simulations and all of the default folder names. Note that all of folder names are parameters with the exception of 'tmp,' and the program creates any that do not exist. If the output folder already exists, the program could over-write some of its content.

```bash
vc-optimus-output/ <--Top-level output directory name, set to vc-optimus-output by default
├── simulatedSamples <-- Folder with all simulations
│   ├── sim1 <-- Outut from a single simulation, name specified by user
|   |   |──processedSamples <-- Output from processing  simulated samples, contains bam alignment and alignment index
|   |   |──vcfs <-- Output from variant calling on simulated sample along with log files from comparison to true VCF
│   ├── sim2
│   ├── sim3
├── processedSamples <-- Output from processing real samples, contains bam alignment and alignment index
└── vcfs <-- Output from variant calling on real sample 
```

**Before starting:**

   * the reference genome should ideally be named <genome code>.fasta
   * input fastq files should already be quality checked and trimmed, and unzipped 
   * <insert advice about how to combine read libraries here for error profiling, is it a simple cat?>
   * insert the path to the gatk<sup>2</sup> executable where indicated at the top of the runVCs.sh script
   * insert the path to the vt<sup>3</sup> executable where indicated at the top of the normAndCombineVcfs.sh script
   * insert the path to the simuG.pl<sup>4</sup> executable where indicated at the top of the simulateSample.sh script (should come packaged in the /Scripts directory)
   * insert the path to the ArtificialFastqGenerator.jar<sup>5</sup> executable where indicated at the top of the simulateSample.sh script (should come packaged in the /Scripts directory)


**Workflow:**

   * use -h switch for any script to list all available options

1. Initiate environment and Perl libraries (must be done at the start of each session)
```
source activate VC-optimus
export PERL5LIB=~/.conda/envs/VC-optimus/lib/5.26.2/
```

2. Process reference
```bash
./processRef.sh -r ref.fa
``` 
   * Writes the dictionary/indices and bed file to the reference file directory
  
3. Process sample

   For paired end reads:
```bash
./processSample.sh -n name -r ref.fa -s read1.fq -s read2.fq [OPTIONS]
```
   * Normalizes read names, aligns reads to reference, and pre-processes data for variant calling (sorting BAM, adding read group names, marking duplicates, and indexing)
   
4. Simulate known variants, call variants using simulated reads, and compare called simulated variants to known simulated variants

```bash
./simulateSample.sh -n name -r ref.fa -s read1.fq -s read2.fq [OPTIONS]
```
   * Simulates a genome with desired sequence divergence from reference (thus generating known variants), and then simulates sequence reads from that genome using error profiles sampled from the real read dataset. Simulated reads are processed and aligned back to the reference, and subsequently analyzed using the variant calling software gatk, Free Bayes, and VarDict (called variants). Outputs a table and figure describing each methods' sensitivity, positive predictive value (PPV), and F1 score (the avg. of sensitivity and ppv, calculated as 2*(sensitivity * ppv)/(sensitivity + ppv) either singly or in combination
   * default divergence proportion is .01 (means 1% of the nucleotide sites in the simulated genome will be different from the reference)
   * creates a directory called 'simulatedSamples' in vc-optimus-output/, or in user-specified output directory 
  
5. Run the optimal combination of variant callers on your real data set

```bash
./runVCs.sh -r <reference> -b <bam file> -o <output directory> -g (include gatk) -f (include free bayes) -v (include vardict) [OPTIONS]
```
   * default output directory is 'snpResults.' program will create the directory if it does not already exist.
   * the flags g, f, and v include gatk, free bayes, and vardict in the analysis.
   * the flag -n normalizes (as per GATK vcf rules) and combines variants from the selected methods into one vcf
      
**References**

1. https://docs.anaconda.com/
2. Van der Auwera, Geraldine A., et al. "From FastQ data to high‐confidence variant calls: the genome analysis toolkit best practices pipeline." Current protocols in bioinformatics 43.1 (2013): 11-10.
3. Tan, Adrian, Gonçalo R. Abecasis, and Hyun Min Kang. "Unified representation of genetic variants." Bioinformatics 31.13 (2015): 2202-2204.
4. Yue, Jia-Xing, and Gianni Liti. "simuG: a general-purpose genome simulator." Bioinformatics 35.21 (2019): 4442-4444.
5. Frampton, Matthew, and Richard Houlston. "Generation of artificial FASTQ files to evaluate the performance of next-generation sequencing pipelines." PLoS One 7.11 (2012): e49110.
6. Broad Institute. “Picard Toolkit.” GitHub Repository (2019). http://broadinstitute.github.io/picard/
7. https://github.com/biopet/vcffilter
8. Li, Heng, et al. "The sequence alignment/map format and SAMtools." Bioinformatics 25.16 (2009): 2078-2079.
9. Li, Heng, and Richard Durbin. "Fast and accurate short read alignment with Burrows–Wheeler transform." bioinformatics 25.14 (2009): 1754-1760.
10. Lai, Zhongwu, et al. "VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research." Nucleic acids research 44.11 (2016): e108-e108.
11. Garrison, Erik, and Gabor Marth. "Haplotype-based variant detection from short-read sequencing." arXiv preprint arXiv:1207.3907 (2012).
  
