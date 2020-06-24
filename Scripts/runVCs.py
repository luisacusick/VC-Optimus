import pandas as pd
import sys
import subprocess
import Bio
from Bio import SeqIO

#Turn verbose off after testing
def sys_call(cmd, verbose=True):
    if verbose:
        print('Running command:', cmd)
    try:
        output = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as sys_call_error:
        print(sys_call_error.output)
        raise sys_call_error('Error running command:', cmd)
        
    if verbose:
        print(output)

def runGatk(refPath, bamPath, outPath):
  command = 'gatk/gatk HaplotypeCaller -R '+ refPath + ' -I ' + bamPath + ' -O ' + outPath 
  sys_call(command)

def createBed(refPath, outPath):

  bed = pd.DataFrame(columns={'contig', 's', 'e', 'label'})
  with open(refPath) as seq:
    sequences = SeqIO.parse(seq, 'fasta')
    for seq in sequences:
      header = seq.id
      bps = len(seq)

      frags = (bps + 4999) // 5000 #ceiling divison 
      for i in range(1, frags):
        if i == 1:
          start = 1
        else:
          start = ((i-1)*5000)-50
        if(i*5000 > bps):
          end = bps
        else:
          end = i*5000

        label = header + "_" + str(i)
        toAdd = {'contig':header, 's':start, 'e':end, 'label':label}
        bed = bed.append(toAdd, ignore_index=True)
  bed['contig','s','e','label'].to_csv(outPath, sep='\t', index=False, header=False)

def runVarDict(refPath, bamPath, outPath, sampleName, bedPath):
   command = 'vardict-java -U -G ' + refPath + ' -f .05 -N ' + sampleName + ' -b ' +  bamPath + ' -c 1 -S 2 -E 3 -g 4 ' + bedPath + ' | teststrandbias.R | var2vcf_valid.pl -N ' + sampleName + ' -E > ' + outPath
   sys_call(command)

def runFreeBayes(refPath, bamPath, outPath):
  command = 'freebayes -f ' + refPath + ' ' + bamPath + ' > ' + outPath
  sys_call(command)

#To-do: write script to parse user input 
def main():
  bamPath = '/pylon5/eb5phrp/luc32/Data/HM038/aln2_mFlag/aln_M_final.bam'
  refPath = '/pylon5/eb5phrp/luc32/Data/Reference/MedTrunRef.fasta'
  #createBed(refPath, outPath='/pylon5/eb5phrp/luc32/Data/MedTrun.bed')
  bedPath = '/pylon5/eb5phrp/luc32/Data/Reference/MedTrunRef.bed'
  runGatk(refPath, bamPath, outPath = '/pylon5/eb5phrp/luc32/Data/HM038/HM038Results_2/gatk.vcf')
  runVarDict(refPath, bamPath, outPath = '/pylon5/eb5phrp/luc32/Data/HM038/HM038Results_2/vardict.vcf', sampleName = 'HM038', bedPath = '/pylon5/eb5phrp/luc32/Data/toyData/toy.bed')
  runFreeBayes(refPath, bamPath, outPath = '/pylon5/eb5phrp/luc32/Data/HM038/HM038Results_2/freeBayes.vcf')
  

if __name__ == "__main__":
  main()
