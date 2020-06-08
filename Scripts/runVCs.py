import sys

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
  command = 'gatk HaplotypeCaller -R +' refPath + ' -I' + bamPath + ' -O' + outPath ' -ERC GVCF'
  sys_call(command)

def runAngsd(bamPath, outPath):
  command = 'angsd -GL 2 -i' + bamPath +  ' -out' + outPath + ' -doBcf 1 -doMajorMinor 1 -doPost 1 -doMaf 1'
  sys_call(command)

def runVarDict(refPath, bamPath, outPath, sampleName, bedPath):
   command = ' vardict -G ' + outPath + ' -f .05 -N ' + HM106 + ' -b ' +  bamPath + ' -d ' ' -c 1 -S 2 -E 3 -g 4 ' + bedPath + ' | teststrandbias.R | var2vcf_valid.pl -N ' + sampleName + ' -E > ' + outPath
   sys_call(command)

def runFreeBayes(refPath, bamPath, outPath):
  command = 'freebayes -f ' + refPath + ' ' + bamPath ' > ' + outPath
  sys_call(command)

#To-do: write script to parse user input 
def main():
  bamPath = '/pylon5/eb5phrp/luc32/Data/toyData/aln/toyAlnRGSort.bam'
  refPath = '/pylon5/eb5phrp/luc32/Data/toyData/aln/toyMedTrunRef.fasta'

  runGatk(refPath, bamPath, outPath = '/pylon5/eb5phrp/luc32/Data/toyData/gatkResults')
  runAngsd(bamPath, outPath = '/pylon5/eb5phrp/luc32/Data/toyData/angsdResults')
  runVarDict(refPath, bamPath, outPath = '/pylon5/eb5phrp/luc32/Data/toyData/vardictResults')
  runFreeBayes(refPath, bamPath, outPath = '/pylon5/eb5phrp/luc32/Data/toyData/freeBayesResults')

if __name__ == "__main__":
  main()
