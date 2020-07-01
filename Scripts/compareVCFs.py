import argparse
import sys
import subprocess

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

def main():
  angsd = '/pylon5/eb5phrp/luc32/Data/toyData/angsdResults/angsdRes.vcf'
  freeBayes = '/pylon5/eb5phrp/luc32/Data/toyData/freeBayesResults/fbRes.vcf'
  gatk = '/pylon5/eb5phrp/luc32/Data/toyData/gatkResults/testOut.vcf'
  vardict = '/pylon5/eb5phrp/luc32/Data/toyData/vardictResults/vardictToy.vcf'
  combineCommand = 'picard MergeVcfs I=' + angsd + ' I=' + freeBayes + ' I=' + gatk + ' O=/pylon5/eb5phrp/luc32/Data/toyData/combined.vcf'
  sys_call(combineCommand)
  
if __name__ == "__main__":
    main() 


