import pandas as pd
import sys
import subprocess
import Bio
from Bio import SeqIO

def createBed(refPath, outPath):

  bed = pd.DataFrame(columns={'1', '2', '3', '4'})
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
        toAdd = {'1':header, '2':start, '3':end, '4':label}
        bed = bed.append(toAdd, ignore_index=True)
  bed = bed[['1', '2', '3', '4']]
  bed.to_csv(outPath, sep='\t', index=False, header=False)
 
def main():
  args = sys.argv
  ref = args[1]
  out = args[2]
  createBed(refPath = ref, outPath = out)
  
if __name__ == "__main__":
  main()
