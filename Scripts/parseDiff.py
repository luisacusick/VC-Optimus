import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import glob 
import os
import argparse
import sys

#plt.switch_backend('agg') #necessary for graphic creation if working on cluster 

def computeSens(tp, fn):
  sens= float(tp)/(float(tp)+float(fn))
  return(sens)

def computePpv(tp, fp):
  return (float(tp)/(float(tp)+float(fp)))

def computeF1(sens, ppv):
  return (2*float(sens)*float(ppv)/(float(sens)+float(ppv)))

def extractInt(line):
  num = [int(s) for s in line.split() if s.isdigit()]
  return num[1]

#ppv, sens, and pointNames are vectors, title is a string
def plotDiff(ppv, sens, pointNames, figName='/pylon5/eb5phrp/luc32/methods.png'):
  print('Called plotDiff')
  assert len(ppv) == len(sens), 'ppv and sensitivity are not the same length'
  assert len(ppv) == len(labels), 'coord len is not equal to labels'
  
  for i, label in enumerate(pointNames):
    plt.scatter(ppv[i], sens[i])
    plt.annotate(label, xy=(ppv[i], sens[i]), textcoords='offset points', ha='right', va= 'bottom')
  
  f=plt.figure()
  f.savefig(figName, bbox_inches='tight', format="png")
  plt.close()

def parseLog(diffFile):
  with open(diffFile) as f: #read in the log file (contains numbers of tps, fps, and fns)
    diff = f.readlines()
  
  diff = [x.strip() for x in diff] #remove whitespace characters at the end of each line

  lineTp = diff[44] #line that contains tp number
  tp = extractInt(lineTp)
  
  lineFp = diff[28]            
  fp = extractInt(lineFp)
  
  lineFn = diff[36]
  fn = extractInt(lineFn)  
  
  sens = computeSens(tp, fn)
  ppv = computePpv(tp, fp)
  f1 = computeF1(sens, ppv)

  return[sens, ppv, f1]

#first arg is the path to the output folder containing logs, second arg is the path to desired output file with run stats
def main():
  
  args = sys.argv
  pathPrefix = args[1]
  output = args[2]
  
  methods = ''
  logFiles = glob.glob(pathPrefix + '/*.log')
   
  for file in logFiles:
    print(file)
    filename, file_extension = os.path.splitext(os.path.basename(file))
    row = parseLog(file)
    methods = methods + filename + ' ' + str(row[0]) + ' ' + str(row[1]) + ' ' + str(row[2]) + '\n'
  
  outF = open(output, 'w')
  outF.writelines('Method' + ' 	Sens' + ' PPV' + ' F1'+ '\n')
  outF.writelines(methods)
  outF.close()
   

if __name__ == "__main__":
  main()
