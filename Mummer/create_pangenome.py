import glob
import pandas as pd
import os
#import sys_call
import subprocess

#First, read in the reference genome
ref = "Medtr_ref.fasta"

#Next, read in the near-ref quality genomes to align against the reference. Here, they are all stored in a directory
queries = glob.glob("Refererences/*.fasta")

#System call to run various mummer programs -- maybe put in a seperate class 
def sys_call(cmd, verbose=False):
    if verbose:
        print('Running command:', cmd, flush=True)
    try:
        output = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as sys_call_error:
        print(sys_call_error.output, file=sys.stderr)
        raise sys_call_error('Error running command:', cmd)
        
    if verbose:
        print(output)

#Ref and query are the path to fasta files with the ref and query sequences to align
#Prefix is the path to the output file 
def align_nucmer(ref, query, prefix):
  seperator = " "
  nucmer_cmd = seperator.join(["nucmer", "-p", prefix, ref, query]) 
  sys_call(nucmer_cmd, verbose=True) #Set to false when done testing

#run show_diff on nucmer output to generate a report of the breaks in mapping
#input is the output from align_nucmer (prefix_nucmer.delta, where prefix is the align_nucmer argument)
def show_diff(delta_file, output): 
  output = output + '_show-diff.out' #add informative extension to output file
  seperator = " "
  showdiff_cmd = seperator.join(["show-diff", "-q", delta_file, ">", output])
  sys_call(showdiff_cmd, verbose=True) #Set to false when done testing
  
#find GAPs and 
#sd_out is the output from show_diff (output_show-diff.out, where output is the argument to show_diff)
#min_gap is an integer represnting the min length of gaps to include
#returns a df where the 
#TO DO: add functionality to the show_diff function to write output to buffer, then read input from buffer here
def parse_show_diff(sd_out_file, min_gap=500, read_from_buffer = False):
  sd_out = pd.read_csv(sd_out_file, sep = "\t", usecols = ['[SEQ]','[TYPE]','[S1]','[E1]','[LEN 1]'], header=2)
  filter_leng = sd_out[sd_out['[LEN 1]'] > min_gap] #select differences of length 500 or greater
  filter_type = filter_leng[(filter_leng['[TYPE]'] == 'BRK') | (filter_leng['[TYPE]'] == 'GAP')] #selects breaks and gaps
  coords_only = filter_type[['[SEQ]', '[S1]', '[E1]']]
  return coords_only

#seq is a single alignment file
#gap_coords is the data frame with 2 columns - the first is the gap start index and the second is the gap end index
def concat_gaps(seq, gap_coords, out_file, format='fasta'):
  
  full_seq = list(SeqIO.parse(seq, format))
  new_seq = ""

  for index, row in gap_coords.iterrows():
    scaffold = row['[SEQ]']
    scaffold_num = scaffold.replace('scaffold_', '')
    scaffold_num = int(scaffold_num)
    start_gap = row['[S1]']
    end_gap = row['[E1]']
    gap_seq = full_seq[scaffold_num][start_gap:end_gap]
    new_seq = new_seq + gap_seq # this combines the 2 sequences but keeps only the first's id info (do we want/need every scaffold number/id?)

  return(new_seq)

 
  
  
  


  




  
  
  

  


  
