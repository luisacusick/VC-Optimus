#!/bin/bash

module load fastqc

declare -a accessions
n=0
while read line
do
    accessions[n++]=${line}
done < "/pylon5/eb5phrp/luc32/missingQR.txt"

counter=1

while [ ${counter} -lt 1157 ]
do
  sra=${accessions[$counter]}
  echo $sra
  fastqc "/pylon5/eb5phrp/luc32/Data/medTrunSRAs/${sra}" --outdir "/pylon5/eb5phrp/luc32/Data/fastqcResults" --noextract 
  let counter=counter+1
done
