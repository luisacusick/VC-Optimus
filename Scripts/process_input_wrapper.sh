#!/bin/bash

PARAM=''
#All user input in param file
while getopts ":hp:" option; do
  case ${option} in
  h) echo 'Usage: process_ref.sh -p <param_file>'
     echo 'See README.md and ref_param.txt for details about the parameter file'
     exit 0;;
  p) PARAM=${OPTARG};;
  esac
done
shift "$((OPTIND -1))"

if [[ -f ${PARAM} ]]
then
  echo "Found ${PARAM}"
else
  echo "Parameter file not found."
  exit 0
fi

#all VC-optimus scripts live in the same directory; save the path to call the others
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source activate VC-optimus
#Add the parameters in the config file to an associative array
declare -A PARAM_ARRAY
n=0
while IFS=$'[ \t]*=[ \t]*' read -r name value
do
  if [ -z "${name}" ]; then #skip blank lines
    continue
  fi
  echo "Read $name = $value" >&2 #set name and value pair based on param file
  PARAM_ARRAY[$name]="$value"
done < ${PARAM}

#To-Do: decide on input verification strategy
#verify that -p (paired flag) is true or false

#If user doesn't include an output directory in the param file then set output to sim/sim.timestamp
if [ ! ${PARAM_ARRAY[OUTPUT_DIRECTORY]+_} ] 
then
  PARAM_ARRAY+=([OUTPUT_DIRECTORY]=vc-optimus-output/)
fi

#If output directory doesn't exist then create it
if [ ! -d ${PARAM_ARRAY[OUTPUT_DIRECTORY]} ]
then
  mkdir ${PARAM_ARRAY[OUTPUT_DIRECTORY]}
fi

#Call processSample for paired reads
if [ ${PARAM_ARRAY[PAIRED]} = true ]
then 
  ${scriptDir}/processSample.sh -r ${PARAM_ARRAY[REFERENCE]} -s ${PARAM_ARRAY[READ1]} -s ${PARAM_ARRAY[READ2]} -o ${PARAM_ARRAY[OUTPUT_DIRECTORY]}
else #call scripts for unpaired reads/single sample
  ${scriptDir}/processSample.sh -r ${PARAM_ARRAY[REFERENCE]} -s ${PARAM_ARRAY[READ1]} -o ${PARAM_ARRAY[OUTPUT_DIRECTORY]} 
fi
