#!/bin/bash

PARAM=''
#All user input in param file
while getopts ":hp:" option; do
  case ${option} in
  h) echo 'Usage: reads-to-variants.sh -p <param_file>'
     echo 'See README.md and sample_param.txt for details about the parameter file'
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
#verify that -p is true or false

#Call processRef
./processRef.sh -r ${PARAM_ARRAY[REFERENCE]}

#Call processSample and simulate.sh for paried reads
if [ ${PARAM_ARRAY[PARIED]} = true ]
then 
  ./processSample.sh -r ${PARAM_ARRAY[REFERENCE]} -s ${PARAM_ARRAY[SAMPLE1]} -s ${PARAM_ARRAY[SAMPLE2]}
  ./simulate.sh -r ${PARAM_ARRAY[REFERENCE]} -d ${PARAM_ARRAY[DIVERGENCE]} -s ${PARAM_ARRAY[SAMPLE1]} -s ${PARAM_ARRAY[SAMPLE2]}
else #call scripts for unpaired reads/single sample
  ./processSample.sh -r ${PARAM_ARRAY[REFERENCE]} -s ${PARAM_ARRAY[SAMPLE1]}
  ./simulate.sh -r ${PARAM_ARRAY[REFERENCE]} -d ${PARAM_ARRAY[DIVERGENCE]} -s ${PARAM_ARRAY[SAMPLE1]}
fi

#Call variants real

#Call norm and combine

#Call parsediff
