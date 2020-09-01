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

#Call processRef and processSim

#Call simulate script

#Call variants real

#Call norm and combine

#Call parsediff
