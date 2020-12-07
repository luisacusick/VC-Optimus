#!/bin/bash

PARAM=''
#All user input in param file
while getopts ":hp:" option; do
  case ${option} in
  h) echo 'Usage: simulate_wrapper.sh -p <param_file>'
     echo 'See README.md and sim_param.txt for details about the parameter file'
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

#If user doesn't include an output directory in the param file then set directory to vc-optimus-output/
if [ ! ${PARAM_ARRAY[OUTPUT_DIRECTORY]+_} ] 
then
  PARAM_ARRAY+=([OUTPUT_DIRECTORY]=vc-optimus-output/)
fi

#If output directory with tmp directory doesn't exist then print a message about running preprocessing steps and exit
if [ ! -d ${PARAM_ARRAY[OUTPUT_DIRECTORY]}/tmp ]
then
  echo "VC-Optimus tmp directory not found. Please execute the pipeline's processing step before continuing. Note that the output directory for each step must be the same."
  exit 0
fi

if [ ! -d ${PARAM_ARRAY[OUTPUT_DIRECTORY]}/simulations ]
then
  mkdir ${PARAM_ARRAY[OUTPUT_DIRECTORY]}/simulations ] 
  echo "Created directory ${PARAM_ARRAY[OUTPUT_DIRECTORY]}/simulations"
fi

simDir=${PARAM_ARRAY[OUTPUT_DIRECTORY]}/simulations/sim.$(date "+%Y.%m.%d-%H.%M.%S")
mkdir ${simDir}

#Call simulate.sh for paried reads
if [ ${PARAM_ARRAY[PAIRED]} = true ]
then 
  ${scriptDir}/simulate.sh -r ${PARAM_ARRAY[REFERENCE]} -d ${PARAM_ARRAY[DIVERGENCE]} -s ${PARAM_ARRAY[SAMPLE1]} -t ${PARAM_ARRAY[SAMPLE2]} -o ${simDir} 1> ${simDir}/simulate.out 2> ${simDir}/simulate.err 
else 
  ${scriptDir}/simulate.sh -r ${PARAM_ARRAY[REFERENCE]} -d ${PARAM_ARRAY[DIVERGENCE]} -s ${PARAM_ARRAY[SAMPLE1]} -o ${simDir} 1> ${simDir}/simulate.out 2> ${simDir}/simulate.err
fi
